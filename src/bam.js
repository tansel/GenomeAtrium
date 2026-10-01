/*
 * Streaming BAM parser. Reads the whole file once (no index needed) and
 * keeps summaries, not reads:
 *   - depth track (mean depth per bin, like samtools depth: skips unmapped,
 *     secondary, QC-fail and duplicate reads)
 *   - discordant read pairs (other contig, or insert far above the median)
 *   - splice junctions (N operations in the CIGAR)
 *   - read length and MAPQ histograms
 * Unaligned BAMs (no @SQ lines, e.g. dorado basecaller output) only get
 * the read length and base quality summaries.
 *
 * BAM positions are 0-based; everything stored is converted to 1-based.
 */
(function (G) {
  var F = { PAIRED: 1, PROPER: 2, UNMAP: 4, MUNMAP: 8, REVERSE: 16, READ1: 64,
            SECONDARY: 256, QCFAIL: 512, DUP: 1024, SUPPLEMENTARY: 2048 };
  var MIN_DISCORDANT = 1000;  // floor for "insert too long" regardless of library
  var PAIR_CLUSTER = 5000;    // discordant pairs are grouped into 5 kb x 5 kb cells
  var MAX_ARC_KEYS = 300000;
  var LEN_BINS = 64;          // read length histogram, log2 scale

  function lenBin(len) { return Math.min(LEN_BINS - 1, Math.floor(Math.log2(Math.max(1, len)) * 3)); }

  // Reads the CG:B,I tag that holds the real CIGAR when it has > 65535 ops.
  function cigarFromCG(view, bytes, p, end) {
    while (p + 3 <= end) {
      var t0 = bytes[p], t1 = bytes[p + 1], ty = String.fromCharCode(bytes[p + 2]);
      p += 3;
      if (ty === 'B') {
        var sub = String.fromCharCode(bytes[p]), n = view.getInt32(p + 1, true);
        p += 5;
        var size = { c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 }[sub];
        if (t0 === 67 && t1 === 71 && sub === 'I') {
          var ops = new Uint32Array(n);
          for (var k = 0; k < n; k++) ops[k] = view.getUint32(p + 4 * k, true);
          return ops;
        }
        p += n * size;
      } else if (ty === 'Z' || ty === 'H') {
        while (p < end && bytes[p] !== 0) p++;
        p++;
      } else {
        p += { A: 1, c: 1, C: 1, s: 2, S: 2, i: 4, I: 4, f: 4 }[ty] || 0;
      }
    }
    return null;
  }

  async function parse(blob, opts) {
    opts = opts || {};
    var gm = G.genome;
    var rd = new G.bgzf.ByteReader(blob, opts);

    if (!(await rd.need(8))) throw new Error('File too short to be BAM');
    var magic = rd.take(4);
    if (magic[0] !== 66 || magic[1] !== 65 || magic[2] !== 77 || magic[3] !== 1) throw new Error('Not a BAM file (bad magic)');
    var lText = new DataView(rd.take(4).slice().buffer).getInt32(0, true);
    await rd.need(lText + 4);
    var headerText = new TextDecoder().decode(rd.take(lText)).replace(/\0+$/, '');
    var nRef = new DataView(rd.take(4).slice().buffer).getInt32(0, true);
    var refs = [];
    for (var r = 0; r < nRef; r++) {
      await rd.need(4);
      var lName = new DataView(rd.take(4).slice().buffer).getInt32(0, true);
      await rd.need(lName + 4);
      var nm = new TextDecoder().decode(rd.take(lName)).replace(/\0+$/, '');
      var lRef = new DataView(rd.take(4).slice().buffer).getInt32(0, true);
      refs.push({ name: nm, length: lRef });
    }

    var genome = new gm.Genome(refs);
    var binSize = gm.chooseBinSize(refs.reduce(function (s, c) { return s + c.length; }, 0) || 1);
    var depth = refs.map(function (c) { return new gm.Track(c.length, binSize, 'mean'); });
    var readsOn = new Float64Array(nRef);

    var stats = { reads: 0, mapped: 0, unmapped: 0, secondary: 0, supplementary: 0, dup: 0, qcfail: 0,
                  paired: 0, properPairs: 0, discordantInter: 0, discordantLong: 0, spliced: 0,
                  cgCigars: 0, bases: 0, aligned: nRef > 0 };
    var lenHist = new Float64Array(LEN_BINS), mapqHist = new Float64Array(61), qualHist = new Float64Array(61);
    var pairCells = new Map(), junctions = new Map(), longCandidates = [];
    var tlenSample = [], tlenSeen = 0;
    var maxReads = opts.maxReads || Infinity;

    while (stats.reads < maxReads) {
      if (opts.signal && opts.signal.aborted) break;
      if (!(await rd.need(4))) break;
      var bs = new DataView(rd.buf.buffer, rd.buf.byteOffset + rd.off, 4).getInt32(0, true);
      if (!(await rd.need(4 + bs))) break;
      rd.take(4);
      var rec = rd.take(bs);
      var v = new DataView(rec.buffer, rec.byteOffset, rec.byteLength);

      var refID = v.getInt32(0, true), pos0 = v.getInt32(4, true);
      var lReadName = rec[8], mapq = rec[9], nCigar = v.getUint16(12, true), flag = v.getUint16(14, true);
      var lSeq = v.getInt32(16, true), nextRef = v.getInt32(20, true), nextPos0 = v.getInt32(24, true), tlen = v.getInt32(28, true);
      var cigarOff = 32 + lReadName;
      stats.reads++;
      stats.bases += lSeq;

      if (flag & F.SECONDARY) { stats.secondary++; continue; }
      if (flag & F.SUPPLEMENTARY) stats.supplementary++;
      else lenHist[lenBin(lSeq)]++;

      if (!stats.aligned || refID < 0 || (flag & F.UNMAP)) {
        stats.unmapped++;
        if (!(flag & F.SUPPLEMENTARY) && lSeq > 0) { // mean base quality of the read
          var qOff = cigarOff + 4 * nCigar + ((lSeq + 1) >> 1), qs = 0;
          if (rec[qOff] !== 255) {
            for (var q = 0; q < lSeq; q++) qs += rec[qOff + q];
            qualHist[Math.min(60, Math.round(qs / lSeq))]++;
          }
        }
        continue;
      }
      stats.mapped++;
      mapqHist[Math.min(60, mapq)]++;
      if (flag & F.QCFAIL) { stats.qcfail++; continue; }
      if (flag & F.DUP) { stats.dup++; continue; }
      readsOn[refID]++;

      var ops = null;
      if (nCigar === 2) {
        var op0 = v.getUint32(cigarOff, true);
        if ((op0 & 15) === 4 && (op0 >>> 4) === lSeq && (v.getUint32(cigarOff + 4, true) & 15) === 3) {
          var tagOff = cigarOff + 8 + ((lSeq + 1) >> 1) + lSeq;
          ops = cigarFromCG(v, rec, tagOff, rec.length);
          if (ops) stats.cgCigars++;
        }
      }
      if (!ops) { ops = new Uint32Array(nCigar); for (var c = 0; c < nCigar; c++) ops[c] = v.getUint32(cigarOff + 4 * c, true); }

      // Walk the CIGAR: M/=/X add depth, N records a junction, D just advances.
      var refPos = pos0 + 1, spliced = false, track = depth[refID];
      for (var i = 0; i < ops.length; i++) {
        var opc = ops[i] & 15, len = ops[i] >>> 4;
        if (opc === 0 || opc === 7 || opc === 8) { track.addSpan(refPos, refPos + len - 1, 1); refPos += len; }
        else if (opc === 2) refPos += len;
        else if (opc === 3) {
          var jk = refID + ':' + (refPos - 1) + ':' + (refPos + len);
          junctions.set(jk, (junctions.get(jk) || 0) + 1);
          refPos += len; spliced = true;
        }
      }
      if (spliced) stats.spliced++;

      if ((flag & F.PAIRED) && !(flag & F.MUNMAP) && !(flag & F.SUPPLEMENTARY) && nextRef >= 0) {
        stats.paired++;
        if (flag & F.PROPER) {
          stats.properPairs++;
          if (tlen > 0) { // reservoir sample of insert sizes
            tlenSeen++;
            if (tlenSample.length < 20000) tlenSample.push(tlen);
            else { var j = Math.floor(Math.random() * tlenSeen); if (j < 20000) tlenSample[j] = tlen; }
          }
        }
        // Count each pair once: from the mate with the smaller coordinate.
        var first = refID < nextRef || (refID === nextRef && (pos0 < nextPos0 || (pos0 === nextPos0 && (flag & F.READ1))));
        if (first) {
          if (nextRef !== refID) {
            stats.discordantInter++;
            addPair(refID, pos0 + 1, nextRef, nextPos0 + 1);
          } else if (Math.abs(tlen) >= MIN_DISCORDANT && !(flag & F.PROPER)) {
            longCandidates.push(refID, pos0 + 1, nextPos0 + 1, Math.abs(tlen));
          }
        }
      }
      if (opts.onRecord) opts.onRecord(stats);
    }

    function addPair(r0, p0, r1, p1) {
      var key = r0 + ':' + Math.floor(p0 / PAIR_CLUSTER) + ':' + r1 + ':' + Math.floor(p1 / PAIR_CLUSTER);
      var cell = pairCells.get(key);
      if (cell) { cell.n++; cell.s0 += p0; cell.s1 += p1; return; }
      if (pairCells.size >= MAX_ARC_KEYS) return;
      pairCells.set(key, { r0: r0, r1: r1, n: 1, s0: p0, s1: p1 });
    }

    // Insert size cut-off: median + 6 * MAD of proper pairs, never below 1 kb.
    tlenSample.sort(function (a, b) { return a - b; });
    var median = tlenSample.length ? tlenSample[tlenSample.length >> 1] : 0;
    var devs = tlenSample.map(function (t) { return Math.abs(t - median); }).sort(function (a, b) { return a - b; });
    var mad = devs.length ? devs[devs.length >> 1] : 0;
    var cutoff = Math.max(MIN_DISCORDANT, median + 6 * mad);
    for (var k = 0; k < longCandidates.length; k += 4) {
      if (longCandidates[k + 3] < cutoff) continue;
      stats.discordantLong++;
      addPair(longCandidates[k], longCandidates[k + 1], longCandidates[k], longCandidates[k + 2]);
    }

    genome.finish(function (c) { return readsOn[refs.findIndex(function (x) { return x.name === c.name; })] > 0; });

    var tracks = {};
    refs.forEach(function (c, i) { depth[i].buildPyramid(); tracks[gm.normName(c.name)] = { depth: depth[i] }; });

    var arcs = [];
    pairCells.forEach(function (cell) {
      var p0 = Math.round(cell.s0 / cell.n), p1 = Math.round(cell.s1 / cell.n);
      var inter = cell.r0 !== cell.r1;
      arcs.push({ c0: refs[cell.r0].name, p0: p0, c1: refs[cell.r1].name, p1: p1,
        type: inter ? 'pair_inter' : 'pair_long', support: cell.n, pass: 1,
        label: cell.n + ' discordant pair' + (cell.n > 1 ? 's' : '') + ': ' + refs[cell.r0].name + ':' + p0.toLocaleString() +
               ' to ' + refs[cell.r1].name + ':' + p1.toLocaleString() + (inter ? '' : ' (' + (p1 - p0).toLocaleString() + ' bp)') });
    });
    junctions.forEach(function (n, key) {
      var parts = key.split(':'), ri = +parts[0], s = +parts[1], e = +parts[2];
      arcs.push({ c0: refs[ri].name, p0: s, c1: refs[ri].name, p1: e, type: 'junction', support: n, pass: 1,
        label: 'splice junction ' + refs[ri].name + ':' + s.toLocaleString() + '-' + e.toLocaleString() +
               ' (' + (e - s - 1).toLocaleString() + ' bp intron, ' + n + ' read' + (n > 1 ? 's' : '') + ')' });
    });
    arcs.sort(function (a, b) { return b.support - a.support; });

    return {
      format: 'bam', genome: genome, build: genome.build, binSize: binSize, header: headerText,
      stats: stats, tracks: tracks, arcs: arcs,
      insert: { median: median, mad: mad, cutoff: cutoff, sampled: tlenSample.length },
      lenHist: lenHist, mapqHist: mapqHist, qualHist: qualHist, lenBinToLength: function (b) { return Math.pow(2, b / 3); },
      aborted: !!(opts.signal && opts.signal.aborted) || stats.reads >= maxReads
    };
  }

  G.bam = { parse: parse, FLAGS: F };
})(globalThis.G = globalThis.G || {});
