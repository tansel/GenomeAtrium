/*
 * Reads for a region, from a BAM file with its .bai index: a local file (read
 * with File.slice, nothing leaves the page) or a URL that serves byte ranges
 * (the public GIAB HG002 BAM in the demo). Only the compressed blocks the
 * index points to for the region are read and inflated.
 *
 *   var src = await G.reads.open(bamSource, baiSource)   // sources: G.reads.fileSource(file) or urlSource(url)
 *   var r = await src.fetch('chr1', 1000000, 1003000)     // 1-based, inclusive
 *   r.reads: [{name, start, end, strand, mapq, flag, cigar: [[op, len]...], seq}]
 * Secondary, QC-failed, duplicate and unmapped reads are left out. At most
 * maxReads (default 6,000) are kept, evenly downsampled, as IGV does; the
 * coverage counts in pileup() use all reads fetched before downsampling.
 * Very deep regions (chrM can carry tens of thousands of reads per base):
 * when the index points to more than maxBytes (default 12 MB), the start and
 * end of the region are found by bisecting the compressed file (probing which
 * read starts at an offset), and if that is still too much, 8 evenly spread
 * slices are read and the result is marked sampled: allele fractions hold,
 * absolute depth does not.
 *
 * pileup(r, ref) builds per-position counts (A, C, G, T, deletion), marks
 * mismatches against the reference bases, and packs reads into rows.
 */
(function (G) {
  var CIGAR_OPS = 'MIDNSHP=X', SEQ_CODES = '=ACMGRSVTWYHKDBN';

  function fileSource(file) {
    return { name: file.name, size: file.size, read: function (a, b) { return file.slice(a, b).arrayBuffer().then(function (x) { return new Uint8Array(x); }); } };
  }
  function urlSource(url) {
    return {
      name: url.split('/').pop(), url: url,
      read: function (a, b) {
        return fetch(url, { headers: { Range: 'bytes=' + a + '-' + (b - 1) } }).then(function (r) {
          if (!(r.ok || r.status === 206)) throw new Error('reads: HTTP ' + r.status + ' for ' + url);
          return r.arrayBuffer();
        }).then(function (x) { return new Uint8Array(x); });
      }
    };
  }

  async function inflateRaw(bytes) {
    var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  function u32(b, p) { return (b[p] | b[p + 1] << 8 | b[p + 2] << 16) + b[p + 3] * 16777216; }
  function i32(b, p) { return b[p] | b[p + 1] << 8 | b[p + 2] << 16 | b[p + 3] << 24; }

  // Inflate whole BGZF blocks from bytes read at file offset base.
  // Returns {data, blocks: [{c: compressed offset, u: offset in data}]}.
  async function inflateBlocks(bytes, base) {
    var parts = [], blocks = [], p = 0, total = 0;
    while (p + 18 <= bytes.length && bytes[p] === 31 && bytes[p + 1] === 139) {
      var xlen = bytes[p + 10] | bytes[p + 11] << 8, q = p + 12, bsize = -1;
      while (q < p + 12 + xlen) { if (bytes[q] === 66 && bytes[q + 1] === 67) bsize = bytes[q + 4] | bytes[q + 5] << 8; q += 4 + (bytes[q + 2] | bytes[q + 3] << 8); }
      if (bsize < 0) throw new Error('not a BGZF file');
      var blen = bsize + 1;
      if (p + blen > bytes.length) break; // a partial block at the end of the read
      var isize = u32(bytes, p + blen - 4);
      blocks.push({ c: base + p, u: total });
      parts.push(isize ? inflateRaw(bytes.subarray(p + 12 + xlen, p + blen - 8)) : Promise.resolve(new Uint8Array(0)));
      total += isize; p += blen;
    }
    var arrs = await Promise.all(parts), data = new Uint8Array(total), o = 0;
    arrs.forEach(function (a) { data.set(a, o); o += a.length; });
    return { data: data, blocks: blocks };
  }

  // A virtual offset from its two 32-bit halves: compressed block offset and offset inside it.
  function voff(b, p) { var lo = u32(b, p), hi = u32(b, p + 4); return { c: hi * 65536 + Math.floor(lo / 65536), u: lo % 65536 }; }

  function parseBai(b) {
    if (b[0] !== 66 || b[1] !== 65 || b[2] !== 73 || b[3] !== 1) throw new Error('not a .bai index');
    var nRef = i32(b, 4), p = 8, refs = [];
    for (var r = 0; r < nRef; r++) {
      var nBin = i32(b, p); p += 4;
      var bins = {};
      for (var k = 0; k < nBin; k++) {
        var bin = u32(b, p), nChunk = i32(b, p + 4); p += 8;
        var chunks = [];
        for (var j = 0; j < nChunk; j++) { chunks.push([voff(b, p), voff(b, p + 8)]); p += 16; }
        bins[bin] = chunks;
      }
      var nIntv = i32(b, p); p += 4;
      var lin = [];
      for (k = 0; k < nIntv; k++) { lin.push(voff(b, p)); p += 8; }
      refs.push({ bins: bins, linear: lin });
    }
    return refs;
  }

  // Bins that may hold reads overlapping 0-based [beg, end) (SAM spec reg2bins).
  function reg2bins(beg, end) {
    var list = [0], k; end--;
    for (k = 1 + (beg >> 26); k <= 1 + (end >> 26); k++) list.push(k);
    for (k = 9 + (beg >> 23); k <= 9 + (end >> 23); k++) list.push(k);
    for (k = 73 + (beg >> 20); k <= 73 + (end >> 20); k++) list.push(k);
    for (k = 585 + (beg >> 17); k <= 585 + (end >> 17); k++) list.push(k);
    for (k = 4681 + (beg >> 14); k <= 4681 + (end >> 14); k++) list.push(k);
    return list;
  }
  function vcmp(a, b) { return a.c - b.c || a.u - b.u; }

  // The first plausible record in inflated data for reference rid: a block size, ref id,
  // position and NUL-terminated name that make sense, and a next record that does too.
  function firstRecord(d, rid, refLen) {
    var ok = function (p) {
      if (p + 36 > d.length) return false;
      var size = i32(d, p), ln = d[p + 12];
      if (size < 32 || size > 1000000 || i32(d, p + 4) !== rid) return false;
      var pos = i32(d, p + 8);
      if (pos < 0 || pos > refLen || ln < 2 || p + 36 + ln > d.length || d[p + 36 + ln - 1] !== 0) return false;
      return true;
    };
    for (var p = 0; p + 36 < d.length; p++) if (ok(p) && (p + 4 + i32(d, p) >= d.length || ok(p + 4 + i32(d, p)))) return { u: p, pos: i32(d, p + 8) };
    return null;
  }

  // Where to start reading in [a, b] (virtual offsets) to reach position beg (0-based) without
  // reading everything before it: bisect on compressed offsets, probing what read starts there.
  ReadsSource.prototype.seek = async function (a, b, rid, beg) {
    var self = this, refLen = this.refs[rid].length, lo = a.c, hi = b.c, best = null;
    var probe = async function (off) {
      var bytes = await self.bam.read(off, off + 160 * 1024);
      for (var q = 0; q + 18 < bytes.length; q++) { // the next BGZF block header
        if (bytes[q] === 31 && bytes[q + 1] === 139 && bytes[q + 2] === 8 && bytes[q + 3] === 4 && bytes[q + 12] === 66 && bytes[q + 13] === 67) {
          var inf = await inflateBlocks(bytes.subarray(q), off + q), rec = firstRecord(inf.data, rid, refLen);
          return rec ? { c: off + q, u: rec.u, pos: rec.pos } : null;
        }
      }
      return null;
    };
    var margin = this.opts.margin === undefined ? 300 : this.opts.margin; // reads overlapping beg start at most this far before it (short reads)
    for (var step = 0; step < 24 && hi - lo > 64 * 1024; step++) {
      var mid = Math.floor((lo + hi) / 2), pr = await probe(mid);
      if (!pr) { hi = mid; continue; }
      if (pr.pos < beg - margin) { best = pr; lo = mid; } else hi = mid;
    }
    return best ? { c: best.c, u: best.u } : a;
  };

  function ReadsSource(bam, refs, header, index, opts) {
    this.bam = bam; this.refs = refs; this.header = header; this.index = index; this.opts = opts || {};
    this.byName = {}; var self = this;
    refs.forEach(function (r, i) { self.byName[G.genome.normName(r.name)] = i; });
    this.build = G.genome.detectBuild(refs);
    this.name = bam.name;
  }

  async function open(bam, bai, opts) {
    var index = parseBai(await bai.read(0, bai.size || 64 * 1024 * 1024));
    // the header: inflate from the start until it is complete
    var want = 256 * 1024, head;
    for (;;) {
      var inf = await inflateBlocks(await bam.read(0, want), 0), d = inf.data;
      if (d.length >= 8 && d[0] === 66 && d[1] === 65 && d[2] === 77 && d[3] === 1) {
        var lText = i32(d, 4), p = 8 + lText;
        if (d.length >= p + 4) {
          var nRef = i32(d, p), refs = [], ok = true; p += 4;
          for (var r = 0; r < nRef; r++) {
            if (p + 4 > d.length) { ok = false; break; }
            var ln = i32(d, p); if (p + 8 + ln > d.length) { ok = false; break; }
            refs.push({ name: new TextDecoder().decode(d.subarray(p + 4, p + 4 + ln - 1)), length: i32(d, p + 4 + ln) }); p += 8 + ln;
          }
          if (ok) { head = { text: new TextDecoder().decode(d.subarray(8, 8 + lText)), refs: refs }; break; }
        }
      } else throw new Error((bam.name || 'file') + ' is not a BAM file');
      if (want > 64 * 1024 * 1024) throw new Error('BAM header too large');
      want *= 4;
    }
    return new ReadsSource(bam, head.refs, head.text, index, opts);
  }

  ReadsSource.prototype.fetch = async function (chrom, start, end) {
    var rid = this.byName[G.genome.normName(chrom)], max = this.opts.maxReads || 6000;
    if (rid === undefined) return { reads: [], total: 0, chrom: chrom, start: start, end: end };
    var ix = this.index[rid], beg = start - 1, chunks = [];
    var minOff = ix.linear[beg >> 14] || { c: 0, u: 0 };
    reg2bins(beg, end).forEach(function (bin) { (ix.bins[bin] || []).forEach(function (ch) { if (vcmp(ch[1], minOff) > 0) chunks.push(ch); }); });
    chunks.sort(function (a, b) { return vcmp(a[0], b[0]); });
    var merged = []; // join chunks that touch
    chunks.forEach(function (ch) { var last = merged[merged.length - 1]; if (last && ch[0].c <= last[1].c + 65536) { if (vcmp(ch[1], last[1]) > 0) last[1] = ch[1]; } else merged.push([ch[0], ch[1]]); });
    var reads = [], seen = new Set(), total = 0, budget = this.opts.maxBytes || 12 * 1024 * 1024, capped = false, sampled = false, self = this;
    // decode the records in d from offset p (stopping at stopU or the region's end)
    var parse = function (d, p, stopU) {
      while (p + 4 <= d.length && p < stopU) {
        var size = i32(d, p);
        if (p + 4 + size > d.length) break;
        var rec = p + 4; p += 4 + size;
        if (i32(d, rec) !== rid) continue;
        var pos = i32(d, rec + 4), lName = d[rec + 8], mapq = d[rec + 9], nCig = d[rec + 12] | d[rec + 13] << 8, flag = d[rec + 14] | d[rec + 15] << 8, lSeq = i32(d, rec + 16);
        if (pos >= end) break;
        if (flag & (4 | 256 | 512 | 1024)) continue;
        var cp = rec + 32 + lName, cigar = [], refLen = 0;
        for (var k = 0; k < nCig; k++) { var v = u32(d, cp + 4 * k), op = CIGAR_OPS[v & 15], len = v >>> 4; cigar.push([op, len]); if ('MDN=X'.indexOf(op) >= 0) refLen += len; }
        if (pos + refLen <= beg) continue; // ends before the region
        var name = new TextDecoder().decode(d.subarray(rec + 32, rec + 32 + lName - 1));
        var key = name + ':' + pos + ':' + (flag & 2048);
        if (seen.has(key)) continue;
        seen.add(key); total++;
        var sp = cp + 4 * nCig, seq = new Array(lSeq);
        for (k = 0; k < lSeq; k++) seq[k] = SEQ_CODES[(d[sp + (k >> 1)] >> ((1 - (k & 1)) * 4)) & 15];
        reads.push({ name: name, start: pos + 1, end: pos + refLen, strand: flag & 16 ? -1 : 1, mapq: mapq, flag: flag, cigar: cigar, seq: seq.join('') });
      }
    };
    for (var m = 0; m < merged.length; m++) {
      var a = merged[m][0], b = merged[m][1];
      if (b.c - a.c > budget) { // deep coverage: find the region's start and end inside the range
        var s0 = await this.seek(a, b, rid, beg), s1 = await this.seek(s0, b, rid, end + (this.opts.margin === undefined ? 300 : this.opts.margin));
        if (s1.c - s0.c > budget) { // still too much to read whole: an even sample of 8 slices
          sampled = true;
          var K = 8, slice = Math.floor(budget / K);
          for (var k = 0; k < K; k++) {
            var off = Math.floor(s0.c + (s1.c - s0.c) * k / K), bytes = await this.bam.read(off, off + slice), q = 0;
            while (q + 18 < bytes.length && !(bytes[q] === 31 && bytes[q + 1] === 139 && bytes[q + 2] === 8 && bytes[q + 3] === 4 && bytes[q + 12] === 66 && bytes[q + 13] === 67)) q++;
            var infK = await inflateBlocks(bytes.subarray(q), off + q), first = firstRecord(infK.data, rid, this.refs[rid].length);
            if (first) parse(infK.data, first.u, infK.data.length);
          }
          continue;
        }
        a = s0; b = s1;
      }
      var inf = await inflateBlocks(await this.bam.read(a.c, b.c + 65536 + 28), a.c);
      var blockAt = {}; inf.blocks.forEach(function (bl) { blockAt[bl.c] = bl.u; });
      parse(inf.data, (blockAt[a.c] || 0) + a.u, blockAt[b.c] !== undefined ? blockAt[b.c] + b.u : inf.data.length);
    }
    reads.sort(function (x, y) { return x.start - y.start; });
    var kept = reads;
    if (reads.length > max) { var step = reads.length / max; kept = []; for (var t = 0; t < max; t++) kept.push(reads[Math.floor(t * step)]); }
    return { reads: kept, all: reads, total: total, downsampled: reads.length > max, sampled: sampled, chrom: chrom, start: start, end: end };
  };

  // Per-position counts over [r.start, r.end] from all fetched reads, mismatches against
  // ref (bases of [r.start, r.end]), and kept reads packed into rows.
  function pileup(r, ref) {
    var n = r.end - r.start + 1, cov = [], i;
    for (i = 0; i < n; i++) cov.push({ A: 0, C: 0, G: 0, T: 0, N: 0, del: 0, ins: 0, depth: 0 });
    var walk = function (rd, cb) {
      var rp = rd.start, qp = 0;
      rd.cigar.forEach(function (c) {
        var op = c[0], len = c[1];
        if (op === 'M' || op === '=' || op === 'X') { for (var k = 0; k < len; k++) cb('base', rp + k, rd.seq[qp + k]); rp += len; qp += len; }
        else if (op === 'D' || op === 'N') { if (op === 'D') for (var k2 = 0; k2 < len; k2++) cb('del', rp + k2); rp += len; }
        else if (op === 'I') { cb('ins', rp, len); qp += len; }
        else if (op === 'S') qp += len;
      });
    };
    (r.all || r.reads).forEach(function (rd) {
      walk(rd, function (kind, p, b) {
        var j = p - r.start; if (j < 0 || j >= n) return;
        if (kind === 'base') { var key = 'ACGT'.indexOf(b) >= 0 ? b : 'N'; cov[j][key]++; cov[j].depth++; }
        else if (kind === 'del') { cov[j].del++; cov[j].depth++; }
        else cov[j].ins++;
      });
    });
    // kept reads: mismatches, deletions and insertions to draw; rows packed greedily
    var rowsEnd = [], rows = [];
    r.reads.forEach(function (rd) {
      var marks = [];
      walk(rd, function (kind, p, b) {
        if (kind === 'base') { var rb = ref ? ref[p - r.start] : null; if (rb && b !== rb && b !== 'N' && rb !== 'N') marks.push(['x', p, b]); }
        else if (kind === 'del') marks.push(['d', p]);
        else marks.push(['i', p, b]);
      });
      rd.marks = marks;
      var row = rowsEnd.findIndex(function (e) { return e < rd.start - 2; });
      if (row < 0) { row = rowsEnd.length; rowsEnd.push(0); rows.push([]); }
      rowsEnd[row] = rd.end; rows[row].push(rd);
    });
    return { cov: cov, rows: rows, start: r.start, end: r.end, ref: ref };
  }

  G.reads = { fileSource: fileSource, urlSource: urlSource, open: open, pileup: pileup, parseBai: parseBai, reg2bins: reg2bins, inflateBlocks: inflateBlocks };
})(globalThis.G = globalThis.G || {});
