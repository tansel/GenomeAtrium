/*
 * Streaming VCF / gVCF parser (plain or bgzipped).
 *
 * Rules carried over from Asclepius (design/decisions D2 to D4):
 *  - A variant is CHROM+POS+REF+ALT, never position alone.
 *  - A missing call (./.) is not a confident 0/0; they get separate states.
 *  - Contig names are compared after normalisation (chr1 == 1, chrM == MT).
 *  - gVCF reference blocks (ALT "." or <NON_REF>) are not variants. They are
 *    kept as a "callable" track so gaps in coverage stay visible.
 */
(function (G) {
  var T = {        // variant type codes
    SNV: 0, MNV: 1, INS: 2, DEL: 3, COMPLEX: 4,
    SV_DEL: 5, SV_DUP: 6, SV_INV: 7, SV_INS: 8, SV_CNV: 9, BND: 10
  };
  var TYPE_NAMES = ['SNV', 'MNV', 'insertion', 'deletion', 'complex',
    'SV deletion', 'SV duplication', 'SV inversion', 'SV insertion', 'CNV', 'breakend'];
  var Z = { UNKNOWN: 0, HET: 1, HOM: 2, REF: 3, MISSING: 4 };
  var ZYG_NAMES = ['no genotype', 'heterozygous', 'homozygous alt', 'reference (0/0)', 'missing call'];

  var LABEL_LIMIT = 1500000; // keep REF>ALT strings for hover below this count
  var MAX_SAMPLES = 16;      // genotypes kept for the first 16 samples
  var BIG_INDEL = 50;        // sequence-resolved indels this long also become arcs

  // Growable typed columns for one contig's variants.
  function Columns() {
    this.n = 0; this.cap = 1024;
    this.pos = new Int32Array(this.cap);
    this.type = new Uint8Array(this.cap);
    this.zyg = new Uint8Array(this.cap);
    this.pass = new Uint8Array(this.cap);
    this.snv = new Uint8Array(this.cap); // 16 | ref << 2 | alt for a plain SNV, else 0
    this.baf = new Uint8Array(this.cap); // alt allele fraction from FORMAT AD, as 1 + round(250 * f); 0 = unknown
    this.dp = new Uint16Array(this.cap); // read depth from FORMAT DP (or the AD sum), capped at 65535
    this.zygS = [];                      // per extra sample (2nd, 3rd, ...): Uint8Array of Z codes
    this.alleleText = new Map();         // index -> "REF>ALT[,ALT]" for everything else
    this.labels = [];
  }
  var BASES = 'ACGT', BASE_CODE = { A: 0, C: 1, G: 2, T: 3 };
  // REF and ALT alleles of variant i, for every variant (labels stop at 1.5M).
  Columns.prototype.alleles = function (i) {
    var c = this.snv[i];
    if (c) return { ref: BASES[(c >> 2) & 3], alts: [BASES[c & 3]] };
    var t = this.alleleText.get(i);
    if (!t) return null;
    var k = t.indexOf('>');
    return { ref: t.slice(0, k), alts: t.slice(k + 1).split(',') };
  };
  // Alt allele fraction of variant i (0..1), or null when the file has no AD.
  Columns.prototype.bafAt = function (i) { var b = this.baf[i]; return b ? (b - 1) / 250 : null; };

  // Zygosity of variant i in sample s (0 = the first sample).
  Columns.prototype.zygOf = function (i, s) { return s ? (this.zygS[s - 1] ? this.zygS[s - 1][i] : Z.UNKNOWN) : this.zyg[i]; };

  Columns.prototype.push = function (pos, type, zyg, pass, label, ref, altField, altFrac, depth, extraZyg) {
    if (this.n === this.cap) {
      this.cap *= 2;
      var grow = function (a, C) { var b = new C(this.cap); b.set(a); return b; }.bind(this);
      this.pos = grow(this.pos, Int32Array);
      this.type = grow(this.type, Uint8Array);
      this.zyg = grow(this.zyg, Uint8Array);
      this.pass = grow(this.pass, Uint8Array);
      this.snv = grow(this.snv, Uint8Array);
      this.baf = grow(this.baf, Uint8Array);
      this.dp = grow(this.dp, Uint16Array);
      for (var zs = 0; zs < this.zygS.length; zs++) this.zygS[zs] = grow(this.zygS[zs], Uint8Array);
    }
    if (ref !== undefined) {
      var rc = BASE_CODE[ref], ac = altField.length === 1 ? BASE_CODE[altField] : undefined;
      if (rc !== undefined && ac !== undefined) this.snv[this.n] = 16 | (rc << 2) | ac;
      else this.alleleText.set(this.n, G.bgzf.own(ref + '>' + altField));
    }
    if (extraZyg) for (var es = 0; es < extraZyg.length; es++) {
      if (!this.zygS[es]) this.zygS[es] = new Uint8Array(this.cap);
      this.zygS[es][this.n] = extraZyg[es];
    }
    if (altFrac != null) this.baf[this.n] = 1 + Math.round(250 * altFrac);
    if (depth) this.dp[this.n] = Math.min(65535, depth);
    this.pos[this.n] = pos; this.type[this.n] = type;
    this.zyg[this.n] = zyg; this.pass[this.n] = pass;
    this.labels.push(G.bgzf.own(label));
    this.n++;
  };

  function infoValue(info, key) {
    var s = ';' + info + ';';
    var i = s.indexOf(';' + key + '=');
    if (i < 0) return null;
    var start = i + key.length + 2;
    return s.slice(start, s.indexOf(';', start));
  }

  function zygosity(gt) {
    if (gt == null || gt === '') return Z.UNKNOWN;
    var al = gt.split(/[\/|]/);
    var missing = 0, zero = 0, first = null, same = true;
    for (var i = 0; i < al.length; i++) {
      if (al[i] === '.') { missing++; continue; }
      if (al[i] === '0') zero++;
      if (first === null) first = al[i]; else if (al[i] !== first) same = false;
    }
    if (missing === al.length) return Z.MISSING;
    if (missing) return Z.MISSING; // half calls (./1) are not trusted either way
    if (zero === al.length) return Z.REF;
    return same ? Z.HOM : Z.HET;
  }

  function classifySeq(ref, alt) {
    if (ref.length === 1 && alt.length === 1) return T.SNV;
    if (ref.length === alt.length) return T.MNV;
    if (ref.length === 1 && alt[0] === ref[0]) return T.INS;
    if (alt.length === 1 && ref[0] === alt[0]) return T.DEL;
    return alt.length > ref.length ? T.INS : ref.length > alt.length ? T.DEL : T.COMPLEX;
  }

  var SYMBOLIC = { DEL: T.SV_DEL, DUP: T.SV_DUP, INV: T.SV_INV, INS: T.SV_INS, CNV: T.SV_CNV };
  var BND_RE = /[\[\]]([^:\[\]]+):(\d+)[\[\]]/;

  async function parse(blob, opts) {
    opts = opts || {};
    var gm = G.genome;
    var genome = new gm.Genome([]);
    var headerContigs = [], samples = [], meta = {};
    var binSize = 0;
    var perContig = {};  // key -> { cols, tracks }
    var arcs = [], arcKeys = {};
    var stats = { lines: 0, variants: 0, refBlocks: 0, filtered: 0, labelsDropped: false,
                  byType: new Array(TYPE_NAMES.length).fill(0), byZyg: new Array(ZYG_NAMES.length).fill(0) };
    var total = 0;
    // ClinVar matching only when the file is on the same build (Asclepius D3).
    var clinvar = opts.clinvar || null, clinvarStatus = clinvar ? null : 'no ClinVar table loaded';
    var clinvarHits = [];

    function contigData(name) {
      var key = gm.normName(name);
      var d = perContig[key];
      if (d) return d;
      if (!binSize) binSize = gm.chooseBinSize(headerContigs.reduce(function (s, c) { return s + c.length; }, 0) || 4e9);
      var c = genome.get(name) || genome.addContig(G.bgzf.own(name), 0);
      var len = Math.max(c.length, binSize);
      d = perContig[key] = {
        contig: c,
        cols: new Columns(),
        snv: new gm.Track(len, binSize), indel: new gm.Track(len, binSize), sv: new gm.Track(len, binSize),
        het: new gm.Track(len, binSize), hom: new gm.Track(len, binSize),
        callable: new gm.Track(len, binSize, 'mean'), lowdp: new gm.Track(len, binSize, 'mean'),
        maxEnd: 0
      };
      return d;
    }

    function addArc(c0, p0, c1, p1, type, label, pass) {
      var a = [gm.normName(c0), p0], b = [gm.normName(c1), p1];
      var key = a[0] < b[0] || (a[0] === b[0] && a[1] <= b[1]) ? a.concat(b).join(':') : b.concat(a).join(':');
      if (arcKeys[key]) return;
      arcKeys[key] = 1;
      var own = G.bgzf.own;
      arcs.push({ c0: own(c0), p0: p0, c1: own(c1), p1: p1, type: type, label: own(label), pass: pass, support: 1 });
    }

    for await (var line of G.bgzf.lines(blob, opts)) {
      if (line.charCodeAt(0) === 35) { // '#'
        if (line.startsWith('##contig=')) {
          var id = /ID=([^,>]+)/.exec(line), ln = /length=(\d+)/.exec(line);
          if (id) { var cn = G.bgzf.own(id[1]); headerContigs.push({ name: cn, length: ln ? +ln[1] : 0 }); genome.addContig(cn, ln ? +ln[1] : 0); }
        } else if (line.startsWith('##reference=') || line.startsWith('##fileformat=') || line.startsWith('##source=')) {
          var eq = line.indexOf('=');
          meta[G.bgzf.own(line.slice(2, eq))] = G.bgzf.own(line.slice(eq + 1));
        } else if (line.startsWith('#CHROM')) {
          samples = line.split('\t').slice(9).map(G.bgzf.own);
        }
        continue;
      }
      if (!line) continue;
      if (clinvar && clinvarStatus === null) {
        var fileBuild = gm.detectBuild(headerContigs);
        clinvarStatus = fileBuild === clinvar.build ? 'matched' :
          fileBuild ? 'not matched: file is ' + fileBuild + ', ClinVar table is ' + clinvar.build :
          'not matched: build unknown (no ##contig lengths)';
        if (clinvarStatus !== 'matched') clinvar = null;
      }
      stats.lines++;
      var f = line.split('\t');
      if (f.length < 8) continue;
      var chrom = f[0], pos = +f[1], ref = f[3], altField = f[4];
      var pass = f[6] === 'PASS' || f[6] === '.' ? 1 : 0;
      var d = contigData(chrom);

      var alts = altField === '.' ? [] : altField.split(',').filter(function (a) {
        return a !== '<NON_REF>' && a !== '<*>' && a !== '*';
      });

      if (!alts.length) { // gVCF reference block or monomorphic site
        var endV = infoValue(f[7], 'END');
        var end = endV ? +endV : pos + ref.length - 1;
        (pass ? d.callable : d.lowdp).addSpan(pos, end, 1);
        if (end > d.maxEnd) d.maxEnd = end;
        stats.refBlocks++;
        continue;
      }

      var gt = null, altFrac = null, depth = 0;
      if (f.length > 9) {
        var fmt = f[8].split(':'), sv = f[9].split(':'), gi = fmt.indexOf('GT');
        if (gi >= 0) gt = sv[gi];
        // B-allele fraction: all non-reference depth over total depth, from AD
        var ai = fmt.indexOf('AD'), di = fmt.indexOf('DP');
        if (ai >= 0 && sv[ai] && sv[ai] !== '.') {
          var ad = sv[ai].split(','), refD = +ad[0] || 0, altD = 0;
          for (var q = 1; q < ad.length; q++) altD += +ad[q] || 0;
          if (refD + altD > 0) { altFrac = altD / (refD + altD); depth = refD + altD; }
        }
        if (di >= 0 && sv[di] && sv[di] !== '.') depth = +sv[di] || depth;
      }
      var zyg = zygosity(gt), extraZyg = null;
      if (f.length > 10) { // further samples, up to MAX_SAMPLES
        var gIdx = f[8].split(':').indexOf('GT');
        extraZyg = [];
        for (var sI = 10; sI < f.length && sI < 9 + MAX_SAMPLES; sI++) extraZyg.push(gIdx >= 0 ? zygosity(f[sI].split(':')[gIdx]) : Z.UNKNOWN);
      }
      if (clinvar) {
        var hits = clinvar.matchRecord(chrom, pos, ref, altField, gt, f[6]);
        if (hits) for (var h = 0; h < hits.length; h++) clinvarHits.push(hits[h]);
      }

      var alt = alts[0], type, svEnd = null, label;
      if (alt[0] === '<') {
        var sym = alt.slice(1, -1).split(':')[0];
        type = SYMBOLIC[sym] !== undefined ? SYMBOLIC[sym] : T.COMPLEX;
        var e = infoValue(f[7], 'END'), svlen = infoValue(f[7], 'SVLEN');
        svEnd = e ? +e : svlen ? pos + Math.abs(parseInt(svlen, 10)) : pos;
      } else if (BND_RE.test(alt)) {
        type = T.BND;
      } else if (alt.indexOf('.') >= 0) {
        type = T.COMPLEX; // single breakend
      } else {
        type = classifySeq(ref, alt);
        if ((type === T.DEL || type === T.INS) && Math.abs(ref.length - alt.length) >= BIG_INDEL) {
          svEnd = type === T.DEL ? pos + ref.length - 1 : pos + 1;
        }
      }

      stats.variants++;
      stats.byType[type]++;
      stats.byZyg[zyg]++;
      if (!pass) stats.filtered++;

      if (stats.variants <= LABEL_LIMIT) {
        var shortRef = ref.length > 12 ? ref.slice(0, 10) + '..(' + ref.length + ')' : ref;
        var shortAlt = altField.length > 24 ? altField.slice(0, 22) + '..' : altField;
        label = (f[2] !== '.' ? f[2] + ' ' : '') + shortRef + '>' + shortAlt + (pass ? '' : ' [' + f[6] + ']') +
                (gt ? ' GT=' + gt : '') + (f[5] !== '.' ? ' Q=' + f[5] : '');
      } else {
        label = null; stats.labelsDropped = true;
      }
      d.cols.push(pos, type, zyg, pass, label, ref, alts.join(','), altFrac, depth, extraZyg);

      if (type === T.SNV || type === T.MNV) d.snv.add(pos, 1);
      else if (type <= T.COMPLEX) d.indel.add(pos, 1);
      else d.sv.add(pos, 1);
      if (zyg === Z.HET) d.het.add(pos, 1); else if (zyg === Z.HOM) d.hom.add(pos, 1);
      d.callable.addSpan(pos, pos + ref.length - 1, 1);
      var last = svEnd && svEnd > pos ? svEnd : pos + ref.length - 1;
      if (last > d.maxEnd) d.maxEnd = last;

      if (type !== T.BND && !(svEnd !== null && svEnd > pos)) continue;
      var svtype = infoValue(f[7], 'SVTYPE');
      var svLabel = TYPE_NAMES[type] + (svtype ? ' (' + svtype + ')' : '') + ' ' + chrom + ':' + pos;
      if (type === T.BND) {
        var m = BND_RE.exec(alt);
        addArc(chrom, pos, m[1], +m[2], type, svLabel + ' to ' + m[1] + ':' + m[2], pass);
        contigData(m[1]); // make sure the mate's contig exists
      } else if (svEnd !== null && svEnd > pos) {
        addArc(chrom, pos, chrom, svEnd, type, svLabel + '-' + svEnd + ' (' + (svEnd - pos + 1).toLocaleString() + ' bp)', pass);
      }
    }

    // Finalise contig lengths, drop empty non-primary contigs, build pyramids.
    Object.keys(perContig).forEach(function (k) {
      var d = perContig[k];
      if (d.maxEnd > d.contig.length) d.contig.length = d.maxEnd;
    });
    genome.finish(function (c) { return !!perContig[c.key]; });
    var tracks = {}, variants = {};
    Object.keys(perContig).forEach(function (k) {
      var d = perContig[k];
      ['snv', 'indel', 'sv', 'het', 'hom', 'callable', 'lowdp'].forEach(function (t) { d[t].buildPyramid(); });
      tracks[k] = { snv: d.snv, indel: d.indel, sv: d.sv, het: d.het, hom: d.hom, callable: d.callable, lowdp: d.lowdp };
      variants[k] = d.cols;
    });

    return {
      format: 'vcf', genome: genome, build: genome.build, binSize: binSize || 1000,
      samples: samples, meta: meta, stats: stats, isGvcf: stats.refBlocks > 0,
      variants: variants, tracks: tracks, arcs: arcs,
      clinvarHits: clinvarHits, clinvarStatus: clinvarStatus || 'not matched: no variant records',
      aborted: !!(opts.signal && opts.signal.aborted)
    };
  }

  G.vcf = { parse: parse, T: T, TYPE_NAMES: TYPE_NAMES, Z: Z, ZYG_NAMES: ZYG_NAMES, zygosity: zygosity, classifySeq: classifySeq };
})(globalThis.G = globalThis.G || {});
