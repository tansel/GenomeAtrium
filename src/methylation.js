/*
 * Methylation overlay, shown only when a methylation file is opened on top
 * of a loaded genome. Formats:
 *  - bedMethyl from modkit (nanopore or PacBio; 18 columns): 0-based start,
 *    mod code in column 4 (m 5mC, h 5hmC, a 6mA), Nvalid_cov in column 10,
 *    Nmod in column 12. Columns 10 and 11 may be space-separated (older modkit).
 *  - ENCODE-style bedMethyl (bisulfite; 11 columns): coverage in column 10,
 *    percent in column 11; counted as 5mC (bisulfite reads 5mC and 5hmC together).
 *  - modbam2bed (older ONT pipelines, e.g. wf-human-variation; 14 columns):
 *    coverage in column 10 counts filtered calls too, so the valid calls are
 *    Ncanonical (column 12) plus Nmod (column 13).
 *  - Bismark coverage (.cov): 1-based start, percent, methylated, unmethylated.
 * A file named with hap1 or hap2 is one haplotype; the two load as a pair, and
 * the checks then compare the copies (imprinting: one methylated, one not).
 * Sites are summed into the genome's bins: methylated calls and valid calls
 * per bin, per mod code, so a region's level is a coverage-weighted mean.
 * Single sites are not kept (a whole-genome file has about 28 million CpGs).
 *
 * The file carries no genome build. It is accepted only when its positions fit
 * the loaded genome's contigs and its name does not name another build.
 */
(function (G) {
  var MIN_COV = 5; // calls per region before a level is shown

  async function parse(blob, data, opts) {
    opts = opts || {};
    var gm = G.genome, genome = data.genome, bs = data.binSize, name = (opts.name || blob.name || '').toLowerCase();
    var hint = /hg19|grch37|b37/.test(name) ? 'GRCh37' : /hg38|grch38/.test(name) ? 'GRCh38' : null;
    if (hint && data.build && hint !== data.build && !(hint === 'GRCh37' && /hg19|GRCh37/.test(data.build)))
      throw new Error('The file name says ' + hint + ' but the loaded genome is ' + data.build + '. Coordinates are never compared across builds.');
    var mods = {}, stats = { lines: 0, sites: 0, skipped: 0, outside: 0, unknownContig: 0, format: null };
    var tracksFor = function (code, c) {
      var m = mods[code] || (mods[code] = {});
      return m[c.key] || (m[c.key] = { meth: new gm.Track(c.length, bs), cov: new gm.Track(c.length, bs), sites: new gm.Track(c.length, bs) });
    };
    for await (var line of G.bgzf.lines(blob, { signal: opts.signal })) {
      if (!line || line[0] === '#' || line.startsWith('track') || line.startsWith('browser')) continue;
      stats.lines++;
      var f = line.split(/\s+/), chrom = f[0], pos, code = 'm', cov, nmod;
      if (f.length === 14 && /^\d+$/.test(f[11]) && /^\d+$/.test(f[12])) { // modbam2bed
        pos = +f[1] + 1; code = /^[a-zA-Z]$/.test(f[3]) ? f[3] : /5hmC/i.test(f[3]) ? 'h' : 'm';
        nmod = +f[12]; cov = +f[11] + nmod;
        stats.format = stats.format || 'modbam2bed';
      } else if (f.length >= 11) {
        // bedMethyl: modkit (18 columns) or ENCODE (11)
        pos = +f[1] + 1; code = /^[a-zA-Z]$/.test(f[3]) ? f[3] : 'm';
        cov = +f[9];
        nmod = f.length >= 18 ? +f[11] : Math.round(cov * (+f[10]) / 100);
        stats.format = stats.format || (f.length >= 18 ? 'modkit bedMethyl' : 'bedMethyl (bisulfite)');
      } else if (f.length === 6) { // Bismark coverage
        pos = +f[1]; nmod = +f[4]; cov = nmod + (+f[5]);
        stats.format = stats.format || 'Bismark coverage';
      } else { stats.skipped++; continue; }
      if (!(cov > 0) || !isFinite(pos)) { stats.skipped++; continue; }
      var c = genome.get(chrom);
      if (!c) { stats.unknownContig++; continue; }
      if (pos > c.length) { stats.outside++; continue; }
      var t = tracksFor(code, c);
      t.meth.add(pos, nmod); t.cov.add(pos, cov); t.sites.add(pos, 1);
      stats.sites++;
      if (opts.onProgress && stats.lines % 200000 === 0) opts.onProgress(stats);
    }
    if (stats.outside > 0.01 * (stats.sites + stats.outside)) throw new Error(stats.outside + ' sites lie past the ends of this genome\'s contigs: the file is probably on another build.');
    if (!stats.sites) throw new Error('No methylation sites matched the loaded genome' + (stats.unknownContig ? ' (contig names did not match)' : '') + '.');
    Object.keys(mods).forEach(function (code) { Object.keys(mods[code]).forEach(function (k) { var t = mods[code][k]; t.meth.buildPyramid(); t.cov.buildPyramid(); t.sites.buildPyramid(); }); });
    var codes = Object.keys(mods).sort(function (a, b) { return (b === 'm') - (a === 'm'); });
    var m = new Methylation(mods, codes, stats, data, opts.name || blob.name);
    var hp = /hap(?:lotype)?[._-]?([12])\b/i.exec(opts.name || blob.name || '');
    m.hap = hp ? +hp[1] : null;
    return m;
  }

  // Two haplotypes as one overlay: tracks summed for the band and genome-wide numbers,
  // the parts kept for the per-copy checks.
  function combine(a, b) {
    var gm = G.genome, mods = {}, data = a.data;
    [a, b].forEach(function (x) {
      Object.keys(x.mods).forEach(function (code) {
        Object.keys(x.mods[code]).forEach(function (k) {
          var src = x.mods[code][k], c = data.genome.contigs.find(function (cc) { return cc.key === k; });
          var dst = (mods[code] = mods[code] || {})[k] || (mods[code][k] = { meth: new gm.Track(c.length, src.cov.binSize), cov: new gm.Track(c.length, src.cov.binSize), sites: new gm.Track(c.length, src.cov.binSize) });
          ['meth', 'cov', 'sites'].forEach(function (t) {
            var from = src[t].levels[0], to = dst[t]._fit(from.length - 1);
            for (var i = 0; i < from.length; i++) to[i] += from[i];
          });
        });
      });
    });
    Object.keys(mods).forEach(function (code) { Object.keys(mods[code]).forEach(function (k) { ['meth', 'cov', 'sites'].forEach(function (t) { mods[code][k][t].buildPyramid(); }); }); });
    var stats = { sites: a.stats.sites + b.stats.sites, format: a.stats.format + ', two haplotypes', lines: a.stats.lines + b.stats.lines };
    var m = new Methylation(mods, Object.keys(mods), stats, data, a.fileName + ' + ' + b.fileName);
    m.parts = a.hap === 2 ? { 1: b, 2: a } : { 1: a, 2: b };
    return m;
  }

  function Methylation(mods, codes, stats, data, fileName) {
    this.mods = mods; this.codes = codes; this.code = codes[0]; this.stats = stats; this.data = data; this.fileName = fileName;
    this.summary = this.globalSummary();
  }

  Methylation.CODE_NAMES = { m: '5mC', h: '5hmC', a: '6mA', c: '4mC', C: 'any C mod', A: 'any A mod' };

  Methylation.prototype.tracks = function (key, code) { var m = this.mods[code || this.code]; return m ? m[key] : null; };

  // Level over [start, end] on contig key: { frac, cov, sites } (frac null below MIN_COV).
  Methylation.prototype.level = function (key, start, end, code) {
    var t = this.tracks(key, code);
    if (!t) return { frac: null, cov: 0, sites: 0 };
    var bs = t.cov.binSize, b0 = Math.max(0, Math.floor((start - 1) / bs)), b1 = Math.floor((end - 1) / bs), m = 0, c = 0, s = 0;
    var M = t.meth.levels[0], C = t.cov.levels[0], S = t.sites.levels[0];
    for (var b = b0; b <= b1 && b < C.length; b++) { m += M[b]; c += C[b]; s += S[b]; }
    return { frac: c >= MIN_COV ? m / c : null, cov: c, sites: s };
  };

  // Genome-wide: share of well-covered bins that are low (<20%), middle, high (>70%).
  Methylation.prototype.globalSummary = function () {
    var self = this, lo = 0, mid = 0, hi = 0, mSum = 0, cSum = 0;
    this.data.genome.contigs.forEach(function (c) {
      var t = self.tracks(c.key);
      if (!t) return;
      var M = t.meth.levels[0], C = t.cov.levels[0];
      for (var b = 0; b < C.length; b++) {
        mSum += M[b]; cSum += C[b];
        if (C[b] < 20) continue;
        var f = M[b] / C[b];
        if (f < 0.2) lo++; else if (f > 0.7) hi++; else mid++;
      }
    });
    var n = lo + mid + hi || 1;
    return { mean: cSum ? mSum / cSum : null, low: lo / n, mid: mid / n, high: hi / n, bins: lo + mid + hi };
  };

  // Places a geneticist checks: imprinted genes, whose control regions sit at or near
  // these promoters and are methylated on one parental copy (about 50% expected), and
  // FMR1, whose promoter is unmethylated on an active X (a full Fragile X expansion
  // methylates it). The window is the promoter, TSS +-1 kb (the H19 control region lies
  // 2 to 4 kb upstream, so H19 uses that window instead). A stand-in for exact DMR
  // coordinates, read from the genome's bins (about 1 kb).
  var CHECKS = [
    ['SNRPN', 'imprinted (15q11-q13: Prader-Willi, Angelman)', 'half'], ['KCNQ1OT1', 'imprinted (11p15 KvDMR: Beckwith-Wiedemann)', 'half'],
    ['H19', 'imprinted (11p15 H19/IGF2 ICR: Beckwith-Wiedemann, Silver-Russell)', 'half', [-4000, -2000]],
    ['MEG3', 'imprinted (14q32: Temple, Kagami-Ogata)', 'half'], ['PEG3', 'imprinted (19q13)', 'half'], ['PLAGL1', 'imprinted (6q24: transient neonatal diabetes)', 'half'],
    ['MEST', 'imprinted (7q32: Silver-Russell)', 'half'], ['GNAS', 'imprinted (20q13: pseudohypoparathyroidism)', 'half'],
    ['FMR1', 'Fragile X promoter: unmethylated on an active X', 'low']
  ];
  Methylation.prototype.checks = function (genes) {
    var self = this, d = this.data;
    if (!genes) return [];
    return CHECKS.map(function (ck) {
      var g = genes.get(ck[0]);
      if (!g) return null;
      var c = d.genome.get(g.chrom);
      if (!c) return null;
      var plus = g.strand !== '-' && g.strand !== -1, tss = plus ? g.start : g.end, w = ck[3] || [-1000, 1000]; // strand is 1/-1 in the genes table
      var a = plus ? tss + w[0] : tss - w[1], b = plus ? tss + w[1] : tss - w[0];
      var lv = self.level(c.key, a, b);
      var row = { gene: ck[0], role: ck[1], chrom: c.name, start: a, end: b, frac: lv.frac, cov: lv.cov, sites: lv.sites };
      if (self.parts && ck[2] === 'half') { // per copy: imprinting shows as one copy methylated, the other not
        var h1 = self.parts[1].level(c.key, a, b).frac, h2 = self.parts[2].level(c.key, a, b).frac;
        row.hap = [h1, h2];
        // The control region need not sit at the gene's outermost start (PLAGL1, MEST and GNAS have
        // several promoters), so also scan the locus for its strongest stretch where the copies differ.
        var scan = self.asmScan(c.key, Math.min(g.start, g.end) - 10000, Math.max(g.start, g.end) + 10000);
        row.scan = scan;
        if (scan && (h1 === null || h2 === null || scan.diff > Math.abs(h1 - h2))) {
          row.hap = scan.hap; row.start = scan.start; row.end = scan.end; row.frac = (scan.hap[0] + scan.hap[1]) / 2;
          h1 = scan.hap[0]; h2 = scan.hap[1]; row.where = 'strongest stretch in the locus, ' + Math.round((scan.end - scan.start + 1) / 1000) + ' kb';
        } else row.where = 'promoter window';
        if (h1 !== null && h2 !== null) {
          row.asm = Math.abs(h1 - h2);
          row.flag = row.asm >= 0.5 ? 'one copy methylated, the other not: the imprinting pattern' : row.asm >= 0.25 ? 'copies differ partly' : 'both copies alike: no imprinting pattern here';
          row.unusual = row.asm < 0.25;
          return row;
        }
      }
      var flag = lv.frac === null ? 'too few calls' : ck[2] === 'half' ? (lv.frac < 0.3 ? 'low for an imprinted region' : lv.frac > 0.7 ? 'high for an imprinted region' : 'about half, as expected')
        : (lv.frac > 0.5 ? 'methylated: check for a full expansion (or an inactive X in a female)' : 'unmethylated or partly, as on an active X');
      row.flag = flag; row.unusual = lv.frac !== null && (ck[2] === 'half' ? (lv.frac < 0.3 || lv.frac > 0.7) : lv.frac > 0.5);
      return row;
    }).filter(Boolean);
  };

  // Strongest allele-specific stretch in [start, end]: consecutive bins where both copies have at
  // least 10 calls and differ by 0.4 or more; returns the run with the largest summed difference.
  Methylation.prototype.asmScan = function (key, start, end) {
    if (!this.parts) return null;
    var p1 = this.parts[1].tracks(key), p2 = this.parts[2].tracks(key);
    if (!p1 || !p2) return null;
    var bs = p1.cov.binSize, b0 = Math.max(0, Math.floor((start - 1) / bs)), b1 = Math.floor((end - 1) / bs);
    var M1 = p1.meth.levels[0], C1 = p1.cov.levels[0], M2 = p2.meth.levels[0], C2 = p2.cov.levels[0];
    var best = null, run = null;
    var close = function () { if (run && (!best || run.score > best.score)) best = run; run = null; };
    for (var b = b0; b <= b1; b++) {
      var ok = C1[b] >= 10 && C2[b] >= 10, f1 = ok ? M1[b] / C1[b] : 0, f2 = ok ? M2[b] / C2[b] : 0;
      if (ok && Math.abs(f1 - f2) >= 0.4 && (!run || Math.sign(f1 - f2) === run.sign)) {
        if (!run) run = { b0: b, sign: Math.sign(f1 - f2), score: 0, m1: 0, c1: 0, m2: 0, c2: 0 };
        run.b1 = b; run.score += Math.abs(f1 - f2); run.m1 += M1[b]; run.c1 += C1[b]; run.m2 += M2[b]; run.c2 += C2[b];
      } else close();
    }
    close();
    if (!best) return null;
    var h = [best.m1 / best.c1, best.m2 / best.c2];
    return { start: best.b0 * bs + 1, end: (best.b1 + 1) * bs, hap: h, diff: Math.abs(h[0] - h[1]), bins: best.b1 - best.b0 + 1 };
  };

  // Blue (unmethylated) to red (methylated).
  Methylation.color = function (f, alpha) {
    var r = Math.round(60 + 195 * f), gr = Math.round(120 - 60 * Math.abs(f - 0.5) * 2 + 40 * (1 - f)), b = Math.round(255 - 205 * f);
    return 'rgba(' + r + ',' + gr + ',' + b + ',' + (alpha === undefined ? 1 : alpha) + ')';
  };

  G.methylation = { parse: parse, combine: combine, Methylation: Methylation, MIN_COV: MIN_COV };
})(globalThis.G = globalThis.G || {});
