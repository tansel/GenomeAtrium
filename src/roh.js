/*
 * Runs of homozygosity (ROH): long stretches with almost no heterozygous
 * calls, from the binned tracks the VCF parser keeps.
 *
 * gVCF (callable track known): the genome is cut into 100 kb windows; each
 * window with at least half its length callable gets a het rate (het calls
 * per callable kb). A window is "low" when its rate is under 15% of the
 * median rate of the autosomes. Runs of low windows join into a segment,
 * bridging up to two windows that are not callable enough (but not a stretch
 * of more than 300 kb without calls, such as a centromere) and one window of
 * moderate rate (under 40% of the median), so a single error does not break
 * a run.
 * Plain VCF (no reference blocks; exomes): 1 Mb windows with at least 8
 * called variants; a window is low when under 10% of its calls are het, and
 * a run needs at least two low windows (2 Mb).
 *
 * Segments of 1 Mb and more are kept. F_ROH is the autosomal length in
 * segments of 1.5 Mb and more over the autosomal callable length (gVCF) or
 * the autosomal length (plain VCF). Long runs (over 5 Mb) point to recent
 * shared ancestry; one chromosome largely in ROH can mean uniparental
 * isodisomy. chrX in a male is hemizygous, so X and Y are listed but left
 * out of F_ROH.
 */
(function (G) {
  var MIN_SEG = 1e6, F_SEG = 1.5e6, LONG = 5e6;

  function isAutosome(name) { return /^(chr)?\d+$/i.test(name); }

  // Sum a track's base level over [b0, b1) bins.
  function sum(track, b0, b1) {
    var a = track.levels[0], s = 0;
    for (var b = Math.max(0, b0); b < Math.min(a.length, b1); b++) s += a[b];
    return s;
  }

  // Per contig: windows of W bp with het count, called variants and callable bp.
  function windows(data, W) {
    var out = [];
    data.genome.contigs.forEach(function (c) {
      var t = data.tracks && data.tracks[c.key];
      if (!t) return;
      var bs = t.het.binSize, per = Math.max(1, Math.round(W / bs)), n = Math.ceil(c.length / (per * bs)), list = [];
      for (var i = 0; i < n; i++) {
        var b0 = i * per, b1 = b0 + per, het = sum(t.het, b0, b1), hom = sum(t.hom, b0, b1);
        list.push({ start: b0 * bs + 1, end: Math.min(c.length, b1 * bs), het: het, called: het + hom, callable: sum(t.callable, b0, b1) * bs });
      }
      out.push({ contig: c, list: list });
    });
    return out;
  }

  function median(a) { if (!a.length) return 0; a = a.slice().sort(function (x, y) { return x - y; }); return a[a.length >> 1]; }

  function call(data) {
    if (!data || !data.tracks || data.format !== 'vcf') return null;
    var gvcf = !!data.isGvcf, W = gvcf ? 1e5 : 1e6, wins = windows(data, W);
    var rates = [];
    if (gvcf) wins.forEach(function (w) {
      if (!isAutosome(w.contig.name)) return;
      w.list.forEach(function (x) { if (x.callable >= W / 2) rates.push(x.het / (x.callable / 1000)); });
    });
    var med = median(rates);
    if (gvcf && !med) return null;
    // classify: 'low', 'mid' (moderate, may sit inside a run once), 'high', or null (cannot tell)
    var cls = function (x) {
      if (gvcf) {
        if (x.callable < W / 2) return null;
        var r = x.het / (x.callable / 1000);
        return r < 0.15 * med ? 'low' : r < 0.4 * med ? 'mid' : 'high';
      }
      if (x.called < 8) return null;
      var f = x.het / x.called;
      return f < 0.1 ? 'low' : f < 0.2 ? 'mid' : 'high';
    };
    var segs = [];
    wins.forEach(function (w) {
      var run = null, gaps = 0, mids = 0, lastLow = null;
      var close = function () {
        // plain VCF windows are 1 Mb and sparse (exomes): a run needs two low windows
        if (run && lastLow) { run.end = lastLow.end; if (run.end - run.start + 1 >= (gvcf ? MIN_SEG : 2 * W)) segs.push(run); }
        run = null; gaps = 0; mids = 0; lastLow = null;
      };
      w.list.forEach(function (x) {
        var k = cls(x);
        if (k === 'low') {
          if (!run) run = { key: w.contig.key, chrom: w.contig.name, start: x.start, end: x.end, het: 0, callable: 0 };
          run.het += x.het; run.callable += x.callable; lastLow = x; gaps = 0;
        } else if (!run) return;
        else if (k === null) { if (++gaps > 2 || x.end - lastLow.end > 3e5) close(); }
        else if (k === 'mid' && mids < 1) { mids++; run.het += x.het; run.callable += x.callable; }
        else close();
      });
      close();
    });
    segs.forEach(function (s) { s.length = s.end - s.start + 1; s.long = s.length > LONG; });
    // F_ROH over the autosomes
    var autoLen = 0, inRoh = 0, perChrom = {};
    wins.forEach(function (w) {
      if (!isAutosome(w.contig.name)) return;
      autoLen += gvcf ? w.list.reduce(function (s, x) { return s + x.callable; }, 0) : w.contig.length;
    });
    segs.forEach(function (s) {
      perChrom[s.chrom] = (perChrom[s.chrom] || 0) + s.length;
      if (isAutosome(s.chrom) && s.length >= F_SEG) inRoh += s.length;
    });
    // chromosomes mostly in ROH (possible isodisomy)
    var mostly = data.genome.contigs.filter(function (c) { return isAutosome(c.name) && (perChrom[c.name] || 0) > 0.5 * c.length; }).map(function (c) { return c.name; });
    return {
      method: gvcf ? 'gvcf' : 'vcf', window: W, medianRate: med, segments: segs,
      fRoh: autoLen ? inRoh / autoLen : 0, autosomalInRoh: inRoh, longCount: segs.filter(function (s) { return s.long && isAutosome(s.chrom); }).length,
      perChrom: perChrom, mostlyRoh: mostly
    };
  }

  // Segments overlapping [start, end] on contig key.
  function overlaps(res, key, start, end) {
    return res ? res.segments.filter(function (s) { return s.key === key && s.start <= end && s.end >= start; }) : [];
  }

  G.roh = { call: call, overlaps: overlaps, isAutosome: isAutosome };
})(globalThis.G = globalThis.G || {});
