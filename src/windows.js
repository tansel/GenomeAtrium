/*
 * Genome windows: the shared model behind the Landscape, the Matrix and the
 * similarity arcs. The genome is cut into about 2,500 windows, each window
 * gets a small feature vector, and features are standardised so windows can
 * be compared.
 *
 * Similarity is the cosine of two standardised feature vectors, so 1 means
 * "same profile" (variant rate, het fraction, filtering, callability), not
 * physical contact. An ordinary window sits near the average profile, where
 * its direction is mostly noise, so the link strength also weights by how
 * far both windows are from average:
 *   strength = cosine * min(1, min(|z_i|, |z_j|) / |z| at the 80th percentile)
 * Strength plays the part of attention scores in moebio.com/attention:
 * unusual windows attend to the windows most like them.
 */
(function (G) {
  var TARGET_WINDOWS = 2500;

  function fmtBp(n) {
    return n >= 1e6 ? (n / 1e6).toFixed(1) + ' Mb' : n >= 1e3 ? (n / 1e3).toFixed(0) + ' kb' : n + ' bp';
  }

  function niceWindow(total, binSize) {
    var raw = Math.max(binSize, total / TARGET_WINDOWS);
    var p = Math.pow(10, Math.floor(Math.log10(raw)));
    var w = [1, 2, 5, 10].map(function (f) { return f * p; }).find(function (v) { return v >= raw; });
    return Math.max(binSize, Math.ceil(w / binSize) * binSize);
  }

  // Sums (or averages) base-level bins into windows of `per` bins.
  function windowed(track, nWin, per, kind) {
    var a = track.levels[0], out = new Float64Array(nWin);
    for (var i = 0; i < a.length; i++) { var w = Math.floor(i / per); if (w < nWin) out[w] += a[i]; }
    if (kind === 'mean') for (var j = 0; j < nWin; j++) out[j] /= per;
    return out;
  }

  // Builds windows and their features from a parsed VCF or BAM result.
  // SV and read-pair arcs are left out of the features on purpose: the
  // similarity arcs should say something the SV arcs do not already say.
  function buildWindows(d) {
    var contigs = d.genome.contigs, total = d.genome.totalLength();
    var win = niceWindow(total, d.binSize), per = Math.round(win / d.binSize);
    var names, windows = [];

    contigs.forEach(function (c, ci) {
      var tr = d.tracks[c.key];
      if (!tr) return;
      var nWin = Math.ceil(c.length / win);
      if (d.format === 'vcf') {
        names = ['SNV per kb', 'indel per kb', 'het fraction', 'filtered fraction'];
        var snv = windowed(tr.snv, nWin, per), indel = windowed(tr.indel, nWin, per);
        var het = new Float64Array(nWin), hom = new Float64Array(nWin), fail = new Float64Array(nWin), all = new Float64Array(nWin);
        var v = d.variants[c.key], Z = G.vcf.Z;
        for (var i = 0; v && i < v.n; i++) {
          var w = Math.floor((v.pos[i] - 1) / win);
          if (w >= nWin) continue;
          all[w]++;
          if (v.zyg[i] === Z.HET) het[w]++; else if (v.zyg[i] === Z.HOM) hom[w]++;
          if (!v.pass[i]) fail[w]++;
        }
        var cal = d.isGvcf ? windowed(tr.callable, nWin, per, 'mean') : null;
        if (cal) names.push('callable fraction');
        for (var j = 0; j < nWin; j++) {
          if (!all[j] && !(cal && cal[j] > 0.01)) continue; // no data: gap, not a point
          var len = Math.min(win, c.length - j * win) / 1000;
          var f = [Math.log1p(snv[j] / len), Math.log1p(indel[j] / len),
            het[j] + hom[j] ? het[j] / (het[j] + hom[j]) : 0.5, all[j] ? fail[j] / all[j] : 0];
          if (cal) f.push(Math.min(1, cal[j]));
          windows.push({ ci: ci, contig: c, j: j, start: j * win + 1, end: Math.min(c.length, (j + 1) * win), f: f, n: all[j] });
        }
      } else {
        names = ['log mean depth', 'depth variation', 'zero-depth fraction'];
        var dp = tr.depth.levels[0];
        for (var k = 0; k < nWin; k++) {
          var s = 0, s2 = 0, zero = 0, m = 0;
          for (var b = k * per; b < Math.min(dp.length, (k + 1) * per); b++) { s += dp[b]; s2 += dp[b] * dp[b]; if (dp[b] < 0.01) zero++; m++; }
          if (!m || s === 0) continue;
          var mean = s / m, sd = Math.sqrt(Math.max(0, s2 / m - mean * mean));
          windows.push({ ci: ci, contig: c, j: k, start: k * win + 1, end: Math.min(c.length, (k + 1) * win),
            f: [Math.log1p(mean), mean ? sd / mean : 0, zero / m], n: mean });
        }
      }
    });
    windows.forEach(function (w, i) { w.index = i; });
    return { windows: windows, names: names || [], win: win };
  }

  // z-scores each feature, then scales every window's vector to unit length,
  // so a dot product is a cosine similarity.
  function standardise(windows, names) {
    var n = windows.length, dim = names.length;
    var mean = new Array(dim).fill(0), sd = new Array(dim).fill(0);
    windows.forEach(function (w) { for (var i = 0; i < dim; i++) mean[i] += w.f[i] / n; });
    windows.forEach(function (w) { for (var i = 0; i < dim; i++) sd[i] += Math.pow(w.f[i] - mean[i], 2) / n; });
    sd = sd.map(Math.sqrt);
    var used = names.map(function (_, i) { return i; }).filter(function (i) { return sd[i] > 1e-9; });
    var m = used.length, unit = new Float32Array(n * m), norms = new Float32Array(n), Z = [];
    windows.forEach(function (w, wi) {
      var z = used.map(function (i) { return (w.f[i] - mean[i]) / sd[i]; });
      Z.push(z);
      var norm = Math.sqrt(z.reduce(function (s, v) { return s + v * v; }, 0));
      norms[wi] = norm;
      for (var k = 0; k < m; k++) unit[wi * m + k] = norm ? z[k] / norm : 0;
    });
    var sorted = Array.from(norms).sort(function (a, b) { return a - b; });
    var normRef = sorted[Math.floor(n * 0.8)] || 1;
    return { Z: Z, used: used, unit: unit, m: m, mean: mean, sd: sd, norms: norms, normRef: normRef };
  }

  function similarity(model, i, j) {
    var m = model.m, u = model.unit, s = 0;
    for (var k = 0; k < m; k++) s += u[i * m + k] * u[j * m + k];
    return s;
  }

  function strength(model, i, j) {
    var w = Math.min(1, Math.min(model.norms[i], model.norms[j]) / model.normRef);
    return similarity(model, i, j) * w;
  }

  // Strongest links, attention style: each window keeps its `perWindow` most
  // similar partners that are not its neighbours on the same chromosome;
  // the best `maxPairs` of those, above `minSim`, become arcs.
  function similarPairs(model, opts) {
    opts = opts || {};
    var perWindow = opts.perWindow || 2, maxPairs = opts.maxPairs || 3000, minSim = opts.minSim || 0.85;
    var gapWindows = opts.gapWindows || 5;
    var W = model.windows, n = W.length, seen = new Set(), pairs = [];
    for (var i = 0; i < n; i++) {
      var best = [];
      for (var j = 0; j < n; j++) {
        if (j === i || (W[j].ci === W[i].ci && Math.abs(W[j].j - W[i].j) <= gapWindows)) continue;
        var s = strength(model, i, j);
        if (s < minSim) continue;
        if (best.length < perWindow) { best.push([j, s]); best.sort(function (a, b) { return b[1] - a[1]; }); }
        else if (s > best[perWindow - 1][1]) { best[perWindow - 1] = [j, s]; best.sort(function (a, b) { return b[1] - a[1]; }); }
      }
      best.forEach(function (b) {
        var key = Math.min(i, b[0]) + ':' + Math.max(i, b[0]);
        if (seen.has(key)) return;
        seen.add(key);
        pairs.push({ i: Math.min(i, b[0]), j: Math.max(i, b[0]), s: b[1] });
      });
    }
    pairs.sort(function (a, b) { return b.s - a.s; });
    return pairs.slice(0, maxPairs);
  }

  function describe(model, w) {
    return model.names.map(function (nm, fi) { return nm + ' ' + w.f[fi].toFixed(2); }).join(', ');
  }

  // Full model for one parsed file, built once.
  function model(d) {
    if (!d || !d.genome || !d.genome.contigs.length) return null;
    if (d._windows !== undefined) return d._windows;
    var built = buildWindows(d);
    if (built.windows.length < 4) { d._windows = null; return null; }
    var st = standardise(built.windows, built.names);
    var mdl = { windows: built.windows, names: built.names, win: built.win, Z: st.Z, used: st.used, unit: st.unit, m: st.m,
      norms: st.norms, normRef: st.normRef };
    mdl.pairs = similarPairs(mdl);
    d._windows = mdl;
    return mdl;
  }

  // Similarity pairs as arcs the Arcs view can draw.
  function pairsToArcs(mdl) {
    return mdl.pairs.map(function (p) {
      var a = mdl.windows[p.i], b = mdl.windows[p.j];
      return {
        c0: a.contig.name, p0: Math.round((a.start + a.end) / 2), c1: b.contig.name, p1: Math.round((b.start + b.end) / 2),
        type: 'similar', support: p.s, pass: 1, wi: p.i, wj: p.j,
        label: 'alike windows (strength ' + p.s.toFixed(2) + ', cosine ' + similarity(mdl, p.i, p.j).toFixed(2) + '): ' + a.contig.name + ':' + fmtBp(a.start) + ' and ' + b.contig.name + ':' + fmtBp(b.start),
        detail: [a.contig.name + ':' + a.start.toLocaleString() + '  ' + describe(mdl, a), b.contig.name + ':' + b.start.toLocaleString() + '  ' + describe(mdl, b)]
      };
    });
  }

  G.windows = { model: model, buildWindows: buildWindows, standardise: standardise, similarity: similarity, strength: strength,
    similarPairs: similarPairs, pairsToArcs: pairsToArcs, niceWindow: niceWindow, describe: describe };
})(globalThis.G = globalThis.G || {});
