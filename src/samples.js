/*
 * Several samples and ancestry.
 *
 *  - concordance: for each pair of samples, the share of sites called in
 *    both where the genotype class (0/0, het, hom) agrees.
 *  - trio: with a child and both parents, Mendelian inconsistencies (the
 *    child's genotype cannot come from the parents') and de novo candidates
 *    (child het, both parents a confident 0/0). Sites missing in any of the
 *    three are skipped: a missing call is never a 0/0 (Asclepius D4).
 *  - ancestry: port of Asclepius app/ancestry.py (genotype_loglik,
 *    ancestry_posterior): Hardy-Weinberg likelihood over its ancestry
 *    informative marker panel, uniform prior, softmax. Illustrative only,
 *    as Asclepius says: about 20 markers, independence assumed, no admixture.
 */
(function (G) {
  function called(z, Z) { return z === Z.HET || z === Z.HOM || z === Z.REF; }

  function concordance(d) {
    var n = Math.min(d.samples.length, 16), Z = G.vcf.Z, same = [], both = [];
    for (var a = 0; a < n; a++) { same.push(new Array(n).fill(0)); both.push(new Array(n).fill(0)); }
    Object.keys(d.variants).forEach(function (k) {
      var c = d.variants[k];
      for (var i = 0; i < c.n; i++) for (var a = 0; a < n; a++) {
        var za = c.zygOf(i, a);
        if (!called(za, Z)) continue;
        for (var b = a + 1; b < n; b++) {
          var zb = c.zygOf(i, b);
          if (!called(zb, Z)) continue;
          both[a][b]++; if (za === zb) same[a][b]++;
        }
      }
    });
    return { n: n, same: same, both: both, rate: function (a, b) { var i = Math.min(a, b), j = Math.max(a, b); return both[i][j] ? same[i][j] / both[i][j] : null; } };
  }

  // Can the child's class arise from the parents' classes (biallelic)?
  function mendelOk(zc, zf, zm, Z) {
    var alt = function (z) { return z === Z.HET || z === Z.HOM; }, ref = function (z) { return z === Z.HET || z === Z.REF; };
    if (zc === Z.HOM) return alt(zf) && alt(zm);
    if (zc === Z.REF) return ref(zf) && ref(zm);
    return (alt(zf) && ref(zm)) || (ref(zf) && alt(zm)); // HET
  }

  function trio(d, child, father, mother) {
    var Z = G.vcf.Z, out = { checked: 0, errors: 0, deNovo: [], errorSites: [] };
    Object.keys(d.variants).forEach(function (k) {
      var c = d.variants[k];
      for (var i = 0; i < c.n; i++) {
        var zc = c.zygOf(i, child), zf = c.zygOf(i, father), zm = c.zygOf(i, mother);
        if (!called(zc, Z) || !called(zf, Z) || !called(zm, Z)) continue;
        out.checked++;
        if (!mendelOk(zc, zf, zm, Z)) {
          out.errors++;
          if (out.errorSites.length < 5000) out.errorSites.push({ key: k, i: i });
          if (zc === Z.HET && zf === Z.REF && zm === Z.REF && c.pass[i] && out.deNovo.length < 2000) out.deNovo.push({ key: k, i: i, pos: c.pos[i] });
        }
      }
    });
    return out;
  }

  // ---- ancestry (Asclepius genotype_loglik / ancestry_posterior)

  function genotypeLoglik(dosage, freq) {
    var p = Math.min(Math.max(+freq, 1e-4), 1 - 1e-4);
    if (dosage === 0) return Math.log((1 - p) * (1 - p));
    if (dosage === 1) return Math.log(2 * p * (1 - p));
    if (dosage === 2) return Math.log(p * p);
    return 0;
  }
  function posterior(ll, pops) {
    var keys = Object.keys(ll);
    if (!keys.length) { var u = {}; pops.forEach(function (p) { u[p] = 1 / pops.length; }); return u; }
    var m = Math.max.apply(null, keys.map(function (k) { return ll[k]; })), e = {}, t = 0;
    keys.forEach(function (k) { e[k] = Math.exp(ll[k] - m); t += e[k]; });
    keys.forEach(function (k) { e[k] /= t; });
    return e;
  }

  // Alt dosage at one marker, or null when the site is not confidently read.
  function dosageAt(d, m) {
    var key = G.genome.normName(m.chrom), c = d.variants[key], Z = G.vcf.Z;
    if (!c) return null;
    var lo = 0, hi = c.n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (c.pos[mid] < m.pos) lo = mid + 1; else hi = mid; }
    var want = G.clinvar.variantKey(m.chrom, m.pos, m.ref, m.alt);
    for (var i = lo; i < c.n && c.pos[i] <= m.pos + 5; i++) {
      var al = c.alleles(i);
      if (!al) continue;
      for (var k = 0; k < al.alts.length; k++) {
        if (G.clinvar.variantKey(m.chrom, c.pos[i], al.ref, al.alts[k]) !== want) continue;
        var z = c.zyg[i];
        return z === Z.HOM ? 2 : z === Z.HET ? 1 : z === Z.REF ? 0 : null;
      }
    }
    // no record of this allele: a confident 0 only where a gVCF says the base was called
    if (!d.isGvcf) return null;
    var tr = d.tracks[key], bs = tr.callable.binSize, b = Math.floor((m.pos - 1) / bs);
    return tr.callable.levels[0][b] >= 0.9 && !(tr.lowdp.levels[0][b] > 0.1) ? 0 : null;
  }

  function ancestry(d, panel) {
    var pops = panel.superpops, ll = {}, used = [];
    pops.forEach(function (p) { ll[p] = 0; });
    panel.markers.forEach(function (m) {
      var dz = dosageAt(d, m);
      if (dz === null) return;
      used.push({ rsid: m.rsid, gene: m.gene, dosage: dz });
      pops.forEach(function (p) { ll[p] += genotypeLoglik(dz, m.freq[p]); });
    });
    return { posterior: posterior(used.length ? ll : {}, pops), used: used, total: panel.markers.length };
  }

  G.samples = { concordance: concordance, trio: trio, mendelOk: mendelOk, ancestry: ancestry, genotypeLoglik: genotypeLoglik, posterior: posterior, dosageAt: dosageAt };
})(globalThis.G = globalThis.G || {});
