/*
 * GWAS Catalog layer: genome-wide significant associations (p <= 5e-8,
 * NHGRI-EBI GWAS Catalog, GRCh38), one entry per SNP and risk allele,
 * exported by tools/fetch_annotations.py.
 *
 * For each SNP the page works out how many copies of the risk allele the
 * sample carries, from the VCF record at that position:
 *   risk allele = ALT -> 1 (het) or 2 (hom); risk allele = REF -> 2 - ALT dosage.
 *   A risk allele matching neither is tried on the other strand, unless the
 *   SNP is A/T or C/G (strand ambiguous; refused, as Asclepius prs.py does).
 *   No record at the site: the sample matches the reference there, but the
 *   catalog gives no reference base, so the count stays unknown.
 * These counts are not a risk score. They say which reported risk alleles
 * this genome carries.
 */
(function (G) {
  var COMP = { A: 'T', T: 'A', C: 'G', G: 'C' };

  function lowerBound(arr, n, v) {
    var lo = 0, hi = n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  function Gwas() { this.byContig = {}; this.traits = []; this.n = 0; this.meta = {}; }

  Gwas.prototype.load = async function (snpBlob, traitsJson) {
    var own = G.bgzf.own, tmp = {};
    for await (var line of G.bgzf.lines(snpBlob)) {
      if (line[0] === '#') { var eq = line.indexOf('='); if (eq > 0) this.meta[line.slice(1, eq)] = own(line.slice(eq + 1)); continue; }
      var f = line.split('\t');
      if (f.length < 8) continue;
      var key = G.genome.normName(f[0]), t = tmp[key] || (tmp[key] = { pos: [], risk: [], lp: [], rsid: [], gene: [], assoc: [], freq: [] });
      t.pos.push(+f[1]); t.risk.push(own(f[3])); t.lp.push(-Math.log10(Math.max(1e-300, +f[5] || 1)));
      t.rsid.push(own(f[2])); t.gene.push(own(f[6])); t.assoc.push(own(f[7])); t.freq.push(f[4] === 'NR' ? null : +f[4]);
      this.n++;
    }
    var self = this;
    Object.keys(tmp).forEach(function (k) {
      var t = tmp[k];
      self.byContig[k] = { n: t.pos.length, pos: Int32Array.from(t.pos), lp: Float32Array.from(t.lp), risk: t.risk, rsid: t.rsid, gene: t.gene, assoc: t.assoc, freq: t.freq };
    });
    this.traits = traitsJson.traits;
    return this;
  };

  // Associations of SNP i on contig k: [{trait, mapped, p, effect, study}], best first.
  Gwas.prototype.associations = function (k, i) {
    var self = this;
    return this.byContig[k].assoc[i].split(';').map(function (a) {
      var f = a.split(':'), tr = self.traits[+f[0]] || ['?', ''];
      return { traitIdx: +f[0], trait: tr[0], mapped: tr[1], p: +f[1], effect: f[2], study: f[3] };
    });
  };

  // Copies of the risk allele in sample s (0 = first), or null when unknown.
  // Returns { dosage, reason }.
  Gwas.prototype.dosage = function (data, k, i, s) {
    var g = this.byContig[k], risk = g.risk[i], c = data && data.variants ? data.variants[k] : null, Z = G.vcf.Z;
    if (!/^[ACGT]$/.test(risk)) return { dosage: null, reason: 'risk allele not a single base' };
    if (!c) return { dosage: null, reason: 'no variants on this chromosome' };
    var j = lowerBound(c.pos, c.n, g.pos[i]);
    if (j >= c.n || c.pos[j] !== g.pos[i]) return { dosage: null, reason: 'matches the reference here; reference base unknown' };
    var al = c.alleles(j);
    if (!al || al.ref.length !== 1) return { dosage: null, reason: 'not a simple SNV in the sample' };
    var z = c.zygOf(j, s || 0), altDose = z === Z.HOM ? 2 : z === Z.HET ? 1 : z === Z.REF ? 0 : null;
    if (altDose === null) return { dosage: null, reason: 'no confident genotype' };
    var alt = al.alts[0];
    var pal = COMP[al.ref] === alt;
    var pick = function (r) { return r === alt ? altDose : r === al.ref ? 2 - altDose : null; };
    var d = pick(risk);
    if (d !== null) return { dosage: d, reason: pal ? 'A/T or C/G SNP: strand taken as reported' : '' };
    if (pal) return { dosage: null, reason: 'A/T or C/G SNP, strand ambiguous' };
    d = pick(COMP[risk]);
    return d === null ? { dosage: null, reason: 'risk allele is neither REF nor ALT' } : { dosage: d, reason: 'risk allele matched on the other strand' };
  };

  // SNPs in [a, b] on contig k: indexes.
  Gwas.prototype.inRange = function (k, a, b) {
    var g = this.byContig[k];
    if (!g) return [];
    var out = [];
    for (var i = lowerBound(g.pos, g.n, a); i < g.n && g.pos[i] <= b; i++) out.push(i);
    return out;
  };

  // Traits whose text contains q (case-insensitive), with their SNP counts.
  Gwas.prototype.searchTraits = function (q, limit) {
    q = q.trim().toLowerCase();
    if (q.length < 3) return [];
    var hits = [];
    this.traits.forEach(function (t, i) { if (t[0].toLowerCase().indexOf(q) >= 0 || (t[1] || '').toLowerCase().indexOf(q) >= 0) hits.push(i); });
    var set = new Set(hits), count = {}, self = this;
    Object.keys(this.byContig).forEach(function (k) {
      self.byContig[k].assoc.forEach(function (a) {
        a.split(';').forEach(function (x) { var ti = +x.slice(0, x.indexOf(':')); if (set.has(ti)) count[ti] = (count[ti] || 0) + 1; });
      });
    });
    return hits.map(function (i) { return { idx: i, trait: self.traits[i][0], n: count[i] || 0 }; })
      .filter(function (t) { return t.n > 0; }).sort(function (a, b) { return b.n - a.n; }).slice(0, limit || 30);
  };

  // All SNPs associated with any of the trait indexes: [{k, i, p}].
  Gwas.prototype.lociFor = function (traitIdxs) {
    var set = new Set(traitIdxs), out = [], self = this;
    Object.keys(this.byContig).forEach(function (k) {
      var g = self.byContig[k];
      for (var i = 0; i < g.n; i++) {
        var parts = g.assoc[i].split(';'), best = null;
        for (var q = 0; q < parts.length; q++) {
          var f = parts[q].split(':');
          if (set.has(+f[0]) && (best === null || +f[1] < best)) best = +f[1];
        }
        if (best !== null) out.push({ k: k, i: i, p: best });
      }
    });
    return out.sort(function (a, b) { return a.p - b.p; });
  };

  G.Gwas = Gwas;
})(globalThis.G = globalThis.G || {});
