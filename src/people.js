/*
 * People: several persons on one genome (a family, or a few samples), each
 * with a colour. A person is one sample of a loaded VCF: a separate file per
 * person, or the samples of a joint (multi-sample) VCF, up to 8 at a time.
 * All must share the build (coordinates are never compared across builds).
 *
 * The first person is primary: the existing views draw that person as
 * before. The others are drawn alongside, in their colour: a het-fraction
 * row each in Arcs, their findings, their runs of homozygosity.
 *
 * Joint-VCF samples beyond the first have no binned tracks of their own in
 * the parser, so het and hom counts are built here from the kept genotypes,
 * in bins of 32 base bins (about 25 kb), and their runs of homozygosity use
 * the plain-VCF method (het share of called sites).
 *
 * Parent check (sharing): of the child's non-reference alleles, how many each
 * parent also carries (same CHROM, POS, REF, ALT; Asclepius D2). An allele
 * absent from a parent's file is "not seen", never "not carried" (D4).
 * Around half from each parent, and few seen in neither, is what a true
 * child of both looks like.
 */
(function (G) {
  var PALETTE = ['#5ad2be', '#ff8a5c', '#b48cff', '#ffd34d', '#5aa0ff', '#ff5c9a', '#9be15a', '#e0e0e0'];
  var MAX = 8;

  function Person(data, sample, name, index) {
    this.data = data; this.sample = sample || 0; this.name = name; this.index = index;
    this.color = PALETTE[index % PALETTE.length]; this.role = ''; this.visible = true;
  }

  // Binned het/hom tracks for one sample of a joint VCF (sample 0 has the parser's own).
  Person.prototype.tracks = function () {
    if (this.sample === 0) return this.data.tracks;
    if (this._tracks) return this._tracks;
    var d = this.data, Z = G.vcf.Z, s = this.sample, bs = d.binSize * 32, out = {};
    d.genome.contigs.forEach(function (c) {
      var col = d.variants[c.key];
      if (!col) return;
      var het = new G.genome.Track(c.length, bs), hom = new G.genome.Track(c.length, bs);
      for (var i = 0; i < col.n; i++) { var z = col.zygOf(i, s); if (z === Z.HET) het.add(col.pos[i], 1); else if (z === Z.HOM) hom.add(col.pos[i], 1); }
      het.buildPyramid(); hom.buildPyramid();
      out[c.key] = { het: het, hom: hom, callable: new G.genome.Track(c.length, bs, 'mean') };
    });
    this._tracks = out;
    return out;
  };

  // Runs of homozygosity for this person (roh.js).
  Person.prototype.roh = function () {
    if (this._roh !== undefined) return this._roh;
    var d = this.data;
    this._roh = this.sample === 0 ? G.roh.call(d) : G.roh.call({ format: 'vcf', isGvcf: false, genome: d.genome, binSize: d.binSize * 32, tracks: this.tracks() });
    return this._roh;
  };

  // ClinVar matches with this person's own genotype, as a findings list.
  Person.prototype.findings = function () {
    if (this._findings) return this._findings;
    var d = this.data, hits = this.sample === 0 ? d.clinvarHits : (d.clinvarHitsBySample || {})[this.sample];
    if (!hits || !hits.length) return (this._findings = []);
    var doc = G.clinvar.toFindings(hits, { name: 'GenomeAtrium', code: 'src/clinvar.js' }, { file: d.fileName, sample: this.name }, d.build, 'ClinVar P/LP from Asclepius');
    this._findings = doc.findings.filter(function (f) { return f.status === 'reported'; });
    return this._findings;
  };

  function lowerBound(arr, n, v) { var lo = 0, hi = n; while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; } return lo; }

  // Does person p carry the allele ref>alt at key:pos (het or hom)? true, false (called
  // without it), or null (not seen in p's file: unknown, never a 0/0).
  function carries(p, key, pos, ref, alt) {
    var c = p.data.variants && p.data.variants[key], Z = G.vcf.Z;
    if (!c) return null;
    for (var j = lowerBound(c.pos, c.n, pos); j < c.n && c.pos[j] === pos; j++) {
      var al = c.alleles(j);
      if (!al || al.ref !== ref || al.alts.indexOf(alt) < 0) continue;
      var z = c.zygOf(j, p.sample);
      if (z === Z.HET || z === Z.HOM) return al.alts.length === 1 ? true : null; // which alt is carried is unknown at multi-allelic sites
      if (z === Z.REF) return false;
      return null;
    }
    return null;
  }

  // Of the child's non-reference alleles (biallelic, PASS), how many each parent carries.
  function sharing(child, mother, father) {
    var Z = G.vcf.Z, out = { n: 0, mother: 0, father: 0, both: 0, neither: 0 };
    child.data.genome.contigs.forEach(function (c) {
      var col = child.data.variants && child.data.variants[c.key];
      if (!col) return;
      for (var i = 0; i < col.n; i++) {
        var z = col.zygOf(i, child.sample);
        if (z !== Z.HET && z !== Z.HOM) continue;
        var al = col.alleles(i);
        if (!al || al.alts.length !== 1) continue;
        var m = carries(mother, c.key, col.pos[i], al.ref, al.alts[0]), f = carries(father, c.key, col.pos[i], al.ref, al.alts[0]);
        out.n++;
        if (m && f) out.both++; else if (m) out.mother++; else if (f) out.father++; else out.neither++;
      }
    });
    return out;
  }

  // Autosomal stretches homozygous in two or more of the people, with who shares them.
  function sharedRoh(people) {
    var segs = [];
    people.forEach(function (p, pi) { var r = p.roh(); if (r) r.segments.forEach(function (s) { if (G.roh.isAutosome(s.chrom)) segs.push({ key: s.key, chrom: s.chrom, start: s.start, end: s.end, who: pi }); }); });
    var out = [];
    for (var a = 0; a < segs.length; a++) for (var b = a + 1; b < segs.length; b++) {
      var x = segs[a], y = segs[b];
      if (x.who === y.who || x.key !== y.key) continue;
      var lo = Math.max(x.start, y.start), hi = Math.min(x.end, y.end);
      if (hi - lo + 1 >= 1e6) out.push({ chrom: x.chrom, key: x.key, start: lo, end: hi, who: [x.who, y.who] });
    }
    return out.sort(function (a, b) { return (b.end - b.start) - (a.end - a.start); });
  }

  G.people = { PALETTE: PALETTE, MAX: MAX, Person: Person, carries: carries, sharing: sharing, sharedRoh: sharedRoh };
})(globalThis.G = globalThis.G || {});
