/*
 * ClinVar P/LP layer: loads the table exported from Asclepius
 * (tools/export_clinvar.py) and matches a sample's variants against it.
 *
 * The rules are ports of Asclepius, kept line for line so the two agree:
 *   normaliseVariant  <- app/match.py normalise_variant (D2)
 *   acmgClass         <- app/acmg.py acmg_class (order matters)
 *   classifySite      <- app/genotype.py classify_site (D4: missing != 0/0)
 *   countAllele       <- app/genotype.py count_allele
 *   classifyZygosity  <- app/vcf.py classify_zygosity
 *   reported          <- app/scan.py scan_variant plus ui.py run_clinical_scan
 * One deliberate difference: Asclepius only looks at the first ALT allele of
 * a record. Here every ALT allele is checked, with its own dosage.
 */
(function (G) {
  var TIERS = ['Pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign'];

  function normaliseVariant(chrom, pos, ref, alt) {
    ref = (ref || '').toUpperCase(); alt = (alt || '').toUpperCase(); pos = +pos;
    while (ref.length > 1 && alt.length > 1 && ref[ref.length - 1] === alt[alt.length - 1]) {
      ref = ref.slice(0, -1); alt = alt.slice(0, -1);
    }
    while (ref.length > 1 && alt.length > 1 && ref[0] === alt[0]) {
      ref = ref.slice(1); alt = alt.slice(1); pos++;
    }
    return { chrom: G.genome.normName(chrom), pos: pos, ref: ref, alt: alt };
  }
  function variantKey(chrom, pos, ref, alt) {
    var v = normaliseVariant(chrom, pos, ref, alt);
    return v.chrom + ':' + v.pos + ':' + v.ref + ':' + v.alt;
  }

  function acmgClass(sig) {
    var s = (sig || '').trim().toLowerCase();
    if (s.indexOf('likely pathogenic') >= 0) return 'Likely pathogenic';
    if (s.indexOf('pathogenic') >= 0) return 'Pathogenic';
    if (s.indexOf('likely benign') >= 0) return 'Likely benign';
    if (s.indexOf('benign') >= 0) return 'Benign';
    return 'Uncertain significance';
  }
  function isReportable(sig) { var c = acmgClass(sig); return c === 'Pathogenic' || c === 'Likely pathogenic'; }

  function alleles(gt) { return gt.trim().split(/[\/|]/); }

  function classifySite(gt, filter) {
    if (gt == null) return 'missing';
    var parts = alleles(gt);
    if (gt.trim() === '' || gt.trim() === '.' || parts.every(function (p) { return p === '.'; })) return 'missing';
    if (filter && filter.indexOf('LowDP') >= 0) return 'lowdp';
    if (parts.every(function (p) { return p === '0'; })) return 'reference';
    return 'variant';
  }

  function countAllele(gt, index) {
    if (gt == null) return null;
    var parts = alleles(gt);
    if (!parts.length || parts.some(function (p) { return p === '.' || p === ''; })) return null;
    var n = 0;
    for (var i = 0; i < parts.length; i++) {
      if (!/^-?\d+$/.test(parts[i].trim())) return null;
      if (parseInt(parts[i], 10) === index) n++;
    }
    return n;
  }

  function classifyZygosity(gt) {
    var g = (gt || '').split(':')[0].trim();
    var al = g.split(/[\/|]/);
    if (al.length < 2 || al.some(function (a) { return a === '' || a === '.'; })) return 'Unknown';
    var alt = al.filter(function (a) { return a !== '0'; }).length;
    if (alt === 0) return 'Reference';
    if (alt === al.length) return 'Homozygous';
    return 'Carrier';
  }

  // Why a ClinVar allele present in the file is or is not a reported finding.
  // Mirrors run_clinical_scan: state must be 'variant', dosage > 0, then
  // scan_variant needs a reportable class and a zygosity that is not
  // Reference or Unknown.
  function judge(entry, gt, filter, alleleIndex) {
    var state = classifySite(gt, filter);
    if (state !== 'variant') return { reported: false, reason: state === 'lowdp' ? 'low depth call (LowDP)' : state === 'missing' ? 'no call (missing genotype)' : 'genotype is reference (0/0)' };
    var dosage = countAllele(gt, alleleIndex);
    if (!dosage) return { reported: false, reason: dosage === null ? 'dosage unreadable' : 'sample carries a different ALT at this site' };
    if (!isReportable(entry.sig)) return { reported: false, reason: 'not P/LP' };
    var z = classifyZygosity(gt);
    if (z === 'Reference' || z === 'Unknown') return { reported: false, reason: 'zygosity ' + z + ' (e.g. haploid call)' };
    return { reported: true, zygosity: z, dosage: dosage };
  }

  function ClinVar() {
    this.index = new Map(); this.byContig = {}; this.meta = {}; this.n = 0; this.tracks = {};
  }

  ClinVar.prototype.load = async function (blob) {
    var byContig = {};
    for await (var line of G.bgzf.lines(blob)) {
      if (line[0] === '#') {
        var eq = line.indexOf('=');
        if (eq > 0) this.meta[line.slice(1, eq)] = line.slice(eq + 1);
        continue;
      }
      var f = line.split('\t');
      if (f.length < 6 || !f[2] || !f[3]) continue; // rows without alleles cannot be matched (D2)
      var own = G.bgzf.own;
      var e = { chrom: own(f[0]), pos: +f[1], ref: own(f[2]), alt: own(f[3]), gene: own(f[4]), sig: own(f[5]), name: own(f[6] || ''), pheno: own(f[7] || '') };
      e.tier = acmgClass(e.sig);
      var k = variantKey(e.chrom, e.pos, e.ref, e.alt);
      var list = this.index.get(k);
      if (list) list.push(e); else this.index.set(k, [e]);
      var ck = G.genome.normName(e.chrom);
      (byContig[ck] = byContig[ck] || []).push(e.pos);
      this.n++;
    }
    for (var c in byContig) this.byContig[c] = Int32Array.from(byContig[c]).sort();
    this.build = this.meta.build || 'GRCh38';
    return this;
  };

  // Called by the VCF parser for each record. Returns hits (reported or not).
  ClinVar.prototype.matchRecord = function (chrom, pos, ref, altField, gt, filter) {
    var hits = null, alts = altField.split(',');
    for (var i = 0; i < alts.length; i++) {
      var a = alts[i];
      if (a[0] === '<' || a === '*' || a === '.' || /[\[\]]/.test(a)) continue;
      var list = this.index.get(variantKey(chrom, pos, ref, a));
      if (!list) continue;
      for (var j = 0; j < list.length; j++) {
        var e = list[j], verdict = judge(e, gt, filter, i + 1);
        var own = G.bgzf.own;
        (hits = hits || []).push({
          chrom: own(chrom), pos: pos, ref: own(ref), alt: own(a), alleleIndex: i + 1, gt: own(gt), filter: own(filter),
          gene: e.gene, significance: e.sig, classification: e.tier, variant_name: e.name, phenotype: e.pheno,
          clinvar: { chrom: e.chrom, pos: e.pos, ref: e.ref, alt: e.alt },
          reported: verdict.reported, zygosity: verdict.zygosity || classifyZygosity(gt), reason: verdict.reason || null
        });
      }
    }
    return hits;
  };

  // ClinVar site density on the sample's contig, binned like its other tracks.
  ClinVar.prototype.trackFor = function (key, length, binSize) {
    var id = key + ':' + binSize;
    if (this.tracks[id]) return this.tracks[id];
    var pos = this.byContig[key];
    if (!pos) return null;
    var t = new G.genome.Track(length, binSize);
    for (var i = 0; i < pos.length; i++) t.add(pos[i], 1);
    t.buildPyramid();
    this.tracks[id] = t;
    return t;
  };

  // How many ClinVar sites fall where a gVCF has no confident call.
  // Approximate: uses the callable track at bin resolution.
  ClinVar.prototype.callability = function (data) {
    var out = { total: 0, called: 0, lowdp: 0, nocall: 0 }, self = this;
    if (!data.isGvcf) return null;
    data.genome.contigs.forEach(function (c) {
      var pos = self.byContig[c.key], tr = data.tracks[c.key];
      if (!pos || !tr) return;
      var cal = tr.callable.levels[0], low = tr.lowdp.levels[0], bs = tr.callable.binSize;
      for (var i = 0; i < pos.length; i++) {
        var b = Math.floor((pos[i] - 1) / bs);
        var cv = b < cal.length ? cal[b] : 0, lv = b < low.length ? low[b] : 0;
        out.total++;
        if (cv >= 0.5) out.called++; else if (lv >= 0.5) out.lowdp++; else out.nocall++;
      }
    });
    return out;
  };

  // Findings contract v1: what one producer says about one sample.
  // See docs/findings-contract.md.
  function toFindings(hits, producer, sample, build, reference) {
    return {
      schema: 'genomeatrium.findings/1',
      producer: producer, reference: reference, build: build, sample: sample,
      created: new Date().toISOString(),
      findings: hits.map(function (h) {
        return {
          chrom: h.clinvar ? h.clinvar.chrom : h.chrom, pos: h.clinvar ? h.clinvar.pos : h.pos,
          ref: h.clinvar ? h.clinvar.ref : h.ref, alt: h.clinvar ? h.clinvar.alt : h.alt,
          gene: h.gene, classification: h.classification, significance: h.significance,
          zygosity: h.zygosity, gt: h.gt, status: h.reported ? 'reported' : 'not_reported',
          reason: h.reason || undefined, variant_name: h.variant_name, phenotype: h.phenotype
        };
      })
    };
  }

  // Reads a findings file and checks the parts the view relies on.
  function parseFindings(text) {
    var d = JSON.parse(text);
    // moebiotobio.findings/1 is the same contract under the project's earlier name
    if (!d || typeof d.schema !== 'string' || !/^(genomeatrium|moebiotobio)\.findings\//.test(d.schema)) throw new Error('Not a findings file (schema genomeatrium.findings/1 expected)');
    if (d.schema !== 'genomeatrium.findings/1' && d.schema !== 'moebiotobio.findings/1') throw new Error('Unsupported findings schema ' + d.schema);
    if (!Array.isArray(d.findings)) throw new Error('findings must be a list');
    d.findings.forEach(function (f, i) {
      if (!f.chrom || !(f.pos > 0) || !f.gene) throw new Error('finding ' + i + ' lacks chrom, pos or gene');
      f.key = variantKey(f.chrom, f.pos, f.ref, f.alt);
    });
    return d;
  }

  G.clinvar = {
    ClinVar: ClinVar, TIERS: TIERS, normaliseVariant: normaliseVariant, variantKey: variantKey,
    acmgClass: acmgClass, isReportable: isReportable, classifySite: classifySite, countAllele: countAllele,
    classifyZygosity: classifyZygosity, judge: judge, toFindings: toFindings, parseFindings: parseFindings
  };
})(globalThis.G = globalThis.G || {});
