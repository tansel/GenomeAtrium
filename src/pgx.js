/*
 * Pharmacogenomics, ported from Asclepius app/pgx.py (KB 2026-06, reviewed 2026-06-14):
 * the CPIC-actionable genes with drug, recommendation, evidence level and
 * link, and the function-defining variants readable from one VCF position
 * (GRCh38). Generated from pgx.py (one em-dash in a recommendation became a
 * colon); its tests are in tests/pgx.test.js.
 *
 * As Asclepius decides (D7, partial lenses are honest and cited), this is
 * carrier detection per variant with gene-level guidance. It does not call
 * diplotypes or metaboliser status: that needs star-allele calling, and for
 * CYP2D6 copy number (deletions, duplications, CYP2D7 hybrids) beyond a VCF.
 *
 * genotypeAt reads one variant on the loaded genome as one of:
 *   copies      0, 1 or 2 copies of the variant allele, called
 *   reference   no record, and the site's bin fully callable in a gVCF: 0 copies
 *   other       a different allele is called at the site (CHROM+POS+REF+ALT, D2)
 *   lowdp / nocall   the site is low depth, or not called (gVCF)
 *   unknown     not in a plain VCF: unknown, never 0 copies (D4)
 */
(function (G) {
  var KB_VERSION = '2026-06', KB_REVIEWED = '2026-06-14', CPIC_LEVELS = ['A', 'B'];
  var PGX_GENES = {
    "DPYD": {
      "drug": "Fluoropyrimidines (5-FU, capecitabine)",
      "cpic_level": "A",
      "recommendation": "Reduce dose or select alternative per DPYD activity score; intensive monitoring.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-fluoropyrimidines-and-dpyd/"
    },
    "TPMT": {
      "drug": "Thiopurines (azathioprine, mercaptopurine)",
      "cpic_level": "A",
      "recommendation": "Reduce starting dose for intermediate/poor metabolizers to avoid myelosuppression.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/"
    },
    "NUDT15": {
      "drug": "Thiopurines (azathioprine, mercaptopurine)",
      "cpic_level": "A",
      "recommendation": "Reduce thiopurine dose for intermediate/poor metabolizers (esp. East-Asian ancestry).",
      "link": "https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/"
    },
    "CYP2C19": {
      "drug": "Clopidogrel, voriconazole, SSRIs",
      "cpic_level": "A",
      "recommendation": "Poor metabolizers: avoid clopidogrel (use prasugrel/ticagrelor); adjust voriconazole/SSRI dose.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-clopidogrel-and-cyp2c19/"
    },
    "CYP2D6": {
      "drug": "Codeine, tramadol, tamoxifen, atomoxetine",
      "cpic_level": "A",
      "recommendation": "Avoid codeine/tramadol in ultra-rapid & poor metabolizers. NOTE: needs star-allele + copy-number calling.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-codeine-and-cyp2d6/"
    },
    "CYP2C9": {
      "drug": "Warfarin, phenytoin, NSAIDs",
      "cpic_level": "A",
      "recommendation": "Combine with VKORC1 for warfarin dosing; reduce phenytoin dose in poor metabolizers.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/"
    },
    "VKORC1": {
      "drug": "Warfarin",
      "cpic_level": "A",
      "recommendation": "Use with CYP2C9 genotype in the CPIC/IWPC warfarin dosing algorithm.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/"
    },
    "SLCO1B1": {
      "drug": "Statins (simvastatin)",
      "cpic_level": "A",
      "recommendation": "Decreased-function: increased myopathy risk; lower dose or alternative statin.",
      "link": "https://cpicpgx.org/guidelines/cpic-guideline-for-statins/"
    },
    "UGT1A1": {
      "drug": "Irinotecan, atazanavir",
      "cpic_level": "A",
      "recommendation": "Reduced-function (*28/*28): increased irinotecan toxicity; consider dose reduction.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-irinotecan-and-ugt1a1/"
    },
    "CYP3A5": {
      "drug": "Tacrolimus",
      "cpic_level": "A",
      "recommendation": "Expressers need a higher starting tacrolimus dose to reach target troughs.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-cyp3a5/"
    },
    "G6PD": {
      "drug": "Rasburicase, primaquine, dapsone",
      "cpic_level": "A",
      "recommendation": "Deficient: AVOID oxidant drugs (haemolysis risk).",
      "link": "https://cpicpgx.org/guidelines/cpic-guideline-for-rasburicase-and-g6pd/"
    },
    "HLA-B": {
      "drug": "Abacavir (B*57:01), carbamazepine/oxcarbazepine (B*15:02)",
      "cpic_level": "A",
      "recommendation": "Allele carriers: AVOID the drug (abacavir hypersensitivity / SJS-TEN). Allele-specific: overlaps HLA Nexus typing.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/"
    },
    "HLA-A": {
      "drug": "Carbamazepine (A*31:01)",
      "cpic_level": "A",
      "recommendation": "A*31:01 carriers: increased risk of carbamazepine hypersensitivity; consider alternative.",
      "link": "https://cpicpgx.org/guidelines/guideline-for-carbamazepine-and-hla-b/"
    }
  };
  var PGX_VARIANTS = [
    {
      "rsid": "rs3918290",
      "gene": "DPYD",
      "star": "*2A",
      "chrom": "1",
      "pos": 97450058,
      "ref": "C",
      "alt": "A"
    },
    {
      "rsid": "rs55886062",
      "gene": "DPYD",
      "star": "*13",
      "chrom": "1",
      "pos": 97515787,
      "ref": "A",
      "alt": "C"
    },
    {
      "rsid": "rs67376798",
      "gene": "DPYD",
      "star": "c.2846A>T",
      "chrom": "1",
      "pos": 97082391,
      "ref": "T",
      "alt": "A"
    },
    {
      "rsid": "rs75017182",
      "gene": "DPYD",
      "star": "HapB3",
      "chrom": "1",
      "pos": 97579893,
      "ref": "G",
      "alt": "A"
    },
    {
      "rsid": "rs1800462",
      "gene": "TPMT",
      "star": "*2",
      "chrom": "6",
      "pos": 18143724,
      "ref": "C",
      "alt": "G"
    },
    {
      "rsid": "rs1800460",
      "gene": "TPMT",
      "star": "*3B",
      "chrom": "6",
      "pos": 18138997,
      "ref": "C",
      "alt": "A"
    },
    {
      "rsid": "rs1142345",
      "gene": "TPMT",
      "star": "*3C",
      "chrom": "6",
      "pos": 18130687,
      "ref": "T",
      "alt": "A"
    },
    {
      "rsid": "rs116855232",
      "gene": "NUDT15",
      "star": "*3",
      "chrom": "13",
      "pos": 48045719,
      "ref": "C",
      "alt": "T"
    },
    {
      "rsid": "rs4244285",
      "gene": "CYP2C19",
      "star": "*2",
      "chrom": "10",
      "pos": 94781859,
      "ref": "G",
      "alt": "A"
    },
    {
      "rsid": "rs4986893",
      "gene": "CYP2C19",
      "star": "*3",
      "chrom": "10",
      "pos": 94780653,
      "ref": "G",
      "alt": "A"
    },
    {
      "rsid": "rs12248560",
      "gene": "CYP2C19",
      "star": "*17 (incr.)",
      "chrom": "10",
      "pos": 94761900,
      "ref": "C",
      "alt": "A"
    },
    {
      "rsid": "rs1799853",
      "gene": "CYP2C9",
      "star": "*2",
      "chrom": "10",
      "pos": 94942290,
      "ref": "C",
      "alt": "A"
    },
    {
      "rsid": "rs1057910",
      "gene": "CYP2C9",
      "star": "*3",
      "chrom": "10",
      "pos": 94981296,
      "ref": "A",
      "alt": "C"
    },
    {
      "rsid": "rs9923231",
      "gene": "VKORC1",
      "star": "-1639A",
      "chrom": "16",
      "pos": 31096368,
      "ref": "C",
      "alt": "A"
    },
    {
      "rsid": "rs4149056",
      "gene": "SLCO1B1",
      "star": "*5",
      "chrom": "12",
      "pos": 21178615,
      "ref": "T",
      "alt": "A"
    },
    {
      "rsid": "rs887829",
      "gene": "UGT1A1",
      "star": "*80(tag*28)",
      "chrom": "2",
      "pos": 233759924,
      "ref": "C",
      "alt": "T"
    },
    {
      "rsid": "rs776746",
      "gene": "CYP3A5",
      "star": "*3",
      "chrom": "7",
      "pos": 99672916,
      "ref": "T",
      "alt": "C"
    },
    {
      "rsid": "rs1050828",
      "gene": "G6PD",
      "star": "A-202A",
      "chrom": "X",
      "pos": 154536002,
      "ref": "C",
      "alt": "T"
    },
    {
      "rsid": "rs3892097",
      "gene": "CYP2D6",
      "star": "*4",
      "chrom": "22",
      "pos": 42128945,
      "ref": "C",
      "alt": "A"
    }
  ];

  function isPgxGene(gene) { return Object.prototype.hasOwnProperty.call(PGX_GENES, String(gene || '').toUpperCase()); }
  function pgxGuidance(gene) {
    var g = String(gene || '').toUpperCase(), e = PGX_GENES[g];
    if (!e) return null;
    var out = { gene: g, standard: 'CPIC' };
    Object.keys(e).forEach(function (k) { out[k] = e[k]; });
    return out;
  }

  function lowerBound(arr, n, v) { var lo = 0, hi = n; while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; } return lo; }

  function genotypeAt(data, v, sample) {
    var Z = G.vcf.Z, c = data.genome.get(v.chrom);
    if (!c) return { state: 'unknown', why: 'chromosome not in this file' };
    var col = data.variants && data.variants[c.key];
    if (col) for (var j = lowerBound(col.pos, col.n, v.pos); j < col.n && col.pos[j] === v.pos; j++) {
      var al = col.alleles(j), z = col.zygOf(j, sample || 0);
      if (!al) continue;
      var hasAlt = al.ref.toUpperCase() === v.ref && al.alts.map(function (a) { return a.toUpperCase(); }).indexOf(v.alt) >= 0;
      if (!hasAlt) return { state: 'other', why: 'called ' + al.ref + '>' + al.alts.join(',') + ' here' };
      if (z === Z.HOM && al.alts.length === 1) return { state: 'copies', copies: 2 };
      if (z === Z.HET && al.alts.length === 1) return { state: 'copies', copies: 1 };
      if (z === Z.REF) return { state: 'copies', copies: 0 };
      return { state: z === Z.HET || z === Z.HOM ? 'other' : 'nocall', why: z === Z.HET || z === Z.HOM ? 'multi-allelic site' : 'genotype missing' };
    }
    if (!data.isGvcf) return { state: 'unknown', why: 'not in this VCF (a plain VCF lists variants only)' };
    var tr = data.tracks && data.tracks[c.key];
    if (!tr) return { state: 'nocall', why: 'nothing called on this chromosome' };
    var b = Math.floor((v.pos - 1) / tr.callable.binSize);
    if (tr.callable.levels[0][b] >= 0.999) return { state: 'reference', copies: 0 };
    if (tr.lowdp.levels[0][b] > 0) return { state: 'lowdp', why: 'low depth here' };
    return { state: 'nocall', why: 'not called here' };
  }

  G.pgx = { KB_VERSION: KB_VERSION, KB_REVIEWED: KB_REVIEWED, CPIC_LEVELS: CPIC_LEVELS, PGX_GENES: PGX_GENES, PGX_VARIANTS: PGX_VARIANTS,
    isPgxGene: isPgxGene, pgxGuidance: pgxGuidance, genotypeAt: genotypeAt };
})(globalThis.G = globalThis.G || {});
