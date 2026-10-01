/*
 * Polygenic scores from the PGS Catalog (pgscatalog.org), computed on the
 * loaded genome, with the rules of Asclepius app/prs.py (contracts
 * complement, match_allele, effect_dosage, prs_score, prs_percentile; their
 * test vectors are in tests/prs.test.js):
 *  - the score's effect allele is matched to REF or ALT; a strand flip is
 *    tried only for SNPs that are not A/T or C/G (palindromes are refused);
 *  - a missing genotype is unknown, never a dosage of 0 (Asclepius D4);
 *  - a percentile is given only against an ancestry-matched reference mean
 *    and sd that the user supplies; otherwise the score is raw and marked
 *    uncalibrated. Clinical judgement stays in Asclepius.
 *
 * Each score site falls in one of three groups:
 *  called    a VCF record at the site whose alleles match: dosage from the genotype;
 *  hom-ref   no record, and the site's bin is fully callable in a gVCF: the
 *            sample is 0/0, so the dosage is 2 when the effect allele is the
 *            reference base, else 0. The reference base comes from GRCh38 itself
 *            (Ensembl REST for scores up to 20,000 sites, fetched for every site
 *            so the request does not reveal which sites this genome lacks; or a
 *            file from tools/pgs_refbases.py for larger scores). Neither the
 *            score's other allele nor "effect = non-reference" can stand in: on
 *            HG002 the effect allele is the reference at 47% of called sites
 *            (PGS000058), and the other allele is the reference at 42 to 53%;
 *  unknown   anything else (not called, low depth, a plain VCF's absent site,
 *            no reference base, an unmatched or palindromic allele): left out and counted.
 * Only the score ID goes to the PGS Catalog and EBI; genotypes stay here.
 */
(function (G) {
  var COMP = { A: 'T', T: 'A', C: 'G', G: 'C' };
  var REST = 'https://www.pgscatalog.org/rest/';

  // ----- Asclepius app/prs.py, ported
  function complement(b) { return COMP[String(b || '').trim().toUpperCase()] || ''; }
  function isStrandAmbiguous(ref, alt) { var r = String(ref || '').toUpperCase(), a = String(alt || '').toUpperCase(); return !!COMP[r] && a === COMP[r]; }
  function matchAllele(effect, ref, alt) {
    var e = String(effect || '').toUpperCase(), r = String(ref || '').toUpperCase(), a = String(alt || '').toUpperCase();
    if (e === r) return 0;
    if (e === a) return 1;
    if (!isStrandAmbiguous(r, a)) {
      if (complement(e) === r) return 0;
      if (complement(e) === a) return 1;
    }
    return null;
  }
  // Copies of allele index idx in a genotype given as allele indices ([0,1], [1,1], ...); null when missing.
  function countAllele(gt, idx) {
    if (!gt || !gt.length || gt.some(function (x) { return x === null || x === undefined || x < 0; })) return null;
    return gt.filter(function (x) { return x === idx; }).length;
  }
  function effectDosage(gt, ref, alt, effect) {
    var idx = matchAllele(effect, ref, alt);
    return idx === null ? null : countAllele(gt, idx);
  }
  function homrefDosage(refBase, effect) { return String(effect || '').toUpperCase() === String(refBase || '').toUpperCase() ? 2 : 0; }
  function prsScore(pairs) { return pairs.reduce(function (s, p) { return s + (+p[0]) * (+p[1]); }, 0); }
  function erf(x) { // Abramowitz and Stegun 7.1.26
    var s = x < 0 ? -1 : 1; x = Math.abs(x);
    var t = 1 / (1 + 0.3275911 * x), y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return s * y;
  }
  function prsPercentile(score, mean, sd) {
    if (!(sd > 0)) return 50;
    var p = 100 * 0.5 * (1 + erf((score - mean) / sd / Math.SQRT2));
    return Math.max(0, Math.min(100, p));
  }

  // ----- the PGS Catalog
  async function getJson(url) { var r = await fetch(url); if (!r.ok) throw new Error('PGS Catalog: HTTP ' + r.status); return r.json(); }
  // Traits matching a term, each with its score IDs.
  async function searchTraits(term) {
    var d = await getJson(REST + 'trait/search?term=' + encodeURIComponent(term) + '&limit=20');
    return (d.results || []).map(function (t) { return { id: t.id, label: t.label, scores: t.associated_pgs_ids || [] }; }).filter(function (t) { return t.scores.length; });
  }
  async function scoreInfo(id) {
    var s = await getJson(REST + 'score/' + encodeURIComponent(id));
    var hm = (s.ftp_harmonized_scoring_files || {}).GRCh38 || {};
    return { id: s.id, name: s.name, trait: s.trait_reported, n: s.variants_number, url: hm.positions, ancestry: s.ancestry_distribution,
      publication: s.publication ? (s.publication.firstauthor + ' ' + (s.publication.date_publication || '').slice(0, 4) + ', ' + s.publication.journal) : '' };
  }

  // A harmonised scoring file: { meta, chr[], pos Int32Array, effect[], other[], weight Float64Array, n }.
  async function parseScoreFile(blob, opts) {
    var own = G.bgzf.own, meta = {}, cols = null, chr = [], pos = [], eff = [], oth = [], w = [], skipped = 0;
    for await (var line of G.bgzf.lines(blob, opts || {})) {
      if (line[0] === '#') { var eq = line.indexOf('='); if (eq > 0) meta[line.slice(1, eq)] = own(line.slice(eq + 1)); continue; }
      var f = line.split('\t');
      if (!cols) { cols = {}; f.forEach(function (c, i) { cols[c] = i; }); continue; }
      var c = cols.hm_chr !== undefined ? f[cols.hm_chr] : f[cols.chr_name], p = +(cols.hm_pos !== undefined ? f[cols.hm_pos] : f[cols.chr_position]);
      var wt = +f[cols.effect_weight];
      if (!c || !(p > 0) || !isFinite(wt)) { skipped++; continue; }
      var o = cols.other_allele !== undefined ? f[cols.other_allele] : '';
      if (!o && cols.hm_inferOtherAllele !== undefined) o = f[cols.hm_inferOtherAllele];
      chr.push(own(c)); pos.push(p); eff.push(own(f[cols.effect_allele] || '')); oth.push(own(o || '')); w.push(wt);
    }
    if (!cols) throw new Error('Not a PGS scoring file');
    return { meta: meta, chr: chr, pos: Int32Array.from(pos), effect: eff, other: oth, weight: Float64Array.from(w), n: pos.length, skipped: skipped };
  }

  // GRCh38 reference bases at every score site: Ensembl REST, 50 regions a request.
  // Returns an array aligned with the score's sites (null where not resolved).
  async function refBasesEnsembl(sf, onProgress) {
    var ref = new Array(sf.n).fill(null), batches = [];
    for (var i = 0; i < sf.n; i += 50) batches.push(i);
    var done = 0, next = 0;
    var worker = async function () {
      while (next < batches.length) {
        var b0 = batches[next++], regions = [];
        for (var i = b0; i < Math.min(sf.n, b0 + 50); i++) {
          var len = Math.max(1, sf.effect[i].length, (sf.other[i] || '').length);
          regions.push(String(sf.chr[i]).replace(/^chr/i, '') + ':' + sf.pos[i] + '..' + (sf.pos[i] + len - 1) + ':1');
        }
        for (var tries = 0; tries < 3; tries++) {
          var r = await fetch('https://rest.ensembl.org/sequence/region/human', { method: 'POST', headers: { 'Content-Type': 'application/json', Accept: 'application/json' }, body: JSON.stringify({ regions: regions }) });
          if (r.status === 429) { await new Promise(function (ok) { setTimeout(ok, 1500 * (tries + 1)); }); continue; }
          if (!r.ok) throw new Error('Ensembl sequence: HTTP ' + r.status);
          var byQuery = {};
          (await r.json()).forEach(function (x) { byQuery[x.query] = x.seq; });
          regions.forEach(function (q, k) { ref[b0 + k] = byQuery[q] ? byQuery[q].toUpperCase() : null; });
          break;
        }
        done++; if (onProgress) onProgress(done / batches.length);
      }
    };
    await Promise.all([worker(), worker(), worker(), worker()]);
    return ref;
  }

  // Reference bases from a tools/pgs_refbases.py file (chrom, pos, bases per line).
  async function refBasesFile(sf, blob) {
    var map = new Map(), ref = new Array(sf.n).fill(null);
    for await (var line of G.bgzf.lines(blob)) {
      if (!line || line[0] === '#') continue;
      var f = line.split('\t'); map.set(String(f[0]).replace(/^chr/i, '') + ':' + f[1], f[2]);
    }
    for (var i = 0; i < sf.n; i++) { var b = map.get(String(sf.chr[i]).replace(/^chr/i, '') + ':' + sf.pos[i]); if (b) ref[i] = b.toUpperCase(); }
    return ref;
  }

  function lowerBound(arr, n, v) { var lo = 0, hi = n; while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; } return lo; }

  // The genotype at record j of contig store c for sample s as allele indices, or null.
  function genotype(c, j, s) {
    var Z = G.vcf.Z, z = c.zygOf(j, s || 0);
    return z === Z.HOM ? [1, 1] : z === Z.HET ? [0, 1] : z === Z.REF ? [0, 0] : null;
  }

  // Score the loaded genome. Returns totals, groups, the other-allele check, per-site
  // contributions for drawing, and per-contig contribution tracks.
  function compute(data, sf, opts) {
    opts = opts || {};
    var s = opts.sample || 0, gvcf = !!data.isGvcf, tracks = {}, sites = [];
    var out = { called: 0, homref: 0, unknown: 0, unknownWhy: { 'no call': 0, 'no reference base': 0, 'allele mismatch': 0, 'palindrome or multi-allelic': 0, 'off this genome': 0 },
      score: 0, scoreCalled: 0, otherIsRef: 0, effectIsRef: 0, otherChecked: 0, n: sf.n };
    for (var i = 0; i < sf.n; i++) {
      var ct = data.genome.get(sf.chr[i]);
      if (!ct || sf.pos[i] > ct.length) { out.unknown++; out.unknownWhy['off this genome']++; continue; }
      var c = data.variants && data.variants[ct.key], dose = null, how = null, p = sf.pos[i];
      if (c && c.n) {
        for (var j = lowerBound(c.pos, c.n, p); j < c.n && c.pos[j] === p; j++) {
          var al = c.alleles(j);
          if (!al || al.alts.length !== 1) { how = how || 'palindrome or multi-allelic'; continue; }
          var idx = matchAllele(sf.effect[i], al.ref, al.alts[0]);
          if (how === null) { out.otherChecked++; if (String(sf.other[i] || '').toUpperCase() === al.ref.toUpperCase()) out.otherIsRef++; if (String(sf.effect[i]).toUpperCase() === al.ref.toUpperCase()) out.effectIsRef++; } // allele orientation, for the record
          if (idx === null) { how = isStrandAmbiguous(al.ref, al.alts[0]) ? 'palindrome or multi-allelic' : 'allele mismatch'; continue; }
          dose = countAllele(genotype(c, j, s), idx);
          how = dose === null ? 'no call' : 'called';
          break;
        }
      }
      if (how === null && gvcf) { // no record: hom-ref when the bin is fully callable
        var tr = data.tracks[ct.key], bin = tr ? Math.floor((p - 1) / tr.callable.binSize) : -1;
        var rb = sf.ref && sf.ref[i];
        if (tr && tr.callable.levels[0][bin] >= 0.999) {
          if (rb) { dose = homrefDosage(rb.slice(0, sf.effect[i].length), sf.effect[i]); how = 'homref'; }
          else how = 'no reference base';
        }
      }
      if (how === 'called' || how === 'homref') {
        var contrib = sf.weight[i] * dose;
        out.score += contrib; if (how === 'called') { out.called++; out.scoreCalled += contrib; } else out.homref++;
        if (contrib) {
          var t = tracks[ct.key] || (tracks[ct.key] = { pos: new G.genome.Track(ct.length, data.binSize), neg: new G.genome.Track(ct.length, data.binSize) });
          (contrib > 0 ? t.pos : t.neg).add(p, Math.abs(contrib));
          sites.push({ key: ct.key, chrom: ct.name, pos: p, effect: sf.effect[i], weight: sf.weight[i], dose: dose, contrib: contrib, how: how });
        }
      } else { out.unknown++; out.unknownWhy[how || 'no call']++; }
    }
    Object.keys(tracks).forEach(function (k) { tracks[k].pos.buildPyramid(); tracks[k].neg.buildPyramid(); });
    sites.sort(function (a, b) { return Math.abs(b.contrib) - Math.abs(a.contrib); });
    out.top = sites.slice(0, 200); out.tracks = tracks; out.coverage = (out.called + out.homref) / Math.max(1, sf.n);
    out.otherRefRate = out.otherChecked ? out.otherIsRef / out.otherChecked : null;
    out.effectRefRate = out.otherChecked ? out.effectIsRef / out.otherChecked : null;
    if (opts.refMean !== undefined && opts.refSd > 0) out.percentile = prsPercentile(out.score, opts.refMean, opts.refSd);
    return out;
  }

  G.prs = { complement: complement, isStrandAmbiguous: isStrandAmbiguous, matchAllele: matchAllele, countAllele: countAllele, effectDosage: effectDosage,
    homrefDosage: homrefDosage, prsScore: prsScore, prsPercentile: prsPercentile, searchTraits: searchTraits, scoreInfo: scoreInfo, parseScoreFile: parseScoreFile, refBasesEnsembl: refBasesEnsembl, refBasesFile: refBasesFile, compute: compute };
})(globalThis.G = globalThis.G || {});
