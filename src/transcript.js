/*
 * A gene's canonical transcript (Ensembl), and what a variant does to its
 * protein. Pure functions over an Ensembl lookup and the coding sequence:
 *
 *   tx = transcript.fromEnsembl(lookupJson)      // canonical transcript, exons, coding positions
 *   transcript.exonsOnProtein(tx)                // residue range each coding exon encodes
 *   transcript.consequence(tx, cds, pos, ref, alt)
 *       -> { kind, residue, short, check }
 * kind: missense, synonymous, nonsense, stop lost, start lost, frameshift,
 * inframe, splice site (within 2 bases of an exon edge, in the intron),
 * UTR or intron. A plain translation on one transcript, for display; not a
 * full consequence predictor (no splice regions beyond 2 bases, no NMD).
 * check is 'ok' when the transcript's base matches the VCF's REF at the
 * site, 'mismatch' otherwise (then no amino acid change is named).
 * chrM genes use the vertebrate mitochondrial code (AGA and AGG stop, ATA
 * methionine, TGA tryptophan).
 */
(function (G) {
  var B = 'TCAG', AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  // vertebrate mitochondrial code: TGA (index 14) Trp, ATA (34) Met, AGA and AGG (46, 47) stop
  var AA_MITO = AA.slice(0, 14) + 'W' + AA.slice(15, 34) + 'M' + AA.slice(35, 46) + '**' + AA.slice(48);
  var THREE = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter' };
  var COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };

  function translate(codon, mito) {
    var i = B.indexOf(codon[0]) * 16 + B.indexOf(codon[1]) * 4 + B.indexOf(codon[2]);
    return i < 0 || codon.length !== 3 || /[^TCAG]/.test(codon) ? 'X' : (mito ? AA_MITO : AA)[i];
  }
  function revcomp(s) { return s.split('').reverse().map(function (c) { return COMP[c] || 'N'; }).join(''); }

  function fromEnsembl(j) {
    var t = (j.Transcript || []).find(function (x) { return x.is_canonical; }) || (j.Transcript || []).find(function (x) { return x.Translation; });
    if (!t || !t.Translation) return null;
    var strand = t.strand || j.strand, exons = (t.Exon || []).map(function (e) { return { start: e.start, end: e.end }; });
    exons.sort(function (a, b) { return strand > 0 ? a.start - b.start : b.start - a.start; }); // transcript order
    var cs = t.Translation.start, ce = t.Translation.end, cds = [], idx = new Map();
    exons.forEach(function (e) {
      if (strand > 0) { for (var p = Math.max(e.start, cs); p <= Math.min(e.end, ce); p++) { idx.set(p, cds.length); cds.push(p); } }
      else { for (var q = Math.min(e.end, ce); q >= Math.max(e.start, cs); q--) { idx.set(q, cds.length); cds.push(q); } }
    });
    return {
      gene: j.display_name, geneId: j.id, chrom: String(j.seq_region_name), strand: strand, id: t.id, name: t.display_name,
      start: t.start, end: t.end, exons: exons, cdsStart: cs, cdsEnd: ce, cds: cds, idx: idx,
      aaLength: t.Translation.length, mito: /^(MT|M|chrM)$/i.test(String(j.seq_region_name))
    };
  }

  // Residue range encoded by each exon (exon numbers in transcript order; non-coding exons left out).
  function exonsOnProtein(tx) {
    var out = [];
    tx.exons.forEach(function (e, i) {
      var first = null, last = null;
      for (var p = e.start; p <= e.end; p++) { var k = tx.idx.get(p); if (k === undefined) continue; if (first === null || k < first) first = k; if (last === null || k > last) last = k; }
      if (first !== null) out.push({ exon: i + 1, aaStart: Math.floor(first / 3) + 1, aaEnd: Math.floor(last / 3) + 1 });
    });
    return out;
  }

  function consequence(tx, cds, pos, ref, alt) {
    ref = String(ref).toUpperCase(); alt = String(alt).toUpperCase();
    var inCds = function (p) { return tx.idx.get(p); };
    var exonEdge = tx.exons.some(function (e) { return (pos >= e.start - 2 && pos <= e.start - 1) || (pos >= e.end + 1 && pos <= e.end + 2); });
    var inExon = tx.exons.some(function (e) { return pos >= e.start && pos <= e.end; });
    if (ref.length !== alt.length) { // an indel: anchored at pos (VCF), affects pos+1 onwards
      var hits = []; for (var p = pos; p < pos + Math.max(ref.length, 2); p++) if (inCds(p) !== undefined) hits.push(inCds(p));
      if (!hits.length) return exonEdge && !inExon ? { kind: 'splice site', check: 'ok' } : { kind: inExon ? 'UTR' : 'intron', check: 'ok' };
      var residue = Math.floor(Math.min.apply(null, hits) / 3) + 1, d = alt.length - ref.length;
      return { kind: d % 3 ? 'frameshift' : 'inframe', residue: residue, short: (d % 3 ? 'frameshift at ' : (d > 0 ? 'insertion at ' : 'deletion at ')) + residue, check: 'ok' };
    }
    if (ref.length !== 1) return { kind: inExon ? 'coding block change' : 'intron', check: 'ok' };
    var k = inCds(pos);
    if (k === undefined) return exonEdge && !inExon ? { kind: 'splice site', check: 'ok' } : { kind: inExon ? 'UTR' : 'intron', check: 'ok' };
    var rb = tx.strand > 0 ? ref : COMP[ref], ab = tx.strand > 0 ? alt : COMP[alt];
    if (!cds || cds[k] !== rb) return { kind: 'coding', residue: Math.floor(k / 3) + 1, check: 'mismatch' };
    var c0 = Math.floor(k / 3) * 3, codon = cds.substr(c0, 3), nc = codon.slice(0, k - c0) + ab + codon.slice(k - c0 + 1);
    var a0 = translate(codon, tx.mito), a1 = translate(nc, tx.mito), res = c0 / 3 + 1;
    var kind = a0 === a1 ? 'synonymous' : a1 === '*' ? 'nonsense' : a0 === '*' ? 'stop lost' : res === 1 && a0 === 'M' ? 'start lost' : 'missense';
    return { kind: kind, residue: res, ref: a0, alt: a1, short: a0 + res + a1, hgvs: 'p.' + (THREE[a0] || a0) + res + (a0 === a1 ? '=' : THREE[a1] || a1), codon: codon + '>' + nc, check: 'ok' };
  }

  G.transcript = { fromEnsembl: fromEnsembl, exonsOnProtein: exonsOnProtein, consequence: consequence, translate: translate, revcomp: revcomp };
})(globalThis.G = globalThis.G || {});
