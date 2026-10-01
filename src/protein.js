/*
 * Protein view: a lollipop plot for one gene. The protein is a bar with its
 * Pfam domains (InterPro API); above it, every ClinVar P/LP variant of the
 * gene placed by the protein change in its HGVS name (p.Gly551Asp -> 551),
 * stems grouped per residue with height by count; the sample's findings in
 * that gene are drawn in red on top. Shows where on the protein the known
 * pathogenic changes cluster, and where the sample's variant falls.
 *
 * Remote calls send only the gene symbol and its UniProt accession
 * (rest.uniprot.org, www.ebi.ac.uk/interpro). Results are cached per gene.
 */
(function (G) {
  var AA3 = { Ala: 'A', Arg: 'R', Asn: 'N', Asp: 'D', Cys: 'C', Gln: 'Q', Glu: 'E', Gly: 'G', His: 'H', Ile: 'I', Leu: 'L', Lys: 'K',
    Met: 'M', Phe: 'F', Pro: 'P', Ser: 'S', Thr: 'T', Trp: 'W', Tyr: 'Y', Val: 'V', Ter: '*' };
  var KIND_COLORS = { missense: 'rgb(120,180,255)', nonsense: 'rgb(255,160,60)', frameshift: 'rgb(235,90,200)', other: 'rgb(170,170,170)' };
  var DOMAIN_COLORS = ['#6f5bd3', '#4e79a7', '#59a14f', '#f28e2b', '#e15759', '#76b7b2', '#edc948'];

  // Protein position and kind from a ClinVar name such as
  // "NM_000492.4(CFTR):c.1652G>A (p.Gly551Asp)".
  function proteinChange(name) {
    var m = /\(p\.([A-Z][a-z]{2})(\d+)([A-Za-z*=]*)/.exec(name || '');
    if (!m) return null;
    var tail = m[3] || '', kind = /fs/.test(tail) ? 'frameshift' : tail === 'Ter' || tail === '*' ? 'nonsense' : AA3[tail] ? 'missense' : 'other';
    return { pos: +m[2], kind: kind, short: (AA3[m[1]] || m[1]) + m[2] + (AA3[tail] || tail) };
  }

  function ProteinView() { this.cache = {}; this.gene = null; }

  ProteinView.prototype.load = function (gene) {
    var self = this;
    if (this.cache[gene]) return this.cache[gene];
    var p = (async function () {
      var up = await fetch('https://rest.uniprot.org/uniprotkb/search?query=gene_exact:' + encodeURIComponent(gene) +
        '+AND+organism_id:9606+AND+reviewed:true&fields=accession,length,protein_name&format=json').then(function (r) { return r.json(); });
      var hit = up.results && up.results[0];
      if (!hit) throw new Error('No reviewed human UniProt entry for ' + gene);
      var acc = hit.primaryAccession, length = hit.sequence ? hit.sequence.length : 0;
      var domains = [];
      try {
        var ip = await fetch('https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/' + acc + '/?format=json').then(function (r) { return r.json(); });
        (ip.results || []).forEach(function (r) {
          r.proteins.forEach(function (pr) { pr.entry_protein_locations.forEach(function (loc) { loc.fragments.forEach(function (f) {
            domains.push({ start: f.start, end: f.end, name: r.metadata.name, acc: r.metadata.accession });
          }); }); });
        });
      } catch (e) { /* domains are optional */ }
      var name = hit.proteinDescription && hit.proteinDescription.recommendedName ? hit.proteinDescription.recommendedName.fullName.value : '';
      return { acc: acc, length: length, domains: domains, name: name };
    })();
    this.cache[gene] = p;
    p.then(function (v) { self.loaded = self.loaded || {}; self.loaded[gene] = v; }, function (e) { self.errors = self.errors || {}; self.errors[gene] = e.message; });
    return p;
  };

  // ClinVar P/LP variants of a gene, grouped by residue.
  ProteinView.prototype.clinvarFor = function (gene) {
    var cv = G.app.clinvar;
    if (!cv) return [];
    if (!this.byGene) {
      var map = new Map();
      cv.index.forEach(function (list) { list.forEach(function (e) {
        String(e.gene).split(/[;,]/).forEach(function (gname) { var l = map.get(gname); if (l) l.push(e); else map.set(gname, [e]); });
      }); });
      this.byGene = map;
    }
    var per = {};
    (this.byGene.get(gene) || []).forEach(function (e) {
      var pc = proteinChange(e.name);
      if (!pc) return;
      var k = pc.pos + ':' + pc.kind, slot = per[k] || (per[k] = { pos: pc.pos, kind: pc.kind, n: 0, names: [] });
      slot.n++; if (slot.names.length < 6) slot.names.push(pc.short);
    });
    return Object.keys(per).map(function (k) { return per[k]; });
  };

  ProteinView.prototype.draw = function (g) {
    var ctx = g.context, view = G.app.view;
    var msg = function (t) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(t, g.cX, g.cY); };
    var findingGenes = [];
    (view.findings || []).forEach(function (f) { String(f.gene).split(/[;,]/).forEach(function (x) { if (x && findingGenes.indexOf(x) < 0) findingGenes.push(x); }); });
    var gene = this.gene || view.focusGene || findingGenes[0];
    if (!gene) return msg('Pick a gene: a finding, search (Ctrl+K), or a click in the Arcs gene lane.');
    this.load(gene);
    if (this.errors && this.errors[gene]) return msg(this.errors[gene]);
    var p = this.loaded && this.loaded[gene];
    if (!p) return msg('Loading ' + gene + ' from UniProt and InterPro...');

    var left = 90, right = g.cW - 60, W = right - left, base = g.cY + 90, bx = function (pos) { return left + (pos - 0.5) / p.length * W; };
    var known = this.clinvarFor(gene), maxN = known.reduce(function (m, k) { return Math.max(m, k.n); }, 1);
    var over = null, bestD = 6;
    // protein bar and domains
    ctx.fillStyle = '#555'; ctx.fillRect(left, base - 6, W, 12);
    var names = [];
    p.domains.forEach(function (dm) {
      var i = names.indexOf(dm.name); if (i < 0) { names.push(dm.name); i = names.length - 1; }
      ctx.fillStyle = DOMAIN_COLORS[i % DOMAIN_COLORS.length];
      ctx.fillRect(bx(dm.start), base - 12, Math.max(2, bx(dm.end) - bx(dm.start)), 24);
      if (g.mY > base - 12 && g.mY < base + 12 && g.mX >= bx(dm.start) && g.mX <= bx(dm.end)) over = { lines: [dm.name + ' (' + dm.acc + ')', 'residues ' + dm.start + '-' + dm.end] };
    });
    g.setText('rgba(255,255,255,0.55)', 10, 'Helvetica, Arial, sans-serif', 'center', 'top');
    for (var t = 0; t <= p.length; t += Math.max(50, Math.round(p.length / 10 / 50) * 50)) { ctx.fillStyle = '#777'; ctx.fillRect(bx(t) , base + 12, 1, 4); g.fText(String(t), bx(t), base + 18); }
    // legend of domains
    g.setText('rgba(255,255,255,0.75)', 11, 'Helvetica, Arial, sans-serif', 'left', 'middle');
    names.forEach(function (n, i) { ctx.fillStyle = DOMAIN_COLORS[i % DOMAIN_COLORS.length]; ctx.fillRect(left + i * 0 , base + 44 + i * 16, 10, 10); g.fText(n, left + 16, base + 49 + i * 16); });

    // ClinVar P/LP lollipops
    known.forEach(function (k) {
      var x = bx(k.pos), h = 20 + 150 * Math.log(1 + k.n) / Math.log(1 + maxN);
      ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, base - 12); ctx.lineTo(x, base - 12 - h); ctx.stroke();
      ctx.fillStyle = KIND_COLORS[k.kind]; ctx.beginPath(); ctx.arc(x, base - 12 - h, 2 + Math.sqrt(k.n), 0, Math.PI * 2); ctx.fill();
      var dd = Math.hypot(g.mX - x, g.mY - (base - 12 - h));
      if (dd < bestD) { bestD = dd; over = { lines: ['residue ' + k.pos + ': ' + k.n + ' ClinVar P/LP ' + k.kind + ' variant' + (k.n > 1 ? 's' : ''), k.names.join(', ')] }; }
    });
    // the sample's findings in this gene
    (view.findings || []).forEach(function (f) {
      if (String(f.gene).split(/[;,]/).indexOf(gene) < 0) return;
      var pc = proteinChange(f.variant_name);
      if (!pc) return;
      var x = bx(pc.pos), top = base - 210;
      ctx.strokeStyle = 'rgb(255,70,70)'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(x, base - 12); ctx.lineTo(x, top); ctx.stroke();
      ctx.fillStyle = 'rgb(255,70,70)'; ctx.beginPath(); ctx.arc(x, top, 7, 0, Math.PI * 2); ctx.fill();
      g.setText('white', 12, 'Helvetica, Arial, sans-serif', 'center', 'bottom'); g.fText(pc.short + ' (this sample, ' + (f.zygosity || '') + ')', x, top - 10);
      if (Math.hypot(g.mX - x, g.mY - top) < 10) over = { lines: [f.variant_name, (f.classification || '') + ', ' + (f.zygosity || '') + ', GT ' + (f.gt || '?')] };
    });

    g.setText('white', 16, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText(gene + '  ' + (p.name || ''), left, 130);
    g.setText('rgba(255,255,255,0.55)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('UniProt ' + p.acc + ', ' + p.length + ' aa. Grey stems: ClinVar P/LP variants by residue (' + known.reduce(function (s, k) { return s + k.n; }, 0) +
      ' with a protein change), dot colour: missense blue, nonsense orange, frameshift pink. Red: this sample.', left, 152);
    if (findingGenes.length) {
      g.fText('Other finding genes: ' + findingGenes.filter(function (x) { return x !== gene; }).join(', ') + ' (pick with search, or click in the findings list)', left, 168);
    }
    if (over) view.drawTooltip(g, over.lines);
  };

  G.proteinChange = proteinChange;
  G.ProteinView = ProteinView; G.proteinChange = proteinChange;
})(globalThis.G = globalThis.G || {});
