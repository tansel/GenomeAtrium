/*
 * Protein view: a lollipop plot for one gene. The protein is a bar with its
 * Pfam domains (InterPro API); above it, every ClinVar P/LP variant of the
 * gene placed by the protein change in its HGVS name (p.Gly551Asp -> 551),
 * stems grouped per residue with height by count; the sample's findings in
 * that gene are drawn in red on top. Shows where on the protein the known
 * pathogenic changes cluster, and where the sample's variant falls.
 *
 * Around the protein:
 *  - the sequence, as letters once zoomed in (wheel to zoom, drag to pan,
 *    double click to reset);
 *  - every coding variant of the sample in the gene, placed by its own
 *    consequence on the canonical transcript (transcript.js), not only the
 *    findings: missense, nonsense, frameshift, in frame, synonymous;
 *  - strips: AlphaFold confidence (pLDDT), AlphaMissense (mean predicted
 *    pathogenicity per residue), and the exon each residue comes from;
 *  - the gene model in genome coordinates: exons (coding thick, UTR thin),
 *    introns, the sample's variants along it, callability (gVCF) and
 *    methylation (when a methylation file is loaded), with the coding exons
 *    joined to the residues they encode;
 *  - buttons to see the protein in 3D, normal and variant side by side, or
 *    in the Protein room.
 *
 * Remote calls send only the gene symbol and its UniProt accession
 * (rest.uniprot.org, www.ebi.ac.uk/interpro, rest.ensembl.org for the
 * canonical transcript and its coding sequence, alphafold.ebi.ac.uk for the
 * confidence and AlphaMissense files). Results are cached per gene.
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

  var CONS_COL = { missense: 'rgb(255,90,90)', nonsense: 'rgb(255,160,60)', 'stop lost': 'rgb(255,160,60)', 'start lost': 'rgb(255,160,60)',
    frameshift: 'rgb(235,90,200)', inframe: 'rgb(200,120,255)', synonymous: 'rgb(150,150,150)', 'splice site': 'rgb(255,220,80)' };
  var AA_COL = { hydrophobic: 'rgb(230,200,90)', polar: 'rgb(120,200,160)', positive: 'rgb(120,160,255)', negative: 'rgb(255,120,120)', special: 'rgb(200,200,200)' };
  var AA_CLASS = { A: 'hydrophobic', V: 'hydrophobic', L: 'hydrophobic', I: 'hydrophobic', M: 'hydrophobic', F: 'hydrophobic', W: 'hydrophobic', C: 'hydrophobic',
    S: 'polar', T: 'polar', N: 'polar', Q: 'polar', Y: 'polar', K: 'positive', R: 'positive', H: 'positive', D: 'negative', E: 'negative', G: 'special', P: 'special' };

  function ProteinView() { this.cache = {}; this.gene = null; this.extra = {}; this.zoom = null; }

  ProteinView.prototype.install = function (canvas) {
    var self = this;
    canvas.addEventListener('wheel', function (e) {
      if (G.app.view.activeMode() !== 'protein' || !self.axis) return;
      var ax = self.axis, f = Math.exp(e.deltaY * 0.0018), r = ax.a + (e.offsetX - ax.left) / ax.W * (ax.b - ax.a);
      var a = r - (r - ax.a) * f, b = r + (ax.b - r) * f;
      self.setZoom(a, b);
    }, { passive: true });
    canvas.addEventListener('dblclick', function () { if (G.app.view.activeMode() === 'protein') self.zoom = null; });
  };
  ProteinView.prototype.setZoom = function (a, b) {
    var L = this.axis ? this.axis.L : 1, span = Math.max(20, Math.min(L, b - a));
    a = Math.max(0.5, Math.min(L + 0.5 - span, a)); this.zoom = span >= L ? null : { a: a, b: a + span };
  };

  ProteinView.prototype.load = function (gene) {
    var self = this;
    if (this.cache[gene]) return this.cache[gene];
    var p = (async function () {
      var up = await G.net.fetchRetry('https://rest.uniprot.org/uniprotkb/search?query=gene_exact:' + encodeURIComponent(gene) +
        '+AND+organism_id:9606+AND+reviewed:true&fields=accession,length,protein_name,sequence&format=json').then(function (r) { return r.json(); });
      var hit = up.results && up.results[0];
      if (!hit) throw new Error('No reviewed human UniProt entry for ' + gene);
      var acc = hit.primaryAccession, length = hit.sequence ? hit.sequence.length : 0;
      var domains = [];
      try {
        var ip = await G.net.fetchRetry('https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/' + acc + '/?format=json').then(function (r) { return r.json(); });
        (ip.results || []).forEach(function (r) {
          r.proteins.forEach(function (pr) { pr.entry_protein_locations.forEach(function (loc) { loc.fragments.forEach(function (f) {
            domains.push({ start: f.start, end: f.end, name: r.metadata.name, acc: r.metadata.accession });
          }); }); });
        });
      } catch (e) { /* domains are optional */ }
      var name = hit.proteinDescription && hit.proteinDescription.recommendedName ? hit.proteinDescription.recommendedName.fullName.value : '';
      return { acc: acc, length: length, domains: domains, name: name, sequence: hit.sequence && hit.sequence.value || '' };
    })();
    this.cache[gene] = p;
    p.then(function (v) { self.loaded = self.loaded || {}; self.loaded[gene] = v; self.loadExtra(gene, v); }, function (e) { // shown for 15 s, then the next draw tries again (a failure is not cached)
      self.errors = self.errors || {}; self.errors[gene] = e.message + ' (trying again shortly)';
      setTimeout(function () { delete self.cache[gene]; delete self.errors[gene]; }, 15000);
    });
    return p;
  };

  // In the background: the canonical transcript and its coding sequence (Ensembl, the genome's
  // build), AlphaFold confidence and AlphaMissense. Each part is optional.
  ProteinView.prototype.loadExtra = function (gene, p) {
    var x = this.extra[gene] = this.extra[gene] || {}, build = G.app.view.data && G.app.view.data.build;
    var host = build === 'GRCh37' ? 'https://grch37.rest.ensembl.org' : 'https://rest.ensembl.org';
    if (!x.txLoading && !x.cds && (build === 'GRCh38' || build === 'GRCh37')) {
      x.txLoading = true; x.build = build;
      x.txPromise = G.net.fetchRetry(host + '/lookup/symbol/homo_sapiens/' + encodeURIComponent(gene) + '?expand=1;content-type=application/json', { signal: AbortSignal.timeout(45000) }).then(function (r) { return r.ok ? r.json() : null; })
        .then(function (j) {
          var tx = j && G.transcript.fromEnsembl(j);
          if (!tx) { x.txError = 'no protein-coding transcript in Ensembl'; return; }
          x.tx = tx; x.exonsAA = G.transcript.exonsOnProtein(tx);
          return G.net.fetchRetry(host + '/sequence/id/' + tx.id + '?type=cds;content-type=text/plain', { signal: AbortSignal.timeout(45000) }).then(function (r) { return r.ok ? r.text() : ''; }).then(function (t) { x.cds = t.trim().toUpperCase(); if (!x.cds) x.cdsError = 'Ensembl gave no coding sequence'; });
        }).catch(function (e) { // Ensembl busy or offline: say so, and try again in 15 s
          if (x.tx) x.cdsError = 'Ensembl coding sequence: ' + e.message + ' (retrying)'; else x.txError = e.message + ' (retrying)';
          setTimeout(function () { x.txLoading = false; x.txError = x.cdsError = null; }, 15000);
        });
    }
    if (!x.afLoading) {
      x.afLoading = true;
      G.net.fetchRetry('https://alphafold.ebi.ac.uk/api/prediction/' + p.acc).then(function (r) { return r.ok ? r.json() : []; }).then(function (list) {
        var e = (list || []).find(function (q) { return q.uniprotAccession === p.acc; });
        if (!e) return;
        if (e.plddtDocUrl) G.net.fetchRetry(e.plddtDocUrl).then(function (r) { return r.json(); }).then(function (c) { x.plddt = c.confidenceScore; }).catch(function () {});
        if (e.amAnnotationsUrl) G.structure.fetchAlphaMissense(e.amAnnotationsUrl).then(function (am) { x.am = am; }).catch(function () {});
      }).catch(function () {});
    }
  };

  // Resolves with the gene's extras once the transcript and coding sequence are read (or failed).
  ProteinView.prototype.ready = function (gene) {
    var self = this;
    return this.load(gene).then(function (p) {
      self.loadExtra(gene, p);
      return self.extra[gene].txPromise;
    }).then(function () { return self.extra[gene]; }, function () { return self.extra[gene] || {}; });
  };

  // The gene the Protein view (and the Atrium's 3D protein) shows, the same rule everywhere:
  //  1. a gene picked here (pv.gene);
  //  2. Arcs zoomed in under 3 Mb: the focused gene if it is in view, else, when the view
  //     moved after the gene was focused, the protein-coding gene nearest the centre;
  //  3. the focused gene (search, a click in a gene lane, Mito, a finding);
  //  4. the first finding.
  ProteinView.prototype.currentGene = function () {
    var view = G.app.view, d = view.data;
    if (this.gene) return this.gene;
    var fg = view.focusGene;
    if (view.visibleSpan && view.zoom > 1.01 && d && d.build === 'GRCh38' && view.genes) { // the gene table is GRCh38
      var sp = view.visibleSpan();
      if (sp && sp.end - sp.start < 3e6) {
        var key = G.genome.normName(sp.chrom), list = view.genes.inRange(key, sp.start, sp.end).filter(function (q) { return q.type === 'protein_coding'; });
        if (fg && list.some(function (q) { return q.name === fg; })) return fg;
        if (list.length && (!fg || (view.lastMove || 0) > (view.focusAt || 0))) {
          var mid = (sp.start + sp.end) / 2;
          list.sort(function (a, b) { return Math.max(0, a.start - mid, mid - a.end) - Math.max(0, b.start - mid, mid - b.end) || (b.end - b.start) - (a.end - a.start); });
          return list[0].name;
        }
      }
    }
    if (fg) return fg;
    var f = (view.findings || []).map(function (q) { return String(q.gene).split(/[;,]/)[0]; }).filter(function (q) { return /^[A-Za-z0-9]/.test(q); })[0];
    return f || null;
  };

  // The sample's variants in the gene (canonical transcript span, plus 2 bases), with consequences.
  ProteinView.prototype.sampleVariants = function (gene) {
    var x = this.extra[gene], d = G.app.view.data;
    if (!x || !x.tx || !x.cds || !d || !d.variants || d.build !== x.build) return []; // wait for the coding sequence
    if (x.sv && x.svFor === d && x.svCds === !!x.cds) return x.sv;
    var c = d.genome.get(x.tx.chrom), col = c && d.variants[c.key], Z = G.vcf.Z, out = [];
    if (col) {
      var lo = 0, hi = col.n, a = x.tx.start - 2;
      while (lo < hi) { var mid = (lo + hi) >> 1; if (col.pos[mid] < a) lo = mid + 1; else hi = mid; }
      for (var i = lo; i < col.n && col.pos[i] <= x.tx.end + 2; i++) {
        var z = col.zygOf(i, 0);
        if (z !== Z.HET && z !== Z.HOM) continue;
        var al = col.alleles(i);
        if (!al) continue;
        var cq = G.transcript.consequence(x.tx, x.cds, col.pos[i], al.ref, al.alts[0]);
        out.push({ pos: col.pos[i], ref: al.ref, alt: al.alts[0], zyg: z === Z.HOM ? 'homozygous' : 'heterozygous', cq: cq, chrom: c.name });
      }
    }
    x.sv = out; x.svFor = d; x.svCds = !!x.cds;
    return out;
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
    var ctx = g.context, view = G.app.view, self = this;
    var msg = function (t) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(t, g.cX, g.cY); };
    var findingGenes = [];
    (view.findings || []).forEach(function (f) { String(f.gene).split(/[;,]/).forEach(function (x) { if (/^[A-Za-z0-9]/.test(x || '') && findingGenes.indexOf(x) < 0) findingGenes.push(x); }); });
    var gene = this.currentGene();
    if (!gene) return msg('Pick a gene: a finding, search (Ctrl+K), or a click in the Arcs gene lane.');
    if (gene !== this.zoomGene) { this.zoom = null; this.zoomGene = gene; }
    this.load(gene);
    if (this.errors && this.errors[gene]) return msg(this.errors[gene]);
    var p = this.loaded && this.loaded[gene];
    if (!p) return msg('Loading ' + gene + ' from UniProt and InterPro...');
    this.loadExtra(gene, p); // no-op once loaded or loading; retries after a failed Ensembl call
    var x = this.extra[gene] || {}, L = p.length || (p.sequence || '').length || 1;

    // the residue axis (zoomable)
    var left = 130, right = g.cW - 60, W = right - left, z = this.zoom || { a: 0.5, b: L + 0.5 };
    this.axis = { left: left, W: W, a: z.a, b: z.b, L: L };
    var rx = function (r) { return left + (r - 0.5 - (z.a - 0.5)) / (z.b - z.a) * W; }, ppr = W / (z.b - z.a);
    // vertical layout from the height available: the page (HTML bars over the top 110 px) or the
    // Atrium's window (1000 x 500, nothing over it). lolH is the lollipop zone above the bar.
    var inPanel = view.mode === 'atrium', T0 = inPanel ? 12 : 118;
    var lolH = Math.max(110, Math.min(210, g.cH - (T0 + 62) - 250)), top = T0 + 62, base = top + lolH, over = null, bestD = 7;
    if (g.MOUSE_PRESSED && g.DX_MOUSE && g.mY > top - 20 && g.mY < base + 80) this.setZoom(z.a - g.DX_MOUSE / ppr, z.b - g.DX_MOUSE / ppr);
    ctx.save(); ctx.beginPath(); ctx.rect(left - 4, 0, W + 8, g.cH); ctx.clip();

    // ClinVar P/LP lollipops (by residue)
    var known = this.clinvarFor(gene), maxN = known.reduce(function (m, k) { return Math.max(m, k.n); }, 1);
    known.forEach(function (k) {
      var xx = rx(k.pos), h = 10 + lolH * 0.42 * Math.log(1 + k.n) / Math.log(1 + maxN);
      ctx.strokeStyle = 'rgba(255,255,255,0.25)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(xx, base - 12); ctx.lineTo(xx, base - 12 - h); ctx.stroke();
      ctx.fillStyle = KIND_COLORS[k.kind]; ctx.beginPath(); ctx.arc(xx, base - 12 - h, 1.5 + Math.sqrt(k.n), 0, Math.PI * 2); ctx.fill();
      var dd = Math.hypot(g.mX - xx, g.mY - (base - 12 - h));
      if (dd < bestD) { bestD = dd; over = ['residue ' + k.pos + ': ' + k.n + ' ClinVar P/LP ' + k.kind + ' variant' + (k.n > 1 ? 's' : ''), k.names.join(', ')]; }
    });
    // the sample's coding variants, by consequence (diamonds above the ClinVar stems)
    if (x.tx && !x.cds) { g.setText('rgba(255,255,255,0.45)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top'); g.fText(x.cdsError || 'Reading the coding sequence from Ensembl to place this sample\'s variants...', left, T0 + 54); }
    var sv = this.sampleVariants(gene), coding = sv.filter(function (v) { return v.cq.residue; });
    coding.forEach(function (v, i) {
      var xx = rx(v.cq.residue), yy = base - lolH * 0.64 - (i % 3) * 11, col = CONS_COL[v.cq.kind] || 'rgb(200,200,200)';
      ctx.strokeStyle = col; ctx.lineWidth = v.cq.kind === 'synonymous' ? 1 : 2;
      ctx.beginPath(); ctx.moveTo(xx, base - 12); ctx.lineTo(xx, yy); ctx.stroke();
      ctx.fillStyle = col; ctx.beginPath(); ctx.moveTo(xx, yy - 6); ctx.lineTo(xx + 5, yy); ctx.lineTo(xx, yy + 6); ctx.lineTo(xx - 5, yy); ctx.closePath(); ctx.fill();
      if (v.cq.kind !== 'synonymous' && (coding.length < 25 || ppr > 3)) { g.setText('white', 11, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText(v.cq.short || v.cq.kind, xx - 8, yy); }
      var dd = Math.hypot(g.mX - xx, g.mY - yy);
      if (dd < bestD) {
        bestD = dd; over = ['this sample: ' + (v.cq.hgvs || v.cq.short || v.cq.kind) + ', ' + v.cq.kind + ' (' + v.zyg + ')', v.chrom + ':' + v.pos.toLocaleString() + ' ' + v.ref + '>' + v.alt + (v.cq.codon ? ', codon ' + v.cq.codon : ''),
          v.cq.check === 'mismatch' ? 'the transcript\'s base differs from the VCF REF here: no change named' : 'click: open in Arcs'];
        g.setCursor('pointer');
        if (g.MOUSE_UP_FAST) { G.app.setMode('arcs'); view.goTo(v.chrom, v.pos - 60, v.pos + 60); }
      }
    });
    // the sample's findings (from the findings file): tall red stems
    (view.findings || []).forEach(function (f) {
      if (String(f.gene).split(/[;,]/).indexOf(gene) < 0) return;
      var pc = proteinChange(f.variant_name);
      if (!pc) return;
      var xx = rx(pc.pos), yy = base - lolH * 0.9;
      ctx.strokeStyle = 'rgb(255,70,70)'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(xx, base - 12); ctx.lineTo(xx, yy); ctx.stroke();
      ctx.fillStyle = 'rgb(255,70,70)'; ctx.beginPath(); ctx.arc(xx, yy, 7, 0, Math.PI * 2); ctx.fill();
      g.setText('white', 12, 'Helvetica, Arial, sans-serif', 'center', 'bottom'); g.fText(pc.short + ' (finding, ' + (f.zygosity || '') + ')', xx, yy - 10);
      if (Math.hypot(g.mX - xx, g.mY - yy) < 10) over = [f.variant_name, (f.classification || '') + ', ' + (f.zygosity || '') + ', GT ' + (f.gt || '?')];
    });
    // the protein bar, domains and the sequence
    ctx.fillStyle = '#4a4d55'; ctx.fillRect(rx(0.5), base - 6, rx(L + 0.5) - rx(0.5), 12);
    var names = [], labelEnd = -1e9;
    p.domains.slice().sort(function (a, b) { return a.start - b.start; }).forEach(function (dm) {
      var i = names.indexOf(dm.name); if (i < 0) { names.push(dm.name); i = names.length - 1; }
      ctx.fillStyle = DOMAIN_COLORS[i % DOMAIN_COLORS.length];
      ctx.fillRect(rx(dm.start - 0.5), base - 11, Math.max(2, rx(dm.end + 0.5) - rx(dm.start - 0.5)), 22);
      var lx = Math.max(left, rx(dm.start - 0.5)) + 4;
      g.setText('white', 10, 'Helvetica, Arial, sans-serif', 'left', 'middle');
      if (lx > labelEnd && rx(dm.end) - lx > g.getTextW(dm.name)) { g.fText(dm.name, lx, base); labelEnd = lx + g.getTextW(dm.name) + 6; }
      if (g.mY > base - 11 && g.mY < base + 11 && g.mX >= rx(dm.start - 0.5) && g.mX <= rx(dm.end + 0.5)) over = [dm.name + ' (' + dm.acc + ')', 'residues ' + dm.start + '-' + dm.end];
    });
    var seq = p.sequence || '';
    if (seq && ppr >= 7) { // letters, coloured by class
      for (var r = Math.max(1, Math.floor(z.a)); r <= Math.min(L, Math.ceil(z.b)); r++) {
        var aa = seq[r - 1]; g.setText(AA_COL[AA_CLASS[aa]] || '#ccc', Math.min(13, ppr * 0.8), 'Menlo, monospace', 'center', 'middle'); g.fText(aa, rx(r), base + 20);
      }
    } else {
      g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'center', 'top');
      var step = [1, 2, 5, 10, 20, 50, 100, 200, 500].find(function (s2) { return s2 * ppr > 60; }) || 1000;
      for (var t = Math.ceil(z.a / step) * step; t <= z.b; t += step) { ctx.fillStyle = '#777'; ctx.fillRect(rx(t), base + 12, 1, 4); g.fText(String(t), rx(t), base + 16); }
    }
    if (g.mY > base - 11 && g.mY < base + 28 && !over && seq) { var rr = Math.round(z.a - 0.5 + (g.mX - left) / ppr + 0.5); if (rr >= 1 && rr <= L) over = ['residue ' + rr + ': ' + seq[rr - 1] + ' (' + (AA_CLASS[seq[rr - 1]] || '') + ')']; }
    // strips: AlphaFold confidence, AlphaMissense, exons
    var y = base + 34;
    if (x.cds && ppr >= 26) { // the codon of each residue (coding strand), past about 26 px a residue
      g.setText('rgba(255,255,255,0.55)', Math.min(11, ppr / 3.2), 'Menlo, monospace', 'center', 'middle');
      for (var rc = Math.max(1, Math.floor(z.a)); rc <= Math.min(L, Math.ceil(z.b)); rc++) g.fText(x.cds.substr((rc - 1) * 3, 3), rx(rc), y + 4);
      ctx.restore(); g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText('codons', left - 6, y + 4);
      ctx.save(); ctx.beginPath(); ctx.rect(left - 4, 0, W + 8, g.cH); ctx.clip();
      y += 16;
    }
    var strip = function (label, get, col) {
      ctx.restore(); g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText(label, left - 6, y + 4); ctx.save(); ctx.beginPath(); ctx.rect(left - 4, 0, W + 8, g.cH); ctx.clip();
      var r0 = Math.max(1, Math.floor(z.a)), r1 = Math.min(L, Math.ceil(z.b)), stepR = Math.max(1, Math.floor((r1 - r0) / W));
      for (var r2 = r0; r2 <= r1; r2 += stepR) { var v2 = get(r2); if (v2 === null || v2 === undefined) continue; ctx.fillStyle = col(v2); ctx.fillRect(rx(r2 - 0.5), y, Math.max(1, ppr * stepR), 8); }
      if (g.mY >= y && g.mY < y + 8) { var rr2 = Math.round(z.a - 0.5 + (g.mX - left) / ppr + 0.5); var vv = get(rr2); if (vv !== null && vv !== undefined) over = [label + ' at residue ' + rr2 + ': ' + (typeof vv === 'number' ? vv.toFixed(2) : vv)]; }
      y += 12;
    };
    if (x.plddt) strip('AlphaFold confidence', function (r3) { return x.plddt[r3 - 1]; }, function (v3) { return G.structure.plddtColor(v3); });
    if (x.am) strip('AlphaMissense', function (r3) { return x.am.mean[r3]; }, function (v3) { var c3 = Math.round(255 * v3); return 'rgb(' + (60 + c3 * 0.75 | 0) + ',' + (110 - Math.abs(v3 - 0.5) * 80 | 0) + ',' + (230 - c3 * 0.75 | 0) + ')'; });
    var exAA = x.exonsAA || [];
    if (exAA.length) strip('exons', function (r3) { var e3 = exAA.find(function (q) { return r3 >= q.aaStart && r3 <= q.aaEnd; }); return e3 ? 'exon ' + e3.exon : null; }, function (v3) { return +v3.split(' ')[1] % 2 ? 'rgb(120,150,200)' : 'rgb(80,100,140)'; });
    var exonY = y - 12;
    ctx.restore();

    // the gene model (genome coordinates), joined to the exons above
    var tx = x.tx, gy = y + 36, endY = gy + 24;
    if (tx && x.build === (view.data && view.data.build)) {
      var g0 = tx.start - 2000, g1 = tx.end + 1000, gx = function (q) { return left + (q - g0) / (g1 - g0) * W; };
      if (tx.strand < 0) { g0 = tx.start - 1000; g1 = tx.end + 2000; }
      ctx.strokeStyle = 'rgba(255,255,255,0.45)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(gx(tx.start), gy + 8); ctx.lineTo(gx(tx.end), gy + 8); ctx.stroke();
      for (var ar = tx.start; ar < tx.end; ar += (tx.end - tx.start) / 12) { var ax2 = gx(ar); ctx.beginPath(); ctx.moveTo(ax2 - 3 * tx.strand, gy + 5); ctx.lineTo(ax2, gy + 8); ctx.lineTo(ax2 - 3 * tx.strand, gy + 11); ctx.stroke(); }
      tx.exons.forEach(function (e, i) {
        var cs = Math.max(e.start, tx.cdsStart), ce = Math.min(e.end, tx.cdsEnd);
        ctx.fillStyle = 'rgba(160,180,220,0.6)'; ctx.fillRect(gx(e.start), gy + 4, Math.max(1, gx(e.end + 1) - gx(e.start)), 8); // UTR height
        if (ce >= cs) { ctx.fillStyle = 'rgb(120,150,200)'; ctx.fillRect(gx(cs), gy, Math.max(1.5, gx(ce + 1) - gx(cs)), 16); }
        var ea = exAA.find(function (q) { return q.exon === i + 1; });
        if (ea && ce >= cs) { // join to the residues it encodes
          var ra = rx(ea.aaStart - 0.5), rb = rx(ea.aaEnd + 0.5), ga = gx(tx.strand > 0 ? cs : ce + 1), gb = gx(tx.strand > 0 ? ce + 1 : cs);
          if (rb > left - 10 && ra < right + 10) {
            ctx.fillStyle = 'rgba(120,150,200,0.12)'; ctx.beginPath(); ctx.moveTo(Math.max(left, ra), exonY + 8); ctx.lineTo(Math.min(right, rb), exonY + 8); ctx.lineTo(gb, gy); ctx.lineTo(ga, gy); ctx.closePath(); ctx.fill();
          }
        }
        if (g.mY >= gy && g.mY <= gy + 16 && g.mX >= gx(e.start) - 1 && g.mX <= gx(e.end + 1) + 1) over = ['exon ' + (i + 1) + ' of ' + tx.exons.length + ' (' + tx.name + ')', tx.chrom + ':' + e.start.toLocaleString() + '-' + e.end.toLocaleString() + ' (' + (e.end - e.start + 1) + ' bp)' + (ea ? ', residues ' + ea.aaStart + '-' + ea.aaEnd : ', non-coding')];
      });
      g.setText('rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle');
      g.fText(tx.name + ' (' + tx.id + ', canonical), chr' + String(tx.chrom).replace(/^chr/i, '') + (tx.strand > 0 ? ' + strand' : ' - strand'), right, gy - 10);
      // the sample's variants along the gene
      var vy = gy + 24;
      sv.forEach(function (v) {
        var xx = gx(v.pos), col = CONS_COL[v.cq.kind] || 'rgba(200,200,200,0.6)';
        ctx.fillStyle = col; ctx.fillRect(xx - 0.5, vy, 1.5, v.zyg === 'homozygous' ? 12 : 7);
        if (Math.abs(g.mX - xx) < 3 && g.mY >= vy && g.mY < vy + 12) over = ['this sample: ' + v.chrom + ':' + v.pos.toLocaleString() + ' ' + v.ref + '>' + v.alt + ' (' + v.zyg + ')', v.cq.kind + (v.cq.hgvs ? ', ' + v.cq.hgvs : '')];
      });
      g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText('variants', left - 6, vy + 5);
      var sy = vy + 18, d = view.data, ckey = d.genome.get(tx.chrom) && d.genome.get(tx.chrom).key, tr = ckey && d.tracks && d.tracks[ckey];
      if (d.isGvcf && tr) { // calls along the gene
        var bs = tr.callable.binSize;
        for (var px = left; px < right; px += 2) {
          var q0 = g0 + (px - left) / W * (g1 - g0), b0 = Math.floor((q0 - 1) / bs), cv = tr.callable.levels[0][b0] || 0, lw = tr.lowdp.levels[0][b0] || 0;
          ctx.fillStyle = lw > 0.05 ? 'rgb(230,170,40)' : cv > 0.5 ? 'rgba(255,255,255,0.35)' : 'rgb(220,70,70)'; ctx.fillRect(px, sy, 2, 4);
        }
        g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText('calls', left - 6, sy + 2); sy += 10;
      }
      var me = G.app.methyl;
      if (me && d.methyl === me && ckey) { // methylation along the gene and its promoter
        var mh = inPanel ? 20 : 30;
        ctx.strokeStyle = 'rgba(255,255,255,0.15)'; ctx.strokeRect(left, sy, W, mh);
        for (var px2 = left; px2 < right; px2 += 3) {
          var q1 = g0 + (px2 - left) / W * (g1 - g0), lv = me.level(ckey, Math.round(q1), Math.round(q1 + (g1 - g0) / W * 3));
          if (lv.frac === null) continue;
          var yy2 = sy + mh - mh * lv.frac; ctx.fillStyle = G.methylation.Methylation.color(lv.frac); ctx.fillRect(px2, yy2 - 1, 3, 2);
        }
        g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle'); g.fText('methylation', left - 6, sy + mh / 2);
        if (g.mY >= sy && g.mY < sy + mh && g.mX >= left && g.mX < right) { var q2 = Math.round(g0 + (g.mX - left) / W * (g1 - g0)), lv2 = me.level(ckey, q2 - 500, q2 + 500); over = ['methylation around ' + tx.chrom + ':' + q2.toLocaleString() + ': ' + (lv2.frac === null ? 'too few calls' : Math.round(100 * lv2.frac) + '%')]; }
        sy += mh + 6;
      }
      endY = sy;
    } else {
      g.setText('rgba(255,255,255,0.45)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
      g.fText(x.txError ? 'Gene model: ' + x.txError : view.data && (view.data.build === 'GRCh38' || view.data.build === 'GRCh37') ? 'Loading the gene model from Ensembl...' : 'The gene model needs a GRCh37 or GRCh38 genome.', left, gy);
    }

    // title, buttons, key
    g.setText('white', 16, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText(gene + '  ' + (p.name || ''), left, T0);
    g.setText('rgba(255,255,255,0.55)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('UniProt ' + p.acc + ', ' + L + ' aa. Stems: ClinVar P/LP by residue (dot colour: missense blue, nonsense orange, frameshift pink). Diamonds: this sample\'s coding variants (' + coding.length + '), red: findings.', left, T0 + 22);
    g.fText('Wheel to zoom (letters, then codons, as you go in), drag to pan, double click to reset.' + (findingGenes.length > 1 ? '  Other finding genes: ' + findingGenes.filter(function (q) { return q !== gene; }).slice(0, 8).join(', ') : ''), left, T0 + 38);
    // the 3D buttons, under the gene model where the eye ends up (also easier to reach in VR)
    var by = Math.min(g.cH - 30, endY + 8);
    var btn = function (label, bx, fn) {
      g.setText('white', 13, 'Helvetica, Arial, sans-serif', 'center', 'middle');
      var w2 = g.getTextW(label) + 26; ctx.fillStyle = 'rgba(90,210,190,0.18)'; ctx.fillRect(bx, by, w2, 26); ctx.strokeStyle = 'rgb(90,210,190)'; ctx.lineWidth = 1; ctx.strokeRect(bx, by, w2, 26);
      g.setText('white', 13, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(label, bx + w2 / 2, by + 13);
      if (g.mX > bx && g.mX < bx + w2 && g.mY > by && g.mY < by + 26) { g.setCursor('pointer'); if (g.MOUSE_UP_FAST) fn(); }
      return bx + w2 + 10;
    };
    var bx2 = btn('3D protein: normal and variant side by side', left, function () { G.app.openProtein3D(gene, false); });
    btn('Protein room', bx2, function () { G.app.openProtein3D(gene, true); });
    if (over) view.drawTooltip(g, over);
  };

  G.proteinChange = proteinChange;
  G.ProteinView = ProteinView;
})(globalThis.G = globalThis.G || {});
