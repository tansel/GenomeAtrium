/*
 * Mito: the mitochondrial genome as the usual circular map. At 16,569 bp it
 * is about 1/190,000 of the nuclear genome, so in the genome-wide views it
 * is less than a pixel wide; here it fills the screen.
 *
 *  - the ring runs clockwise from position 1 at the top (rCRS, the GRCh38
 *    chrM); genes sit outside the ring when on the heavy strand (+) and
 *    inside when on the light strand (-): protein-coding genes, rRNAs and
 *    tRNAs (one-letter amino acid) from GENCODE, and the control region
 *    (D-loop, 16,024 to 576) that holds the origin of replication;
 *  - the sample's variants as stems out from the ring, their length the
 *    heteroplasmy (the share of reads carrying the change, from FORMAT AD):
 *    near 100% is homoplasmic, a lower share means a mixture of mtDNA copies;
 *  - ClinVar P/LP mitochondrial entries as ticks inside the ring;
 *  - with a gVCF, a thin ring for calls: confidently called, low depth,
 *    or not called (kept apart, Asclepius D4).
 * Hover reads a gene or variant; a click on a protein-coding gene opens its
 * Protein view. Positions are only comparable on rCRS (16,569 bp): an hg19
 * chrM (16,571 bp, a different sequence) is refused.
 */
(function (G) {
  var LEN = 16569, TAU = Math.PI * 2;
  var TYPE_COL = { protein_coding: 'rgb(90,210,190)', Mt_rRNA: 'rgb(255,120,170)', Mt_tRNA: 'rgb(240,200,80)' };
  var AA1 = { F: 1, V: 1, L: 1, I: 1, Q: 1, M: 1, W: 1, A: 1, N: 1, C: 1, Y: 1, S: 1, D: 1, K: 1, G: 1, R: 1, H: 1, E: 1, T: 1, P: 1 };

  function MitoView() {}

  MitoView.prototype.contig = function (d) {
    var c = d && d.genome && d.genome.get('MT');
    return c && c.length === LEN ? c : null;
  };

  MitoView.prototype.genes = function () {
    if (this._genes) return this._genes;
    var t = G.app.view.genes, out = [];
    if (!t || !t.byContig) return out;
    var mt = t.byContig.MT; // keyed by normName
    (mt ? mt.list : []).forEach(function (g) { out.push(g); });
    if (out.length) this._genes = out;
    return out;
  };

  MitoView.prototype.clinvarSites = function () {
    var cv = G.app.clinvar;
    if (!cv) return [];
    if (this._cv && this._cvFor === cv) return this._cv;
    var out = [];
    cv.index.forEach(function (list) { list.forEach(function (e) { if (G.genome.normName(e.chrom) === 'MT') out.push(e); }); });
    this._cv = out; this._cvFor = cv;
    return out;
  };

  MitoView.prototype.draw = function (g) {
    var ctx = g.context, view = G.app.view, d = view.data, self = this;
    var msg = function (t) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(t, g.cX, g.cY); };
    var c = this.contig(d);
    if (!c) return msg(d && d.genome && d.genome.get('MT') ? 'This file\'s chrM is not the rCRS (16,569 bp; hg19 uses a different chrM): positions cannot be compared.' : 'This file has no mitochondrial contig (chrM or MT).');
    var R = Math.min(g.cW * 0.36, (g.cH - 230) * 0.33), cx = g.cX + 60, cy = 140 + (g.cH - 170) / 2; // room for the outer labels and stems
    var ang = function (p) { return -Math.PI / 2 + TAU * (p - 1) / LEN; };
    var pt = function (p, r) { var a = ang(p); return [cx + Math.cos(a) * r, cy + Math.sin(a) * r]; };
    var hover = null, best = 1e9, mx = g.mX, my = g.mY;
    var mouseAng = Math.atan2(my - cy, mx - cx), mouseR = Math.hypot(mx - cx, my - cy);
    var mousePos = Math.round(((mouseAng + Math.PI / 2 + TAU) % TAU) / TAU * LEN) + 1;

    // calls ring (gVCF): called, low depth, not called
    var tr = d.tracks && d.tracks[c.key];
    if (d.isGvcf && tr) {
      var N = 360;
      for (var i = 0; i < N; i++) {
        var a0 = 1 + LEN * i / N, a1 = 1 + LEN * (i + 1) / N, bs = tr.callable.binSize;
        var b0 = Math.floor((a0 - 1) / bs), b1 = Math.floor((a1 - 1) / bs), call = 0, low = 0, n = 0;
        for (var b = b0; b <= b1; b++) { call += tr.callable.levels[0][b] || 0; low += tr.lowdp.levels[0][b] || 0; n++; }
        call /= n; low /= n;
        ctx.strokeStyle = low > 0.05 ? 'rgba(230,170,40,0.9)' : call > 0.5 ? 'rgba(255,255,255,0.35)' : 'rgba(255,80,80,0.8)';
        ctx.lineWidth = 3; ctx.beginPath(); ctx.arc(cx, cy, R - 30, ang(a0), ang(a1)); ctx.stroke();
      }
    }
    // the ring and the control region
    ctx.strokeStyle = 'rgba(255,255,255,0.5)'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TAU); ctx.stroke();
    ctx.strokeStyle = 'rgba(160,160,180,0.6)'; ctx.lineWidth = 10; ctx.beginPath(); ctx.arc(cx, cy, R, ang(16024), ang(LEN + 576)); ctx.stroke();
    var dl = pt(16024 + (LEN - 16024 + 576) / 2 - LEN, R + 14);
    g.setText('rgba(200,200,220,0.8)', 11, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText('D-loop', dl[0], dl[1] - 20);
    if (mouseR > R - 8 && mouseR < R + 8 && (mousePos >= 16024 || mousePos <= 576)) hover = ['Control region (D-loop), 16,024 to 576', 'holds the origin of heavy-strand replication and the promoters'];
    // position ticks every 1,000
    g.setText('rgba(255,255,255,0.4)', 10, 'Helvetica, Arial, sans-serif', 'center', 'middle');
    for (var p = 1000; p < LEN; p += 1000) { var q = pt(p, R - 14); g.fText((p / 1000) + 'k', q[0], q[1]); }
    // genes: outside the ring on the heavy strand, inside on the light strand
    this.genes().forEach(function (gn) {
      var out = gn.strand > 0, r = out ? R + 16 : R - 16, w = gn.type === 'Mt_tRNA' ? 8 : 14;
      ctx.strokeStyle = TYPE_COL[gn.type] || 'rgb(180,180,180)'; ctx.lineWidth = w;
      ctx.beginPath(); ctx.arc(cx, cy, r, ang(gn.start), ang(gn.end + 1)); ctx.stroke();
      var mid = (gn.start + gn.end) / 2, name = gn.name.replace(/^MT-/, '');
      if (gn.type === 'Mt_tRNA') name = AA1[name.slice(1, 2)] ? name.slice(1) : name.slice(1);
      var lp = pt(mid, out ? R + (gn.type === 'Mt_tRNA' ? 32 : 40) : R - (gn.type === 'Mt_tRNA' ? 32 : 42));
      g.setText(gn.type === 'Mt_tRNA' ? 'rgba(240,200,80,0.9)' : 'white', gn.type === 'Mt_tRNA' ? 10 : 12, 'Helvetica, Arial, sans-serif', 'center', 'middle');
      g.fText(name, lp[0], lp[1]);
      if (Math.abs(mouseR - r) < w / 2 + 3 && mousePos >= gn.start && mousePos <= gn.end) {
        hover = [gn.name + '  ' + { protein_coding: 'protein-coding', Mt_rRNA: 'ribosomal RNA', Mt_tRNA: 'transfer RNA' }[gn.type] + ', ' + (out ? 'heavy' : 'light') + ' strand', 'chrM:' + gn.start.toLocaleString() + '-' + gn.end.toLocaleString() + ' (' + (gn.end - gn.start + 1) + ' bp)'];
        if (gn.type === 'protein_coding') {
          hover.push('click: its protein'); g.setCursor('pointer');
          if (g.MOUSE_UP_FAST) { G.app.focusGene(gn.name, false); G.app.setMode('protein'); }
        }
      }
    });
    // ClinVar P/LP entries: ticks inside the ring
    this.clinvarSites().forEach(function (e) {
      var a = pt(e.pos, R - 44), b = pt(e.pos, R - 54);
      ctx.strokeStyle = 'rgba(255,80,80,0.6)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(a[0], a[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
    });
    // the sample's variants: stems out from the ring, length = heteroplasmy
    var col = d.variants && d.variants[c.key], nVar = 0, nHet = 0, Z = G.vcf.Z;
    if (col) for (var j = 0; j < col.n; j++) {
      var z = col.zygOf(j, 0);
      if (z !== Z.HET && z !== Z.HOM) continue;
      nVar++;
      var f = col.bafAt(j), het = f !== null && f < 0.95, len = 20 + 70 * (f === null ? 1 : f);
      if (het && f >= 0.03) nHet++;
      var p0 = pt(col.pos[j], R + 52), p1 = pt(col.pos[j], R + 52 + len);
      ctx.strokeStyle = het ? 'rgb(255,160,60)' : 'rgb(120,180,255)'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.moveTo(p0[0], p0[1]); ctx.lineTo(p1[0], p1[1]); ctx.stroke();
      ctx.fillStyle = ctx.strokeStyle; ctx.beginPath(); ctx.arc(p1[0], p1[1], 3, 0, TAU); ctx.fill();
      var dd = Math.hypot(mx - p1[0], my - p1[1]);
      if (dd < 8 && dd < best) {
        best = dd;
        var al = col.alleles(j), gene = self.genes().find(function (x) { return col.pos[j] >= x.start && col.pos[j] <= x.end; });
        var cvHits = G.app.clinvar && al ? G.app.clinvar.index.get(G.clinvar.variantKey('MT', col.pos[j], al.ref, al.alts[0])) || [] : [];
        hover = ['chrM:' + col.pos[j].toLocaleString() + ' ' + (al ? al.ref + '>' + al.alts.join(',') : '') + (gene ? '  in ' + gene.name : '  between genes'),
          'heteroplasmy ' + (f === null ? 'unknown (no AD)' : Math.round(100 * f) + '%' + (het ? ' (a mixture of mtDNA copies)' : ' (homoplasmic)')) + (col.dp[j] ? ', ' + col.dp[j] + ' reads' : ''),
          cvHits.length ? 'ClinVar: ' + cvHits[0].sig + ', ' + (cvHits[0].pheno || '').split('|')[0] : ''].filter(Boolean);
      }
    }
    // title and key
    g.setText('white', 15, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('Mitochondrial genome (chrM, 16,569 bp, rCRS)', 20, 150);
    g.setText('rgba(255,255,255,0.6)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText(nVar + ' variants in this sample, ' + nHet + ' heteroplasmic (3 to 95% of reads); ' + this.clinvarSites().length + ' ClinVar P/LP sites (red ticks inside).', 20, 172);
    g.fText('Genes outside the ring: heavy strand; inside: light strand. Teal protein-coding, pink rRNA, yellow tRNA (one-letter amino acid).', 20, 188);
    g.fText('Stems: blue homoplasmic, orange heteroplasmic, length = share of reads.' + (d.isGvcf ? ' Inner ring: called (white), low depth (amber), not called (red).' : ''), 20, 204);
    if (!col || !nVar) g.fText(d.format === 'vcf' ? 'No mitochondrial calls in this file: many pipelines skip chrM by default.' : '', 20, 222);
    if (hover) view.drawTooltip(g, hover);
  };

  G.MitoView = MitoView;
})(globalThis.G = globalThis.G || {});
