/*
 * Pileup in Arcs (reads.js and fasta.js do the reading): when the view is
 * zoomed to 4 kb or less on one chromosome, the band area shows
 *  - the reference bases (letters when they fit, coloured bars otherwise);
 *  - with reads attached: a coverage histogram, coloured by base where a
 *    non-reference base passes 20% of the depth, and the reads packed into
 *    rows: forward reads pink, reverse blue, mismatches in base colours
 *    (A green, C blue, G orange, T red), deletions as dark lines, insertions
 *    as purple marks. High depth is downsampled for drawing (the counts use
 *    every read fetched); very deep regions are sampled and say so.
 * Regions are fetched after the view settles, with a margin either side.
 * Without reads, only a loaded FASTA gives reference letters (at base-level
 * zoom, above the axis): Ensembl is asked only when reads need it, since a
 * request there tells it which region is on screen.
 */
(function (G) {
  var LIMIT = 4000, BASE = { A: 'rgb(0,190,0)', C: 'rgb(50,120,255)', G: 'rgb(230,150,30)', T: 'rgb(230,40,40)', N: 'rgb(150,150,150)' };

  var V = G.View.prototype;

  // The single-chromosome window on screen, if zoomed in enough: {key, chrom, start, end}.
  V.pileupWindow = function () {
    if (!this.visibleSpan || this.zoom <= 1.01) return null;
    var sp = this.visibleSpan();
    if (!sp || sp.end - sp.start > LIMIT) return null;
    var c = this.data.genome.get(sp.chrom);
    return c ? { key: c.key, chrom: c.name, start: Math.max(1, sp.start), end: Math.min(c.length, sp.end) } : null;
  };

  // Fetch the reference (and reads, when attached) for the window, once the view settles.
  V.pileupWant = function (w) {
    var st = this.pileupState || (this.pileupState = {}), self = this;
    if (st.have && st.have.key === w.key && st.have.start <= w.start && st.have.end >= w.end) return;
    if (st.loading && st.loading.key === w.key && st.loading.start <= w.start && st.loading.end >= w.end) return;
    var stamp = w.key + ':' + w.start + ':' + w.end;
    if (st.waitFor !== stamp) { st.waitFor = stamp; st.waitSince = Date.now(); return; }
    if (Date.now() - st.waitSince < 350) return;
    var span = w.end - w.start + 1, c = this.data.genome.get(w.chrom);
    var want = { key: w.key, chrom: w.chrom, start: Math.max(1, w.start - Math.round(span / 2)), end: Math.min(c.length, w.end + Math.round(span / 2)) };
    st.loading = want; st.error = null;
    var build = this.data.build, reads = this.reads;
    Promise.all([
      G.refseq.get(build, want.chrom, want.start, want.end).catch(function (e) { return null; }),
      reads ? reads.fetch(want.chrom, want.start, want.end) : Promise.resolve(null)
    ]).then(function (res) {
      if (st.loading !== want) return;
      st.loading = null; st.have = want; st.ref = res[0]; st.r = res[1];
      st.pile = res[1] ? G.reads.pileup(res[1], res[0]) : null;
    }, function (err) {
      if (st.loading !== want) return;
      st.loading = null; st.error = err.message; console.error(err);
    });
  };

  // No reads: the bases of a loaded FASTA, as letters just above the axis, at base-level zoom.
  V.drawRefOnly = function (g, axisY, w) {
    var f = G.refseq.fasta, s = this.segByKey[w.key];
    if (!f || G.refseq.source(this.data.build) !== 'fasta' || !s) return;
    var sc = this.segScreen(s), ppb = sc.w / s.contig.length;
    if (ppb < 8) return;
    this.pileupWant(w);
    var st = this.pileupState || {}, have = st.have && st.have.key === w.key ? st.have : null;
    if (!have || !st.ref) return;
    for (var p = Math.max(have.start, w.start - 1); p <= Math.min(have.end, w.end + 1); p++) {
      var b = st.ref[p - have.start], x = sc.x + (p - 1) / s.contig.length * sc.w;
      g.setText(BASE[b] || '#999', Math.min(14, ppb * 0.9), 'Helvetica, Arial, sans-serif', 'center', 'bottom'); g.fText(b, x + ppb / 2, axisY - 2);
    }
  };

  // Draws the pileup in the band area; returns true when it did (so the bands are skipped).
  V.drawPileup = function (g, axisY) {
    var w = this.pileupWindow();
    if (!w) return false;
    if (!this.reads) { this.drawRefOnly(g, axisY, w); return false; }
    if (!G.refseq.source(this.data.build)) return false;
    this.pileupWant(w);
    var st = this.pileupState || {}, ctx = g.context, s = this.segByKey[w.key];
    if (!s) return false;
    var sc = this.segScreen(s), L = s.contig.length, ppb = sc.w / L;
    var xOf = function (p) { return sc.x + (p - 1) / L * sc.w; };
    var top = axisY + (this.genesOn() ? 66 : 46), H = g.cH - top - 62, y = top, self = this;
    var have = st.have && st.have.key === w.key ? st.have : null, hover = null;
    // reference bases
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
    g.fText('reference' + (G.refseq.fasta && G.refseq.source(this.data.build) === 'fasta' ? ' (' + G.refseq.fasta.name + ')' : ' (Ensembl)'), 12, y + 1);
    y += 4;
    if (have && st.ref) {
      for (var p = Math.max(have.start, w.start - 1); p <= Math.min(have.end, w.end + 1); p++) {
        var b = st.ref[p - have.start], x = xOf(p);
        if (ppb >= 8) { g.setText(BASE[b] || '#999', Math.min(14, ppb * 0.9), 'Helvetica, Arial, sans-serif', 'center', 'top'); g.fText(b, x + ppb / 2, y); }
        else { ctx.fillStyle = BASE[b] || '#999'; ctx.fillRect(x, y + 2, Math.max(1, ppb), 8); }
      }
    }
    y += 16;
    var pile = have ? st.pile : null;
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
    var r = have && st.r, note = st.error ? 'reads: ' + st.error : st.loading ? 'loading reads...' : r ? r.total.toLocaleString() + ' reads' + (r.downsampled ? ', ' + r.reads.length.toLocaleString() + ' drawn' : '') + (r.sampled ? '; very deep, so sampled: base shares hold, depth does not' : '') : '';
    g.fText('reads: ' + this.reads.name + (note ? '  (' + note + ')' : ''), 12, y + 12);
    y += 16;
    if (!pile) return true;
    // coverage
    var CH = 40, maxD = 1;
    for (p = w.start; p <= w.end; p++) { var cv = pile.cov[p - pile.start]; if (cv && cv.depth > maxD) maxD = cv.depth; }
    for (p = Math.max(pile.start, w.start - 1); p <= Math.min(pile.end, w.end + 1); p++) {
      var cvp = pile.cov[p - pile.start], rb = pile.ref ? pile.ref[p - pile.start] : null;
      if (!cvp || !cvp.depth) continue;
      var x0 = xOf(p), bw = Math.max(1, ppb - (ppb > 3 ? 1 : 0)), h = CH * cvp.depth / maxD, base = y + CH;
      var alt = 0; ['A', 'C', 'G', 'T'].forEach(function (k) { if (k !== rb) alt += cvp[k]; }); alt += cvp.del;
      if (rb && alt / cvp.depth > 0.2) { // stacked by base
        var yy = base;
        ['A', 'C', 'G', 'T'].forEach(function (k) { var hh = CH * cvp[k] / maxD; ctx.fillStyle = BASE[k]; ctx.fillRect(x0, yy - hh, bw, hh); yy -= hh; });
        if (cvp.del) { var hd = CH * cvp.del / maxD; ctx.fillStyle = 'rgb(70,70,70)'; ctx.fillRect(x0, yy - hd, bw, hd); }
      } else { ctx.fillStyle = 'rgba(180,180,180,0.7)'; ctx.fillRect(x0, base - h, bw, h); }
      if (g.mX >= x0 && g.mX < x0 + Math.max(1, ppb) && g.mY >= y && g.mY <= base) hover = [w.chrom + ':' + p.toLocaleString() + '  reference ' + (rb || '?'), 'depth ' + cvp.depth + ': A ' + cvp.A + ', C ' + cvp.C + ', G ' + cvp.G + ', T ' + cvp.T + (cvp.del ? ', deleted ' + cvp.del : '') + (cvp.ins ? ', insertion after ' + cvp.ins : '')];
    }
    g.setText('rgba(255,255,255,0.4)', 10, 'Helvetica, Arial, sans-serif', 'right', 'top'); g.fText('max depth ' + maxD, g.cW - 12, y);
    y += CH + 6;
    // reads
    var rowH = ppb >= 8 ? 9 : ppb >= 2 ? 5 : 3, gap = 1, maxRows = Math.max(1, Math.floor((top + H - y) / (rowH + gap)));
    pile.rows.slice(0, maxRows).forEach(function (row, ri) {
      var ry = y + ri * (rowH + gap);
      row.forEach(function (rd) {
        if (rd.end < w.start - 1 || rd.start > w.end + 1) return;
        var x0 = xOf(rd.start), x1 = xOf(rd.end + 1);
        ctx.fillStyle = rd.strand > 0 ? 'rgba(215,165,165,' + (rd.mapq ? 0.6 : 0.25) + ')' : 'rgba(165,180,220,' + (rd.mapq ? 0.6 : 0.25) + ')';
        ctx.fillRect(x0, ry, Math.max(1, x1 - x0), rowH);
        rd.marks.forEach(function (m) {
          var mx = xOf(m[1]);
          if (m[0] === 'x') { ctx.fillStyle = BASE[m[2]] || '#999'; ctx.fillRect(mx, ry, Math.max(1, ppb), rowH); }
          else if (m[0] === 'd') { ctx.fillStyle = 'rgb(30,30,34)'; ctx.fillRect(mx, ry, Math.max(1, ppb), rowH); ctx.fillStyle = 'rgba(255,255,255,0.6)'; ctx.fillRect(mx, ry + rowH / 2, Math.max(1, ppb), 1); }
          else { ctx.fillStyle = 'rgb(170,60,230)'; ctx.fillRect(mx - 1, ry - 1, 2, rowH + 2); }
        });
        if (g.mY >= ry && g.mY < ry + rowH && g.mX >= x0 && g.mX < x1) hover = [rd.name, w.chrom + ':' + rd.start.toLocaleString() + '-' + rd.end.toLocaleString() + (rd.strand > 0 ? ' forward' : ' reverse') + ', MAPQ ' + rd.mapq,
          'CIGAR ' + rd.cigar.map(function (c) { return c[1] + c[0]; }).join('').slice(0, 60) + ', ' + rd.marks.filter(function (m) { return m[0] === 'x'; }).length + ' mismatches'];
      });
    });
    if (pile.rows.length > maxRows) { g.setText('rgba(255,255,255,0.45)', 10, 'Helvetica, Arial, sans-serif', 'right', 'bottom'); g.fText('+' + (pile.rows.length - maxRows) + ' more rows (zoom the page out, or a taller window)', g.cW - 12, top + H); }
    if (hover) this.hover = { lines: hover };
    return true;
  };
})(globalThis.G = globalThis.G || {});
