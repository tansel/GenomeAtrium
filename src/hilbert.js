/*
 * Hilbert map: the whole genome folded into a square along a Hilbert curve,
 * which keeps neighbours on the genome next to each other in 2D (Anders
 * 2009; HilbertCurve, Gu et al. 2016). At order 9 the square is 512 x 512
 * cells, each about 12 kb of a human genome: a whole genome at kilobase
 * detail in one picture.
 *
 * Layers: callability (gVCF), variant density, het fraction, ClinVar P/LP
 * density, depth (BAM). Chromosome outlines, findings and the shared
 * selection are drawn on top. Wheel zooms, drag pans, click flies to the
 * cell in Arcs.
 */
(function (G) {
  var ORDER = 9, N = 1 << ORDER, CELLS = N * N;

  // Hilbert index d -> [x, y] (the classic iterative form).
  function d2xy(d) {
    var x = 0, y = 0, t = d;
    for (var s = 1; s < N; s *= 2) {
      var rx = 1 & (t >> 1), ry = 1 & (t ^ rx);
      if (ry === 0) { if (rx === 1) { x = s - 1 - x; y = s - 1 - y; } var tmp = x; x = y; y = tmp; }
      x += s * rx; y += s * ry; t >>= 2;
    }
    return [x, y];
  }
  var XY = null; // precomputed cell -> x, y
  function table() {
    if (XY) return XY;
    XY = { x: new Uint16Array(CELLS), y: new Uint16Array(CELLS), inv: new Int32Array(CELLS) };
    for (var d = 0; d < CELLS; d++) { var p = d2xy(d); XY.x[d] = p[0]; XY.y[d] = p[1]; XY.inv[p[1] * N + p[0]] = d; }
    return XY;
  }

  var LAYERS = {
    callable: { label: 'callability', need: function (d) { return d.format === 'vcf' && d.isGvcf; } },
    density: { label: 'variant density', need: function (d) { return d.format === 'vcf'; } },
    het: { label: 'het fraction', need: function (d) { return d.format === 'vcf'; } },
    clinvar: { label: 'ClinVar P/LP sites', need: function (d) { return d.format === 'vcf' && d.clinvarStatus === 'matched'; } },
    depth: { label: 'depth', need: function (d) { return d.format === 'bam'; } }
  };

  function Hilbert() { this.zoom = 1; this.ox = 0; this.oy = 0; this.layer = null; this.img = null; }

  Hilbert.prototype.setData = function (d) {
    this.data = d; this.img = null; this.zoom = 1; this.ox = this.oy = 0;
    if (!d || !d.genome || !d.genome.contigs.length) return;
    var avail = Object.keys(LAYERS).filter(function (k) { return LAYERS[k].need(d); });
    if (!this.layer || avail.indexOf(this.layer) < 0) this.layer = avail[0];
    this.layers = avail;
    this.build();
  };

  Hilbert.prototype.setLayer = function (k) { this.layer = k; this.build(); };

  // Genome position <-> cell. The concatenated genome is spread over all cells.
  Hilbert.prototype.build = function () {
    var d = this.data, t = table();
    var total = d.genome.totalLength(), bpc = total / CELLS, offsets = {}, acc = 0;
    d.genome.contigs.forEach(function (c) { offsets[c.key] = acc; acc += c.length; });
    this.bpc = bpc; this.offsets = offsets;
    var val = new Float32Array(CELLS), owner = new Int16Array(CELLS).fill(-1), layer = this.layer, cv = G.app && G.app.clinvar;
    var self = this;
    d.genome.contigs.forEach(function (c, ci) {
      var tr = d.tracks[c.key];
      var d0 = Math.floor(offsets[c.key] / bpc), d1 = Math.min(CELLS - 1, Math.floor((offsets[c.key] + c.length - 1) / bpc));
      var lv = {}, sample = function (name, a, b, kind) { var L = lv[name] || (lv[name] = tr[name].levelFor(bpc)); return self.sample(L, a, b, kind); };
      var cvt = layer === 'clinvar' && cv ? cv.trackFor(c.key, c.length, d.binSize) : null, cvL = cvt ? cvt.levelFor(bpc) : null;
      for (var k = d0; k <= d1; k++) {
        owner[k] = ci;
        if (!tr) continue;
        var a = Math.max(0, k * bpc - offsets[c.key]), b = Math.min(c.length, (k + 1) * bpc - offsets[c.key]);
        // callability: 2 = mostly low depth, else the called fraction (0..1)
        if (layer === 'callable') { var low = sample('lowdp', a, b, 'mean'); val[k] = low > 0.3 ? 2 : Math.min(1, sample('callable', a, b, 'mean')); }
        else if (layer === 'density') val[k] = sample('snv', a, b, 'sum') + sample('indel', a, b, 'sum');
        else if (layer === 'het') { var he = sample('het', a, b, 'sum'), ho = sample('hom', a, b, 'sum'); val[k] = he + ho >= 3 ? he / (he + ho) : -1; }
        else if (layer === 'clinvar') val[k] = cvL ? self.sample(cvL, a, b, 'sum') : 0;
        else if (layer === 'depth') val[k] = sample('depth', a, b, 'mean');
      }
    });
    // colour scale: 98th percentile for counts
    var sorted = Array.from(val).filter(function (v) { return v > 0; }).sort(function (x, y) { return x - y; });
    var hi = sorted.length ? sorted[Math.floor(sorted.length * 0.98)] : 1;
    var cvs = document.createElement('canvas'); cvs.width = N; cvs.height = N;
    var ctx = cvs.getContext('2d'), im = ctx.createImageData(N, N), px = im.data;
    for (var k2 = 0; k2 < CELLS; k2++) {
      var o = (t.y[k2] * N + t.x[k2]) * 4, v = val[k2], col;
      if (owner[k2] < 0) col = [20, 20, 20];
      else if (layer === 'callable') col = v >= 2 ? [230, 170, 40] : v > 0.5 ? [70 + 120 * Math.min(1, v), 70 + 120 * Math.min(1, v), 90 + 120 * Math.min(1, v)] : [150, 40, 40];
      else if (layer === 'het') col = v < 0 ? [35, 35, 35] : [40 + 150 * v, 30 + 110 * v, 80 + 175 * v];
      else { var f = Math.sqrt(Math.min(1, v / (hi || 1))); col = layer === 'clinvar' ? [40 + 215 * f, 30 + 50 * f, 30 + 50 * f] : [30 + 60 * f, 40 + 170 * f, 60 + 130 * f]; }
      px[o] = col[0]; px[o + 1] = col[1]; px[o + 2] = col[2]; px[o + 3] = 255;
    }
    // chromosome outlines: a cell next to a cell of another chromosome
    for (var y = 0; y < N; y++) for (var x = 0; x < N; x++) {
      var dd = t.inv[y * N + x], own = owner[dd];
      if (own < 0) continue;
      var edge = (x + 1 < N && owner[t.inv[y * N + x + 1]] !== own) || (y + 1 < N && owner[t.inv[(y + 1) * N + x]] !== own);
      if (edge) { var q = (y * N + x) * 4; px[q] = px[q + 1] = px[q + 2] = 235; }
    }
    ctx.putImageData(im, 0, 0);
    this.img = cvs; this.owner = owner; this.val = val; this.hi = hi;
    // label points: mean x/y of each chromosome's cells
    var sx = {}, sy = {}, n = {};
    for (var k3 = 0; k3 < CELLS; k3++) { var ow = owner[k3]; if (ow < 0) continue; sx[ow] = (sx[ow] || 0) + t.x[k3]; sy[ow] = (sy[ow] || 0) + t.y[k3]; n[ow] = (n[ow] || 0) + 1; }
    this.labels = Object.keys(n).map(function (ci) { return { name: d.genome.contigs[ci].name, x: sx[ci] / n[ci], y: sy[ci] / n[ci], n: n[ci] }; });
  };

  Hilbert.prototype.sample = function (lv, a, b, kind) {
    var bs = lv.binSize, arr = lv.data, i0 = Math.max(0, Math.floor(a / bs)), i1 = Math.min(arr.length - 1, Math.floor((b - 1e-6) / bs));
    if (i1 < i0) return 0;
    var s = 0;
    if (kind === 'mean') { for (var i = i0; i <= i1; i++) s += arr[i]; return s / (i1 - i0 + 1); }
    for (var j = i0; j <= i1; j++) { var lo = Math.max(a, j * bs), hi = Math.min(b, (j + 1) * bs); if (hi > lo) s += arr[j] * (hi - lo) / bs; }
    return s;
  };

  Hilbert.prototype.cellOf = function (chrom, pos) {
    var k = G.genome.normName(chrom);
    if (this.offsets[k] === undefined) return -1;
    return Math.min(CELLS - 1, Math.floor((this.offsets[k] + pos - 1) / this.bpc));
  };

  Hilbert.prototype.draw = function (g) {
    var ctx = g.context;
    if (!this.img) { g.setText('rgba(255,255,255,0.6)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText('Open a file first.', g.cX, g.cY); return; }
    var inWin = G.app.view.mode === 'atrium'; // the Atrium's window: no page chrome, fill the height
    var side = inWin ? Math.min(g.cW - 40, g.cH - 50) : Math.min(g.cW - 80, g.cH - 200), base = side / N;
    if (g.WHEEL_CHANGE || this.pendingWheel) { /* handled by the wheel listener */ }
    var s = base * this.zoom, left = (g.cW - side) / 2, top = inWin ? 12 : 130;
    if (g.MOUSE_PRESSED) { this.ox += g.DX_MOUSE; this.oy += g.DY_MOUSE; }
    this.ox = Math.min(0, Math.max(this.ox, side - N * s)); this.oy = Math.min(0, Math.max(this.oy, side - N * s));
    var X = left + this.ox, Y = top + this.oy, t = table(), self = this;
    this.geom = { left: left, top: top, side: side, s: s, X: X, Y: Y, base: base };
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, side, side); ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(this.img, X, Y, N * s, N * s);
    // selection and findings on the curve
    var sel = G.app.view.selection;
    if (sel) {
      ctx.fillStyle = 'rgba(120,200,255,0.8)';
      var c0 = this.cellOf(sel.chrom, sel.start), c1 = this.cellOf(sel.chrom, sel.end);
      for (var k = c0; k >= 0 && k <= c1; k++) ctx.fillRect(X + t.x[k] * s, Y + t.y[k] * s, Math.max(1, s), Math.max(1, s));
    }
    (G.app.view.findings || []).forEach(function (f) {
      var k = self.cellOf(f.chrom, f.pos);
      if (k < 0) return;
      var cx = X + (t.x[k] + 0.5) * s, cy = Y + (t.y[k] + 0.5) * s;
      ctx.strokeStyle = f.classification === 'Pathogenic' ? 'rgb(255,70,70)' : 'rgb(255,160,60)'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(cx, cy, 5, 0, 2 * Math.PI); ctx.stroke();
    });
    g.setText('rgba(255,255,255,0.85)', 12, 'Helvetica, Arial, sans-serif', 'center', 'middle');
    this.labels.forEach(function (l) {
      if (l.n * s * s < 900) return;
      ctx.fillStyle = 'rgba(0,0,0,0.55)'; var w = g.getTextW(l.name) + 8;
      ctx.fillRect(X + l.x * s - w / 2, Y + l.y * s - 9, w, 18);
      g.setText('white', 12, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(l.name, X + l.x * s, Y + l.y * s);
    });
    ctx.restore();
    ctx.strokeStyle = 'rgba(255,255,255,0.25)'; ctx.strokeRect(left - 0.5, top - 0.5, side + 1, side + 1);

    // hover
    var cx0 = Math.floor((g.mX - X) / s), cy0 = Math.floor((g.mY - Y) / s);
    if (g.mX >= left && g.mX <= left + side && g.mY >= top && g.mY <= top + side && cx0 >= 0 && cy0 >= 0 && cx0 < N && cy0 < N) {
      var dd = t.inv[cy0 * N + cx0], ow = this.owner[dd];
      if (ow >= 0) {
        var c = this.data.genome.contigs[ow], pos = Math.max(1, Math.round(dd * this.bpc - this.offsets[c.key] + this.bpc / 2));
        var v = this.val[dd];
        var vt = this.layer === 'callable' ? (v >= 2 ? 'low depth' : v > 0.5 ? 'called ' + Math.round(Math.min(1, v) * 100) + '%' : 'not called') :
          this.layer === 'het' ? (v < 0 ? 'too few genotypes' : 'het fraction ' + v.toFixed(2)) :
          this.layer === 'depth' ? 'mean depth ' + v.toFixed(1) + 'x' : v.toFixed(1) + (this.layer === 'clinvar' ? ' ClinVar P/LP sites' : ' variants');
        ctx.strokeStyle = 'white'; ctx.lineWidth = 1.5; ctx.strokeRect(X + cx0 * s - 1, Y + cy0 * s - 1, Math.max(3, s + 2), Math.max(3, s + 2));
        G.app.view.drawTooltip(g, [c.name + ':' + pos.toLocaleString() + '  (cell of ' + G.fmtBp(this.bpc) + ')', vt, 'click to open in Arcs']);
        g.setCursor('pointer');
        if (g.MOUSE_UP_FAST) { G.app.setMode('arcs'); G.app.view.goTo(c.name, Math.max(1, pos - 50000), pos + 50000); }
      }
    }
    // legend and layer buttons are HTML (app.js); a caption here
    g.setText('rgba(255,255,255,0.5)', 11, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText('Hilbert map: the genome folded into a square, neighbours stay neighbours. ' + N + ' x ' + N + ' cells of ' + G.fmtBp(this.bpc) + '. Wheel zooms, drag pans, click opens in Arcs.', g.cX, top + side + 10);
  };

  Hilbert.prototype.install = function (canvas) {
    var self = this;
    canvas.addEventListener('wheel', function (e) {
      if (G.app.view.activeMode() !== 'hilbert') return;
      self.zoomAt(e.offsetX, e.offsetY, Math.exp(-e.deltaY * 0.0018));
    }, { passive: true });
  };

  // Zoom by factor k keeping canvas point (x, y) in place.
  Hilbert.prototype.zoomAt = function (x, y, k) {
    if (!this.geom) return;
    var gm = this.geom, z0 = this.zoom, z1 = Math.max(1, Math.min(64, z0 * k));
    var ux = (x - gm.left - this.ox) / (gm.base * z0), uy = (y - gm.top - this.oy) / (gm.base * z0);
    this.zoom = z1; this.ox = x - gm.left - ux * gm.base * z1; this.oy = y - gm.top - uy * gm.base * z1;
  };

  G.Hilbert = Hilbert; G.HILBERT_LAYERS = LAYERS; G.hilbertD2xy = d2xy;
})(globalThis.G = globalThis.G || {});
