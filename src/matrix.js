/*
 * Matrix mode, after the matrix view of moebio.com/attention: genome
 * windows along both axes, and each cell is how strongly two windows share
 * an unusual profile (strength in src/windows.js: cosine weighted by how far
 * both windows are from the average window).
 *
 * The diagonal is bright only where a window is unusual (an ordinary window
 * has a weak profile, so even its self-link is weak). Bright blocks are regions
 * that behave alike: chrX and chrY in a male, runs of homozygosity, or
 * regions that are all hard to call. SV breakends and discordant read pairs
 * are drawn as circles at the cell of their two ends, and findings as red
 * marks along the edges.
 *
 * Wheel zooms around the cursor, drag pans, double click resets, click a
 * cell to open its column window in Arcs.
 */
(function (G) {
  var MARGIN = 70, TOP = 140;
  var FLOOR = 0.3; // strengths below this are drawn as background, to cut speckle

  // black -> violet -> orange -> pale yellow, for similarity 0 to 1
  var RAMP = [[20, 20, 20], [70, 30, 120], [200, 60, 90], [250, 160, 60], [255, 245, 200]];
  function ramp(t) {
    t = Math.max(0, Math.min(1, t)) * (RAMP.length - 1);
    var i = Math.min(RAMP.length - 2, Math.floor(t)), f = t - i, a = RAMP[i], b = RAMP[i + 1];
    return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
  }

  function Matrix() {
    this.zoom = 1; this.zoomTarget = 1; this.ox = 0; this.oy = 0; this.anchor = [0, 0];
    var self = this;
    G.app.view.g.canvas.addEventListener('wheel', function (e) {
      if (G.app.view.activeMode() !== 'matrix') return;
      self.zoomAt(e.offsetX, e.offsetY, Math.exp(-e.deltaY * 0.0018));
    }, { passive: true });
    G.app.view.g.canvas.addEventListener('dblclick', function () { if (G.app.view.mode === 'matrix') self.reset(); });
  }

  Matrix.prototype.zoomAt = function (x, y, k) {
    this.zoomTarget = Math.max(1, Math.min(60, this.zoomTarget * k));
    this.anchor = [x, y];
  };

  Matrix.prototype.reset = function () { this.zoom = this.zoomTarget = 1; this.ox = this.oy = 0; };

  Matrix.prototype.setData = function (d) {
    this.data = d; this.mdl = null; this.img = null; this.reset();
    var mdl = G.windows.model(d);
    if (!mdl) return;
    var n = mdl.windows.length, m = mdl.m, u = mdl.unit, norms = mdl.norms, ref = mdl.normRef;
    var cv = document.createElement('canvas');
    cv.width = n; cv.height = n;
    var cx = cv.getContext('2d'), im = cx.createImageData(n, n), px = im.data;
    for (var i = 0; i < n; i++) {
      for (var j = i; j < n; j++) {
        var s = 0;
        for (var k = 0; k < m; k++) s += u[i * m + k] * u[j * m + k];
        s *= Math.min(1, Math.min(norms[i], norms[j]) / ref);
        var c = ramp((s - FLOOR) / (1 - FLOOR)); // below the floor stays background
        var p1 = (i * n + j) * 4, p2 = (j * n + i) * 4;
        px[p1] = px[p2] = c[0]; px[p1 + 1] = px[p2 + 1] = c[1]; px[p1 + 2] = px[p2 + 2] = c[2]; px[p1 + 3] = px[p2 + 3] = 255;
      }
    }
    cx.putImageData(im, 0, 0);
    this.img = cv;
    this.mdl = mdl;
    // contig blocks along the axes
    var blocks = [], cur = null;
    mdl.windows.forEach(function (w, i) {
      if (!cur || cur.ci !== w.ci) { cur = { ci: w.ci, name: w.contig.name, from: i, to: i }; blocks.push(cur); }
      cur.to = i;
    });
    this.blocks = blocks;
    // SV and read-pair arcs placed on cells
    var idx = {};
    mdl.windows.forEach(function (w, i) { idx[w.contig.key + ':' + w.j] = i; });
    var cell = function (c, p) { return idx[G.genome.normName(c) + ':' + Math.floor((p - 1) / mdl.win)]; };
    this.marks = (d.arcs || []).map(function (a) {
      var i = cell(a.c0, a.p0), j = cell(a.c1, a.p1);
      // an SV inside one window sits on the diagonal and says nothing here
      return i === undefined || j === undefined || i === j ? null : { i: i, j: j, a: a };
    }).filter(Boolean);
    this.cellOf = cell;
  };

  Matrix.prototype.draw = function (g) {
    var ctx = g.context, mdl = this.mdl, self = this;
    if (!mdl) {
      g.setText('rgba(255,255,255,0.6)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle');
      g.fText('Not enough windows with data for a matrix.', g.cX, g.cY);
      return;
    }
    var n = mdl.windows.length;
    var base = Math.min(g.cW - MARGIN * 2, g.cH - TOP - 56) / n;
    var left = Math.max(MARGIN, (g.cW - base * n) / 2), top = TOP;
    // zoom around the anchor, drag to pan, keep the matrix on screen
    var before = base * this.zoom;
    var ux = (this.anchor[0] - left - this.ox) / before, uy = (this.anchor[1] - top - this.oy) / before;
    this.zoom += (this.zoomTarget - this.zoom) * 0.2;
    var s = base * this.zoom;
    this.ox = this.anchor[0] - left - ux * s; this.oy = this.anchor[1] - top - uy * s;
    if (g.MOUSE_PRESSED) { this.ox += g.DX_MOUSE; this.oy += g.DY_MOUSE; this.zoomTarget = this.zoom; }
    var span = n * s, viewW = base * n;
    this.ox = Math.min(0, Math.max(this.ox, viewW - span)); this.oy = Math.min(0, Math.max(this.oy, viewW - span));
    var X = left + this.ox, Y = top + this.oy;

    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, viewW, viewW); ctx.clip();
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(this.img, X, Y, span, span);
    // contig grid
    ctx.strokeStyle = 'rgba(255,255,255,0.18)'; ctx.lineWidth = 1;
    this.blocks.forEach(function (b) {
      var p = b.from * s;
      ctx.beginPath(); ctx.moveTo(X + p, top); ctx.lineTo(X + p, top + viewW); ctx.moveTo(left, Y + p); ctx.lineTo(left + viewW, Y + p); ctx.stroke();
    });
    // SV / discordant pair cells
    ctx.strokeStyle = 'rgba(120,255,220,0.9)'; ctx.lineWidth = 1.2;
    this.marks.forEach(function (mk) {
      [[mk.i, mk.j], [mk.j, mk.i]].forEach(function (c) {
        ctx.beginPath(); ctx.arc(X + (c[1] + 0.5) * s, Y + (c[0] + 0.5) * s, Math.max(3, s * 0.6), 0, 2 * Math.PI); ctx.stroke();
      });
    });
    // shared selection (alt+drag in Arcs, or a click here or in the Landscape)
    var sel = this.selection;
    if (sel) {
      ctx.fillStyle = 'rgba(120,200,255,0.18)';
      mdl.windows.forEach(function (w, i) {
        if (G.genome.normName(sel.chrom) !== w.contig.key || w.end < sel.start || w.start > sel.end) return;
        ctx.fillRect(X + i * s, top, Math.max(1, s), viewW); ctx.fillRect(left, Y + i * s, viewW, Math.max(1, s));
      });
    }
    ctx.restore();

    // findings along the top and left edges
    ctx.fillStyle = 'rgb(255,70,70)';
    (G.app.view.findings || []).forEach(function (f) {
      var i = self.cellOf(f.chrom, f.pos);
      if (i === undefined) return;
      var p = (i + 0.5) * s;
      if (X + p >= left && X + p <= left + viewW) ctx.fillRect(X + p - 1, top - 8, 2, 6);
      if (Y + p >= top && Y + p <= top + viewW) ctx.fillRect(left - 8, Y + p - 1, 6, 2);
    });

    // contig labels
    g.setText('rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'center', 'bottom');
    this.blocks.forEach(function (b) {
      var c = X + (b.from + b.to + 1) / 2 * s, w = (b.to - b.from + 1) * s;
      if (w < 14 || c < left || c > left + viewW) return;
      g.fText(b.name, c, top - 12);
    });
    g.setText('rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'right', 'middle');
    this.blocks.forEach(function (b) {
      var c = Y + (b.from + b.to + 1) / 2 * s, w = (b.to - b.from + 1) * s;
      if (w < 12 || c < top || c > top + viewW) return;
      g.fText(b.name, left - 12, c);
    });

    // hover a cell: crosshair and both windows
    var j = Math.floor((g.mX - X) / s), i = Math.floor((g.mY - Y) / s);
    if (g.mX >= left && g.mX <= left + viewW && g.mY >= top && g.mY <= top + viewW && i >= 0 && j >= 0 && i < n && j < n) {
      ctx.fillStyle = 'rgba(255,255,255,0.12)';
      ctx.fillRect(left, Y + i * s, viewW, Math.max(1, s)); ctx.fillRect(X + j * s, top, Math.max(1, s), viewW);
      var wi = mdl.windows[i], wj = mdl.windows[j];
      var lines = ['strength ' + G.windows.strength(mdl, i, j).toFixed(2) + '  (cosine ' + G.windows.similarity(mdl, i, j).toFixed(2) + ')',
        'column ' + wj.contig.name + ':' + wj.start.toLocaleString() + '-' + wj.end.toLocaleString(), '  ' + G.windows.describe(mdl, wj),
        'row ' + wi.contig.name + ':' + wi.start.toLocaleString() + '-' + wi.end.toLocaleString(), '  ' + G.windows.describe(mdl, wi)];
      this.marks.forEach(function (mk) { if ((mk.i === i && mk.j === j) || (mk.i === j && mk.j === i)) lines.push(mk.a.label); });
      lines.push('click to open the column window in Arcs');
      G.app.view.drawTooltip(g, lines);
      g.setCursor('pointer');
      if (g.MOUSE_UP_FAST) {
        G.app.select({ chrom: wj.contig.name, start: wj.start, end: wj.end });
        G.app.setMode('arcs'); G.app.view.goTo(wj.contig.name, wj.start, wj.end);
      }
    }

    // legend
    var lx = 12, ly = g.cH - 40;
    for (var t = 0; t <= 100; t++) { var c = ramp(t / 100); ctx.fillStyle = 'rgb(' + c.map(Math.round).join(',') + ')'; ctx.fillRect(lx + t * 1.6, ly, 2, 8); }
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('strength ' + FLOOR + ' to 1: shared profile, weighted by how unusual both windows are   circles: SV and read-pair ends in two windows   red marks: findings', lx, ly + 11);
    g.fText(n.toLocaleString() + ' windows of ' + G.fmtBp(mdl.win) + '; features: ' + mdl.names.join(', '), lx, ly + 24);
  };

  if (typeof mo !== 'undefined' && G.app) G.matrix = new Matrix();
})(globalThis.G = globalThis.G || {});
