/*
 * Landscape mode, after moebio.com/mind: the genome is cut into windows,
 * each window gets a small feature vector, and PCA places the windows in
 * 3D. Consecutive windows of a chromosome are joined, so each chromosome
 * is a thread through the space. Windows that look alike end up close
 * together wherever they sit on the genome.
 *
 * Drawn with Moebio's Engine3D. Drag rotates, wheel zooms, clicking a
 * window opens that region in the arcs view.
 */
(function (G) {
  var CUBE = 260;

  // Eigen decomposition of a small symmetric matrix (Jacobi rotations).
  function eigenSym(M) {
    var n = M.length, A = M.map(function (r) { return r.slice(); });
    var V = A.map(function (_, i) { return A.map(function (_, j) { return i === j ? 1 : 0; }); });
    for (var sweep = 0; sweep < 60; sweep++) {
      var off = 0;
      for (var p = 0; p < n; p++) for (var q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
      if (off < 1e-12) break;
      for (p = 0; p < n; p++) for (q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-15) continue;
        var th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        var t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        var c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (var k = 0; k < n; k++) {
          var akp = A[k][p], akq = A[k][q];
          A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq;
        }
        for (k = 0; k < n; k++) {
          var apk = A[p][k], aqk = A[q][k];
          A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk;
        }
        for (k = 0; k < n; k++) {
          var vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq;
        }
      }
    }
    var vals = A.map(function (r, i) { return r[i]; });
    var order = vals.map(function (_, i) { return i; }).sort(function (a, b) { return vals[b] - vals[a]; });
    return order.map(function (i) { return { value: vals[i], vector: V.map(function (r) { return r[i]; }) }; });
  }

  // PCA on the standardised features of the shared window model.
  function pca(mdl) {
    var names = mdl.names, used = mdl.used, Z = mdl.Z, n = Z.length, m = used.length;
    var C = used.map(function () { return new Array(m).fill(0); });
    Z.forEach(function (z) { for (var a = 0; a < m; a++) for (var b = 0; b < m; b++) C[a][b] += z[a] * z[b] / n; });
    var eig = m ? eigenSym(C) : [];
    while (eig.length < 3) eig.push({ value: 0, vector: new Array(m).fill(0) });
    var totalVar = eig.reduce(function (s, e) { return s + Math.max(0, e.value); }, 0) || 1;
    var coords = Z.map(function (z) {
      return [0, 1, 2].map(function (k) { return z.reduce(function (s, v, i) { return s + v * eig[k].vector[i]; }, 0); });
    });
    // Scale so the cloud fills the cube; clip far outliers so they do not squash the rest.
    [0, 1, 2].forEach(function (k) {
      var vals = coords.map(function (c) { return Math.abs(c[k]); }).sort(function (a, b) { return a - b; });
      var q = vals[Math.floor(vals.length * 0.98)] || 1;
      coords.forEach(function (c) { c[k] = Math.max(-1.3, Math.min(1.3, c[k] / q)) * CUBE; });
    });
    var axes = [0, 1, 2].map(function (k) {
      var load = used.map(function (fi, i) { return { name: names[fi], w: eig[k].vector[i] }; })
        .sort(function (a, b) { return Math.abs(b.w) - Math.abs(a.w); }).slice(0, 2);
      return { share: Math.max(0, eig[k].value) / totalVar, load: load };
    });
    return { coords: coords, axes: axes };
  }

  function contigColor(i, n, a) { return G.contigColor(i, n, a); } // shared with the Arcs view

  function Landscape() {
    this.e3D = new mo.Engine3D({ lens: 900 });
    this.e3D.setAngles(new mo.Point3D(0.35, -0.5, 0));
    this.zoom = 1; this.zoomTarget = 1; this.spin = 0.0025;
    this.data = null;
    var self = this;
    G.app.view.g.canvas.addEventListener('wheel', function (e) {
      if (G.app.view.activeMode() !== '3d') return;
      self.zoomTarget = Math.max(0.3, Math.min(6, self.zoomTarget * Math.exp(-e.deltaY * 0.0015)));
    }, { passive: true });
  }

  Landscape.prototype.setData = function (d) {
    this.data = d;
    this.model = null;
    if (!d || !d.genome || !d.genome.contigs.length) return;
    var t0 = performance.now();
    var mdl = G.windows.model(d);
    if (!mdl) return;
    var p = pca(mdl);
    mdl.windows.forEach(function (w, i) { w.p = p.coords[i]; });
    this.model = { windows: mdl.windows, names: mdl.names, win: mdl.win, axes: p.axes, ms: performance.now() - t0 };
  };

  Landscape.prototype.draw = function (g) {
    var ctx = g.context, m = this.model;
    if (!m) {
      g.setText('rgba(255,255,255,0.6)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle');
      g.fText(this.data ? 'Not enough windows with data for a landscape.' : 'Open a file first.', g.cX, g.cY);
      return;
    }
    var e = this.e3D;
    if (g.MOUSE_PRESSED) e.applyRotation(new mo.Point(g.DX_MOUSE * 0.006, g.DY_MOUSE * 0.006));
    else e.applyRotation(new mo.Point(this.spin, 0));
    this.zoom += (this.zoomTarget - this.zoom) * 0.15;
    var k = Math.min(g.cW, g.cH) / 1000 * this.zoom, cx = g.cX, cy = g.cY + 20;
    var nC = G.app.view.data.genome.contigs.length;

    // cube, as in /mind
    var h = CUBE * 1.3, corners = [], self = this;
    [-h, h].forEach(function (x) { [-h, h].forEach(function (y) { [-h, h].forEach(function (z) { corners.push(e.projectCoordinates(x, y, z)); }); }); });
    ctx.strokeStyle = 'rgba(255,255,255,0.12)'; ctx.lineWidth = 1;
    [[0, 1], [0, 2], [0, 4], [1, 3], [1, 5], [2, 3], [2, 6], [3, 7], [4, 5], [4, 6], [5, 7], [6, 7]].forEach(function (ed) {
      var a = corners[ed[0]], b = corners[ed[1]];
      ctx.beginPath(); ctx.moveTo(cx + a.x * k, cy + a.y * k); ctx.lineTo(cx + b.x * k, cy + b.y * k); ctx.stroke();
    });

    var proj = m.windows.map(function (w) { var q = e.projectCoordinates(w.p[0], w.p[1], w.p[2]); return { w: w, x: cx + q.x * k, y: cy + q.y * k, z: q.z }; });

    // hover: nearest window on screen
    var over = null, best = 14;
    proj.forEach(function (p) { var dd = Math.hypot(p.x - g.mX, p.y - g.mY); if (dd < best) { best = dd; over = p; } });
    var focusCi = over ? over.w.ci : null;

    // threads: consecutive windows of one contig
    for (var i = 1; i < proj.length; i++) {
      var a = proj[i - 1], b = proj[i];
      if (a.w.ci !== b.w.ci || b.w.j !== a.w.j + 1) continue;
      var on = focusCi === null || a.w.ci === focusCi;
      ctx.strokeStyle = contigColor(a.w.ci, nC, (on ? 0.35 : 0.05) * Math.min(1, (a.z + b.z) / 2));
      ctx.lineWidth = focusCi === a.w.ci ? 1.4 : 0.8;
      ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
    }
    // points, far to near; the shared selection is drawn white and on top
    var sel = this.selection, selKey = sel ? G.genome.normName(sel.chrom) : null;
    var inSel = function (w) { return sel && w.contig.key === selKey && w.end >= sel.start && w.start <= sel.end; };
    proj.slice().sort(function (a, b) { return a.z - b.z; }).forEach(function (p) {
      var on = focusCi === null || p.w.ci === focusCi;
      if (inSel(p.w)) {
        ctx.fillStyle = 'white'; ctx.beginPath(); ctx.arc(p.x, p.y, 5 * Math.max(0.6, p.z), 0, 2 * Math.PI); ctx.fill();
        return;
      }
      ctx.fillStyle = contigColor(p.w.ci, nC, (on ? 0.85 : 0.12) * Math.min(1, p.z));
      ctx.beginPath(); ctx.arc(p.x, p.y, Math.max(0.8, 2.2 * p.z * Math.sqrt(self.zoom)), 0, 2 * Math.PI); ctx.fill();
    });

    // axis explanation (what each principal component is made of)
    g.setText('rgba(255,255,255,0.5)', 11, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
    m.axes.forEach(function (ax, i) {
      var txt = 'PC' + (i + 1) + ' (' + Math.round(ax.share * 100) + '%): ' + ax.load.map(function (l) { return (l.w >= 0 ? '+' : '-') + l.name; }).join(', ');
      g.fText(txt, 12, g.cH - 30 - (2 - i) * 15);
    });
    g.fText(m.windows.length.toLocaleString() + ' windows of ' + G.fmtBp(m.win) + ', coloured by chromosome, joined along each chromosome', 12, g.cH - 12);

    if (over) {
      ctx.strokeStyle = 'white'; ctx.lineWidth = 2;
      ctx.beginPath(); ctx.arc(over.x, over.y, 6, 0, 2 * Math.PI); ctx.stroke();
      var w = over.w;
      var lines = [w.contig.name + ':' + w.start.toLocaleString() + '-' + w.end.toLocaleString()];
      m.names.forEach(function (nm, fi) {
        var v = w.f[fi];
        var shown = /per kb|depth$|ends/.test(nm) ? (Math.expm1(v)).toFixed(2) : v.toFixed(2);
        lines.push(nm.replace('log ', '') + ': ' + shown);
      });
      lines.push('click to open in Arcs');
      G.app.view.drawTooltip(g, lines);
      g.setCursor('pointer');
      if (g.MOUSE_UP_FAST) {
        G.app.select({ chrom: w.contig.name, start: w.start, end: w.end });
        G.app.setMode('arcs');
        G.app.view.goTo(w.contig.name, w.start, w.end);
      }
    }
  };

  if (typeof mo !== 'undefined' && G.app) G.landscape = new Landscape();
  G.Landscape = { pca: pca, eigenSym: eigenSym };
})(globalThis.G = globalThis.G || {});
