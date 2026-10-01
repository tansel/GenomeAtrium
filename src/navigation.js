/*
 * Navigation for the Arcs view: overview strip, animated fly-to, history,
 * keyboard and trackpad, drag to zoom, and the shared selection.
 *
 *  - Overview strip (bottom): every chromosome with its cytogenetic bands,
 *    a box for what is on screen, finding ticks, and the selection. Click
 *    to fly there, drag across it to fly to that stretch.
 *  - Fly-to: zooms out, travels, zooms in (after van Wijk and Nuij's smooth
 *    zooming, simplified), so a jump between chromosomes stays readable.
 *  - History: every place the view settles is recorded; Alt+Left/Right or
 *    the toolbar arrows go back and forward.
 *  - Keys: arrows pan and zoom, + and - zoom, 0 resets. Trackpad: two-finger
 *    sideways swipe pans, pinch zooms.
 *  - Shift+drag zooms to the dragged stretch. Alt+drag selects it; the
 *    selection is shared with the Matrix and the Landscape.
 */
(function (G) {
  var View = G.View;
  var FLIGHT_MS = 900, SETTLE_MS = 700, STRIP_H = 14;
  var STAIN = {
    gneg: 'rgba(255,255,255,0.10)', gpos25: 'rgba(255,255,255,0.25)', gpos50: 'rgba(255,255,255,0.40)',
    gpos75: 'rgba(255,255,255,0.55)', gpos100: 'rgba(255,255,255,0.70)', acen: 'rgba(255,90,90,0.65)',
    gvar: 'rgba(150,150,255,0.40)', stalk: 'rgba(255,255,255,0.18)'
  };

  // ----- coordinates

  // Screen x -> { contig, pos }, clamped to the nearest contig edge.
  View.prototype.screenToGenome = function (x) {
    if (!this.segs || !this.segs.length) return null;
    var best = null, bestD = Infinity;
    for (var i = 0; i < this.segs.length; i++) {
      var sc = this.segScreen(this.segs[i]), c = this.segs[i].contig;
      if (x >= sc.x && x <= sc.x + sc.w) return { contig: c, pos: Math.max(1, Math.min(c.length, Math.round((x - sc.x) / sc.w * c.length))) };
      var d = x < sc.x ? sc.x - x : x - sc.x - sc.w;
      if (d < bestD) { bestD = d; best = { contig: c, pos: x < sc.x ? 1 : c.length }; }
    }
    return best;
  };

  // Visible stretch as { a: {contig,pos}, b: {contig,pos} }.
  View.prototype.visibleRange = function () {
    var a = this.screenToGenome(0), b = this.screenToGenome(this.g.cW);
    return a && b ? { a: a, b: b } : null;
  };

  // Visible stretch on the chromosome under the centre of the screen:
  // { chrom, start, end }. The screen edges can fall in the gap next to a
  // neighbouring chromosome; this keeps to one.
  View.prototype.visibleSpan = function () {
    var c = this.screenToGenome(this.g.cW / 2), a = this.screenToGenome(0), b = this.screenToGenome(this.g.cW);
    if (!c) return null;
    return { chrom: c.contig.name, start: a.contig === c.contig ? a.pos : 1, end: b.contig === c.contig ? b.pos : c.contig.length };
  };

  View.prototype.centerU = function () {
    return (this.g.cW / 2 - this.margins().l - this.x0) / this.scale();
  };

  // ----- fly-to

  // Replaces the instant goTo: same arguments and return value, animated.
  View.prototype.goTo = function (name, start, end, opts) {
    if (!this.data || !this.data.genome) return false;
    var c = this.data.genome.get(name);
    if (!this.segs) this.computeLayout();
    if (!c || !this.segByKey[c.key]) return false;
    start = Math.max(1, start || 1); end = Math.min(c.length, end || c.length);
    if (end <= start) { start = Math.max(1, start - 50); end = start + 100; }
    this.focus = null;
    for (var i = 0; i < this.weights.length; i++) this.weights[i] = 1;
    this.computeLayout();
    var s = this.segByKey[c.key];
    var u0 = s.x + (start - 1) / c.length * s.w, u1 = s.x + end / c.length * s.w;
    var zoom1 = Math.max(1, (this.layoutWidth / (u1 - u0)) * 0.96);
    this.startFlight(zoom1, (u0 + u1) / 2, opts);
    return true;
  };

  View.prototype.startFlight = function (zoom1, uc1, opts) {
    opts = opts || {};
    var z0 = this.zoom, uc0 = this.centerU();
    if (opts.instant || !isFinite(uc0)) { this.setCenter(zoom1, uc1); this.flight = null; this.afterMove(opts); return; }
    // How far out to go so both ends are on screen together (in log zoom).
    var dist = Math.abs(uc1 - uc0) / this.layoutWidth; // 0..1 of the genome
    var fitZoom = Math.max(1, Math.min(z0, zoom1) / Math.max(1, dist * Math.min(z0, zoom1) * 1.2));
    this.flight = { t0: Date.now(), z0: z0, z1: zoom1, u0: uc0, u1: uc1, zMid: fitZoom, opts: opts };
  };

  View.prototype.setCenter = function (zoom, uc) {
    this.zoom = this.zoomTarget = zoom;
    var k = (this.g.cW - this.margins().l - this.margins().r) / this.layoutWidth * zoom;
    this.x0 = this.g.cW / 2 - this.margins().l - uc * k;
  };

  function ease(t) { return t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2; }

  // Called from navigate(); returns true while a flight owns the view.
  View.prototype.stepFlight = function () {
    var f = this.flight;
    if (!f) return false;
    var t = Math.min(1, (Date.now() - f.t0) / FLIGHT_MS), e = ease(t);
    var lz0 = Math.log(f.z0), lz1 = Math.log(f.z1), lzm = Math.log(f.zMid);
    var lz = lz0 + (lz1 - lz0) * e;
    var dip = Math.max(0, Math.min(lz0, lz1) - lzm) * Math.sin(Math.PI * e); // zoom out on the way
    // travel faster when zoomed out: move the centre in proportion to 1/zoom
    var uc = f.u0 + (f.u1 - f.u0) * e;
    this.setCenter(Math.exp(lz - dip), uc);
    if (t >= 1) { this.setCenter(f.z1, f.u1); this.flight = null; this.afterMove(f.opts); }
    return true;
  };

  View.prototype.afterMove = function (opts) {
    this.lastMove = Date.now();
    this.historyPending = !(opts && opts.fromHistory);
  };

  View.prototype.reset = function () {
    this.focus = null;
    if (this.layoutWidth) this.startFlight(1, this.layoutWidth / 2);
    else { this.zoom = this.zoomTarget = 1; this.x0 = 0; }
  };

  // ----- history and bookmarks

  View.prototype.recordHistory = function () {
    var sp = this.visibleSpan();
    if (!sp) return;
    var e = { chrom: sp.chrom, start: sp.start, end: sp.end, whole: this.zoom <= 1.01 };
    var h = this.history = this.history || [];
    var cur = h[this.histIdx];
    if (cur && cur.chrom === e.chrom && Math.abs(cur.start - e.start) < (e.end - e.start) * 0.05 && Math.abs(cur.end - e.end) < (e.end - e.start) * 0.05) return;
    h.length = (this.histIdx == null ? -1 : this.histIdx) + 1;
    h.push(e);
    if (h.length > 100) h.shift();
    this.histIdx = h.length - 1;
    if (G.app && G.app.onHistory) G.app.onHistory();
  };

  View.prototype.historyGo = function (step) {
    var h = this.history || [], i = (this.histIdx || 0) + step;
    if (i < 0 || i >= h.length) return false;
    this.histIdx = i;
    var e = h[i];
    if (e.whole) this.startFlight(1, this.layoutWidth / 2, { fromHistory: true });
    else this.goTo(e.chrom, e.start, e.end, { fromHistory: true });
    if (G.app && G.app.onHistory) G.app.onHistory();
    return true;
  };

  // ----- input: wheel (pan and pinch), drag modes, keys

  View.prototype.installNavigation = function () {
    var self = this, cv = this.g.canvas;
    this.panVel = 0;
    cv.addEventListener('wheel', function (e) {
      if (self.mode !== 'arcs') return;
      self.flight = null;
      self.lastMove = Date.now(); self.historyPending = true;
      if (Math.abs(e.deltaX) > Math.abs(e.deltaY) && !e.ctrlKey) { self.panVel -= e.deltaX; } // two-finger sideways swipe
      else if (e.ctrlKey) { self.zoomTarget *= Math.exp(-e.deltaY * 0.012); self.anchorX = e.offsetX; } // pinch
    }, { passive: true });
    cv.addEventListener('mousedown', function (e) {
      if (self.mode !== 'arcs') return;
      self.flight = null;
      var inStrip = e.offsetY > self.g.cH - STRIP_H - 14;
      self.drag = { mode: inStrip ? 'overview' : e.shiftKey ? 'zoom' : e.altKey ? 'select' : 'pan', x0: e.offsetX, x1: e.offsetX };
    });
    window.addEventListener('mousemove', function (e) {
      if (self.drag && self.drag.mode !== 'pan') self.drag.x1 = e.clientX - cv.getBoundingClientRect().left;
    });
    window.addEventListener('mouseup', function () {
      var d = self.drag;
      self.drag = null;
      if (!d || self.mode !== 'arcs') return;
      if (d.mode === 'overview') return self.overviewRelease(d);
      if (d.mode !== 'zoom' && d.mode !== 'select') return;
      var lo = Math.min(d.x0, d.x1), hi = Math.max(d.x0, d.x1);
      if (hi - lo < 4) { if (d.mode === 'select') self.setSelection(null); return; }
      var a = self.screenToGenome(lo), b = self.screenToGenome(hi);
      if (!a || !b) return;
      if (a.contig !== b.contig) b = { contig: a.contig, pos: a.contig.length }; // one chromosome at a time
      if (d.mode === 'zoom') self.goTo(a.contig.name, a.pos, b.pos);
      else self.setSelection({ chrom: a.contig.name, start: a.pos, end: b.pos });
    });
    document.addEventListener('keydown', function (e) {
      if (self.mode !== 'arcs' || !self.data || !self.segs) return;
      var tag = (e.target && e.target.tagName) || '';
      if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA') return;
      var W = self.g.cW;
      if (e.altKey && e.key === 'ArrowLeft') { self.historyGo(-1); e.preventDefault(); return; }
      if (e.altKey && e.key === 'ArrowRight') { self.historyGo(1); e.preventDefault(); return; }
      var zoomBy = function (f) { self.flight = null; self.anchorX = W / 2; self.zoomTarget = Math.max(1, self.zoomTarget * f); self.lastMove = Date.now(); self.historyPending = true; };
      switch (e.key) {
        case 'ArrowLeft': self.flight = null; self.panVel += W * 0.25; break;
        case 'ArrowRight': self.flight = null; self.panVel -= W * 0.25; break;
        case 'ArrowUp': case '+': case '=': zoomBy(2); break;
        case 'ArrowDown': case '-': case '_': zoomBy(0.5); break;
        case '0': self.reset(); break;
        default: return;
      }
      e.preventDefault();
    });
  };

  // Pan velocity decays each frame: smooth arrow keys and swipes.
  View.prototype.applyPan = function () {
    if (Math.abs(this.panVel) < 0.5) { this.panVel = 0; return; }
    var step = this.panVel * 0.25;
    this.x0 += step; this.panVel -= step;
    this.lastMove = Date.now(); this.historyPending = true;
  };

  // ----- selection (shared with Matrix and Landscape)

  View.prototype.setSelection = function (sel) {
    this.selection = sel;
    if (G.app && G.app.onSelection) G.app.onSelection(sel);
  };

  View.prototype.drawSelection = function (g) {
    var s = this.selection, ctx = g.context;
    if (s) {
      var key = G.genome.normName(s.chrom), x0 = this.bpToX(key, s.start), x1 = this.bpToX(key, s.end);
      if (x0 !== null) {
        ctx.fillStyle = 'rgba(120,200,255,0.10)'; ctx.fillRect(x0, 0, Math.max(2, x1 - x0), g.cH - STRIP_H - 16);
        ctx.fillStyle = 'rgba(120,200,255,0.6)'; ctx.fillRect(x0, 0, 1, g.cH - STRIP_H - 16); ctx.fillRect(x1, 0, 1, g.cH - STRIP_H - 16);
      }
    }
    var d = this.drag;
    if (d && (d.mode === 'zoom' || d.mode === 'select')) {
      var lo = Math.min(d.x0, d.x1), hi = Math.max(d.x0, d.x1);
      ctx.fillStyle = d.mode === 'zoom' ? 'rgba(255,255,255,0.12)' : 'rgba(120,200,255,0.18)';
      ctx.fillRect(lo, 0, hi - lo, g.cH - STRIP_H - 16);
      g.setText('rgba(255,255,255,0.8)', 11, 'Helvetica, Arial, sans-serif', 'center', 'top');
      g.fText(d.mode === 'zoom' ? 'release to zoom here' : 'release to select (shared with Matrix and Landscape)', (lo + hi) / 2, 60);
    }
  };

  // ----- overview strip

  View.prototype.overviewLayout = function () {
    var cs = this.data.genome.contigs, W = this.g.cW, l = 40, r = 40, gap = 2;
    var total = cs.reduce(function (s, c) { return s + c.length; }, 0);
    var k = (W - l - r - gap * (cs.length - 1)) / total, x = l, out = {};
    cs.forEach(function (c) { out[c.key] = { x: x, w: c.length * k, c: c }; x += c.length * k + gap; });
    return out;
  };

  View.prototype.drawOverview = function (g) {
    var ctx = g.context, self = this, y = g.cH - STRIP_H - 6, lay = this.overviewLayout();
    this.ovLayout = lay;
    var bands = this.cytobands;
    Object.keys(lay).forEach(function (k) {
      var o = lay[k], list = bands && bands[k];
      if (list) list.forEach(function (b) {
        ctx.fillStyle = STAIN[b[3]] || STAIN.gneg;
        ctx.fillRect(o.x + (b[0] / o.c.length) * o.w, y, Math.max(0.5, (b[1] - b[0]) / o.c.length * o.w), STRIP_H);
      });
      else { ctx.fillStyle = 'rgba(255,255,255,0.15)'; ctx.fillRect(o.x, y, o.w, STRIP_H); }
      if (o.w > 14) {
        g.setText('rgba(255,255,255,0.55)', 9, 'Helvetica, Arial, sans-serif', 'center', 'bottom');
        g.fText(o.c.name.replace(/^chr/i, ''), o.x + o.w / 2, y - 1);
      }
    });
    // findings as red ticks
    ctx.fillStyle = 'rgb(255,70,70)';
    (this.findings || []).forEach(function (f) {
      var o = lay[G.genome.normName(f.chrom)];
      if (o) ctx.fillRect(o.x + f.pos / o.c.length * o.w - 0.5, y - 3, 1.5, STRIP_H + 3);
    });
    // selection
    var s = this.selection;
    if (s) {
      var os = lay[G.genome.normName(s.chrom)];
      if (os) { ctx.fillStyle = 'rgba(120,200,255,0.5)'; ctx.fillRect(os.x + s.start / os.c.length * os.w, y - 2, Math.max(2, (s.end - s.start) / os.c.length * os.w), STRIP_H + 4); }
    }
    // viewport box
    var vr = this.visibleRange();
    if (vr && this.zoom > 1.01) {
      var oa = lay[vr.a.contig.key], ob = lay[vr.b.contig.key];
      var xa = oa.x + vr.a.pos / oa.c.length * oa.w, xb = ob.x + vr.b.pos / ob.c.length * ob.w;
      ctx.strokeStyle = 'white'; ctx.lineWidth = 1.5;
      ctx.strokeRect(Math.min(xa, xb - 3) - 1, y - 3, Math.max(3, xb - xa) + 2, STRIP_H + 6);
    }
    // drag preview and hover
    var d = this.drag;
    if (d && d.mode === 'overview' && Math.abs(d.x1 - d.x0) > 3) {
      ctx.fillStyle = 'rgba(255,255,255,0.25)'; ctx.fillRect(Math.min(d.x0, d.x1), y - 3, Math.abs(d.x1 - d.x0), STRIP_H + 6);
    }
    if (g.mY > y - 12 && g.mY < y + STRIP_H + 4) {
      var hit = this.overviewAt(g.mX);
      if (hit) {
        g.setCursor('pointer');
        var band = bands && bands[hit.c.key] ? bands[hit.c.key].find(function (b) { return hit.pos > b[0] && hit.pos <= b[1]; }) : null;
        this.hover = { lines: [hit.c.name + (band ? band[2] : '') + '  ' + hit.pos.toLocaleString(), 'click to fly here, drag to fly to a stretch'] };
      }
    }
  };

  View.prototype.overviewAt = function (x) {
    var lay = this.ovLayout;
    if (!lay) return null;
    for (var k in lay) {
      var o = lay[k];
      if (x >= o.x - 1 && x <= o.x + o.w + 1) return { c: o.c, pos: Math.max(1, Math.min(o.c.length, Math.round((x - o.x) / o.w * o.c.length))) };
    }
    return null;
  };

  View.prototype.overviewRelease = function (d) {
    var a = this.overviewAt(Math.min(d.x0, d.x1)), b = this.overviewAt(Math.max(d.x0, d.x1));
    if (!a) return;
    if (Math.abs(d.x1 - d.x0) > 3 && b) {
      if (b.c !== a.c) b = { c: a.c, pos: a.c.length };
      this.goTo(a.c.name, a.pos, b.pos);
    } else { // click: keep the current span, centre on the point
      var vr = this.visibleRange(), span = vr && vr.a.contig === vr.b.contig ? vr.b.pos - vr.a.pos : a.c.length / 4;
      span = Math.min(span, a.c.length);
      this.goTo(a.c.name, Math.round(a.pos - span / 2), Math.round(a.pos + span / 2));
    }
  };

  // Cytobands: UCSC cytoBandIdeo text (chrom, start, end, name, stain).
  View.prototype.setCytobands = function (text) {
    var out = {};
    text.split('\n').forEach(function (line) {
      var f = line.split('\t');
      if (f.length < 5) return;
      var k = G.genome.normName(f[0]);
      (out[k] = out[k] || []).push([+f[1], +f[2], f[3], f[4]]);
    });
    this.cytobands = out;
  };

  G.navigation = { STRIP_H: STRIP_H };
})(globalThis.G = globalThis.G || {});
