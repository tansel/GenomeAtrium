/*
 * Genome arcs view, drawn with Moebio Framework's mo.Graphics.
 *
 * Modelled on moebio.com/attention: the genome is laid out along one line
 * the way that piece lays out words, and relations between two places are
 * arcs above the line. Chromosomes can be opened with a fisheye (click a
 * label), and the wheel zooms smoothly around the cursor down to single
 * bases. Below the line sits a density band (VCF) or depth band (BAM).
 */
(function (G) {
  var BG = 'rgb(20,20,20)';
  var COLORS = {
    snv: 'rgb(120,180,255)', indel: 'rgb(255,170,80)', sv: 'rgb(235,90,200)',
    depth: 'rgb(90,210,190)', callable: 'rgba(255,255,255,0.28)', lowdp: 'rgb(230,170,40)', nocall: 'rgba(200,60,60,0.55)',
    arc: {
      5: 'rgb(255,95,95)', 6: 'rgb(90,150,255)', 7: 'rgb(255,200,60)', 8: 'rgb(120,220,120)',
      9: 'rgb(180,140,255)', 10: 'rgb(235,90,200)', 3: 'rgb(255,95,95)', 2: 'rgb(120,220,120)',
      pair_inter: 'rgb(235,90,200)', pair_long: 'rgb(255,95,95)', junction: 'rgb(120,220,120)'
    }
  };
  var TICK_COLORS = ['rgb(120,180,255)', 'rgb(120,180,255)', 'rgb(120,220,120)', 'rgb(255,120,100)', 'rgb(255,170,80)'];
  var MAX_ARCS = 5000;
  var TICK_MAX_BP_PER_PX = 5000, TICK_MAX_VISIBLE = 15000;
  // Tick colours built once: [type group][pass] -> rgba string. Building them
  // per tick per frame made thousands of throwaway strings (GC pauses).
  var TICK_FILL = TICK_COLORS.concat([COLORS.sv]).map(function (c) { return [withAlpha(c, 0.3), withAlpha(c, 0.95)]; });
  var MIN_BP_PER_PX = 1 / 24; // deepest zoom: 24 px per base

  function withAlpha(rgb, a) { return rgb.replace('rgb(', 'rgba(').replace(')', ',' + a.toFixed(3) + ')'); }
  function fmtBp(n) {
    if (n >= 1e6) return (n / 1e6).toFixed(n >= 1e7 ? 1 : 2).replace(/\.?0+$/, '') + ' Mb';
    if (n >= 1e3) return (n / 1e3).toFixed(n >= 1e4 ? 0 : 1).replace(/\.0$/, '') + ' kb';
    return Math.round(n) + ' bp';
  }
  function lowerBound(arr, n, v) {
    var lo = 0, hi = n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  function View(containerSelector) {
    var self = this;
    this.data = null;
    this.layers = { arcs: true, similar: true, snv: true, indel: true, het: true, depth: true, ticks: true, callable: true, clinvar: true,
      genes: true, regulatory: true, gwas: true, antisense: true, lncRNA: true, smallRNA: true, pseudogene: false,
      findings: true, landscape: true, labels: true, roh: true, panelRing: true, methyl: true, prs: true, atriumSimilar: false, snapTurn: false, vignette: true };
    // switches survive reloads (per browser)
    try { var saved = JSON.parse(localStorage.getItem('genomeatrium.layers') || '{}'); for (var k in saved) if (k in this.layers) this.layers[k] = !!saved[k]; } catch (e) { /* not kept */ }
    this.focusGene = null;
    this.findings = [];
    this.zoom = 1; this.zoomTarget = 1; this.x0 = 0; this.anchorX = 0;
    this.focus = null; this.weights = [];
    this.gain = {}; this.status = null;

    this.g = new mo.Graphics({
      container: containerSelector,
      init: function () {},
      cycle: function () { self.cycle(this); }
    });
    this.g.setBackgroundColor(BG);
    this._hiDPI();

    // Standard wheel event (Moebio listens to the old "mousewheel" only).
    this.g.canvas.addEventListener('wheel', function (e) {
      e.preventDefault();
      if (self.mode !== 'arcs') return;
      if (e.ctrlKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) return; // pinch and sideways swipes: navigation.js
      self.flight = null; self.lastMove = Date.now(); self.historyPending = true;
      var d = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY;
      self.zoomTarget *= Math.exp(-d * 0.0018);
      self.anchorX = e.offsetX;
    }, { passive: false });
    this.g.canvas.addEventListener('dblclick', function () { if (self.mode === 'arcs') self.reset(); });
    if (this.installNavigation) this.installNavigation();
  }

  // Draw at device resolution so text and thin arcs stay sharp on Retina.
  View.prototype._hiDPI = function () {
    var g = this.g, orig = g._adjustCanvas.bind(g), self = this;
    g._adjustCanvas = function (dim) {
      orig(dim);
      var dpr = self.forceDpr || window.devicePixelRatio || 1; // 1 while drawn into the VR panel
      g.canvas.width = g.cW * dpr; g.canvas.height = g.cH * dpr;
      g.canvas.style.width = g.cW + 'px'; g.canvas.style.height = g.cH + 'px';
      g.context.setTransform(dpr, 0, 0, dpr, 0, 0);
    };
    g._adjustCanvas();
  };

  View.prototype.setData = function (data) {
    this.data = data;
    this.gain = {};
    this.focus = null;
    this.hover = null;
    if (data && data.genome) {
      this.weights = data.genome.contigs.map(function () { return 1; });
      this.visibleContigs = {};
      var self = this;
      data.genome.contigs.forEach(function (c) { self.visibleContigs[c.key] = c; });
      // Arcs whose both ends are laid out, most supported first, capped.
      this.arcs = (data.arcs || []).filter(function (a) {
        return self.visibleContigs[G.genome.normName(a.c0)] && self.visibleContigs[G.genome.normName(a.c1)];
      }).slice(0, MAX_ARCS);
      this.maxSupport = this.arcs.reduce(function (m, a) { return Math.max(m, a.support); }, 1);
      // Similarity arcs from the shared window model (src/windows.js).
      var mdl = G.windows ? G.windows.model(data) : null;
      this.windowModel = mdl;
      this.simArcs = mdl ? G.windows.pairsToArcs(mdl) : [];
      this.winIndex = {};
      if (mdl) mdl.windows.forEach(function (w, i) { self.winIndex[w.contig.key + ':' + w.j] = i; });
    }
    this.reset();
  };

  View.prototype.setStatus = function (s) { this.status = s; };
  View.prototype.setFindings = function (list) { this.findings = list || []; };

  View.prototype.reset = function () {
    this.zoom = this.zoomTarget = 1;
    this.x0 = 0;
    this.focus = null;
  };

  // ----- layout: contig -> [x, w] in layout units, then to screen pixels

  View.prototype.margins = function () { return { l: 40, r: 40 }; };

  View.prototype.computeLayout = function () {
    var cs = this.data.genome.contigs, total = 0, i;
    for (i = 0; i < cs.length; i++) total += cs[i].length;
    // Fisheye target: the focused contig takes about 55% of the line.
    var targets = cs.map(function () { return 1; });
    if (this.focus !== null) {
      var lf = cs[this.focus].length, rest = total - lf;
      targets[this.focus] = lf / total >= 0.55 ? 1 : 0.55 * rest / (0.45 * lf);
    }
    var sum = 0;
    for (i = 0; i < cs.length; i++) {
      this.weights[i] = this.weights[i] * 0.88 + targets[i] * 0.12;
      sum += cs[i].length * this.weights[i];
    }
    var gap = cs.length > 1 ? sum * 0.004 : 0;
    var x = 0, segs = [];
    for (i = 0; i < cs.length; i++) {
      var w = cs[i].length * this.weights[i];
      segs.push({ contig: cs[i], x: x, w: w, i: i });
      x += w + gap;
    }
    this.layoutWidth = x - gap;
    this.segs = segs;
    this.segByKey = {};
    for (i = 0; i < segs.length; i++) this.segByKey[segs[i].contig.key] = segs[i];
  };

  View.prototype.scale = function () {
    var m = this.margins();
    return (this.g.cW - m.l - m.r) / this.layoutWidth * this.zoom;
  };
  View.prototype.segScreen = function (s) {
    var k = this.scale(), l = this.margins().l;
    return { x: l + this.x0 + s.x * k, w: s.w * k };
  };
  View.prototype.bpToX = function (key, pos) {
    var s = this.segByKey[key];
    if (!s) return null;
    var sc = this.segScreen(s);
    return sc.x + (pos - 0.5) / s.contig.length * sc.w;
  };

  // Moves the view so contig:start-end fills the screen.
  View.prototype.goTo = function (name, start, end) {
    if (!this.data || !this.data.genome) return false;
    var c = this.data.genome.get(name);
    if (!c || !this.segByKey[c.key]) return false;
    start = Math.max(1, start || 1); end = Math.min(c.length, end || c.length);
    if (end <= start) { start = Math.max(1, start - 50); end = start + 100; }
    this.focus = null;
    for (var i = 0; i < this.weights.length; i++) this.weights[i] = 1;
    this.computeLayout();
    var s = this.segByKey[c.key], m = this.margins(), avail = this.g.cW - m.l - m.r;
    var u0 = s.x + (start - 1) / c.length * s.w, u1 = s.x + end / c.length * s.w;
    this.zoom = this.zoomTarget = Math.max(1, (this.layoutWidth / (u1 - u0)) * 0.96);
    var k = avail / this.layoutWidth * this.zoom;
    this.x0 = avail * 0.02 - u0 * k;
    return true;
  };

  // ----- frame

  // One bad frame must not take the page down: report it and keep drawing.
  View.prototype.cycle = function (g) {
    try { this.frame(g); this.lastError = null; }
    catch (err) {
      if (!this.lastError || this.lastError.message !== err.message) console.error('frame failed', err);
      this.lastError = err;
      g.setText('rgb(255,120,100)', 12, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
      g.fText('Drawing error: ' + err.message + ' (details in the browser console)', 12, g.cH - 40);
    }
  };

  // The view the canvas draws: the page's mode, or inside the Atrium the window's tab.
  View.prototype.activeMode = function () {
    return this.mode === 'atrium' && this.panelMode ? this.panelMode : this.mode;
  };

  View.prototype.frame = function (g) {
    var d = this.data, mode = this.activeMode();
    if (this.status) this.drawStatus(g);
    if (!d) { if (!this.status) this.drawEmpty(g); return; }
    if (d.format === 'bam' && !d.stats.aligned) { this.drawUnaligned(g); return; }
    if (!d.genome.contigs.length) { this.drawMessage(g, 'No contigs with data in this file.'); return; }
    if (mode === '3d' && G.landscape) { G.landscape.draw(g); return; }
    if (mode === 'matrix' && G.matrix) { G.matrix.draw(g); return; }
    if (mode === 'hilbert' && G.app.hilbert) { G.app.hilbert.draw(g); return; }
    if (mode === 'circos' && G.app.circos) { G.app.circos.draw(g); return; }
    if (mode === 'gene' && G.app.geneView) { G.app.geneView.draw(g); return; }
    if (mode === 'protein' && G.app.proteinView) { G.app.proteinView.draw(g); return; }
    if (mode === 'hic' && G.app.hicView) { G.app.hicView.draw(g); return; }
    if (mode === 'pathways' && G.app.pathways) { G.app.pathways.draw(g); return; }

    this.computeLayout();
    this.navigate(g);
    var axisY = Math.round(g.cH * (this.bandRows().length > 1 ? 0.5 : 0.6));
    this.axisY = axisY;
    this.hover = null;

    if (this.drawSelection) this.drawSelection(g);
    this.drawBand(g, axisY);
    this.drawAxis(g, axisY);
    this.drawArcs(g, axisY);
    if (this.regOn()) this.drawRegulatory(g, axisY); else this.regVisible = 0;
    if (this.genesOn()) this.drawGenes(g, axisY);
    if (this.genesOn() && this.layers.antisense) this.drawAntisense(g, axisY);
    if (this.gwas && this.layers.gwas && this.grch38()) this.drawGwas(g, axisY);
    if (this.layers.ticks && d.format === 'vcf') this.drawTicks(g, axisY);
    if (this.findings.length && this.layers.findings) this.drawFindings(g, axisY);
    if (this.drawOverview) this.drawOverview(g);
    this.drawReadout(g, axisY);
    if (this.hover) this.drawTooltip(g, this.hover.lines);
  };

  // Wheel zoom around the cursor, drag to pan, keep at least some genome on screen.
  View.prototype.navigate = function (g) {
    var m = this.margins(), avail = g.cW - m.l - m.r;
    if (this.stepFlight && this.stepFlight()) return;
    if (this.historyPending && Date.now() - (this.lastMove || 0) > 700 && !g.MOUSE_PRESSED) { this.historyPending = false; this.recordHistory(); }
    var maxZoom = this.layoutWidth / (avail * MIN_BP_PER_PX) * (this.layoutWidth / this.data.genome.totalLength());
    this.zoomTarget = Math.max(1, Math.min(this.zoomTarget, Math.max(1, maxZoom)));
    var before = this.scale();
    var u = (this.anchorX - m.l - this.x0) / before;
    this.zoom += (this.zoomTarget - this.zoom) * 0.22;
    var after = this.scale();
    this.x0 = this.anchorX - m.l - u * after;
    if (g.MOUSE_PRESSED && (!this.drag || this.drag.mode === 'pan')) {
      if (g.DX_MOUSE) { this.x0 += g.DX_MOUSE; this.lastMove = Date.now(); this.historyPending = true; }
      this.zoomTarget = this.zoom;
    }
    if (this.applyPan) this.applyPan();
    var span = this.layoutWidth * after;
    this.x0 = Math.min(avail * 0.5, Math.max(this.x0, avail * 0.5 - span));
    if (this.zoom <= 1.0005 && !g.MOUSE_PRESSED) this.x0 *= 0.85; // settle when fully zoomed out
  };

  View.prototype.drawAxis = function (g, y) {
    var ctx = g.context, self = this;
    var clickedLabel = null;
    this.segs.forEach(function (s, i) {
      var sc = self.segScreen(s);
      if (sc.x > g.cW || sc.x + sc.w < 0) return;
      var over = g.mY > y - 8 && g.mY < y + 30 && g.mX >= sc.x && g.mX <= sc.x + sc.w;
      ctx.fillStyle = i === self.focus ? 'rgba(255,255,255,0.95)' : over ? 'rgba(255,255,255,0.8)' : i % 2 ? 'rgba(255,255,255,0.35)' : 'rgba(255,255,255,0.55)';
      ctx.fillRect(Math.max(sc.x, -2), y - 1, Math.min(sc.w, g.cW + 4), 2);

      var name = s.contig.name;
      var fs = Math.max(9, Math.min(13, sc.w * 0.3));
      g.setText(over || i === self.focus ? 'white' : 'rgba(255,255,255,0.6)', fs, 'Helvetica, Arial, sans-serif', 'center', 'top');
      var tw = g.getTextW(name);
      if (tw < sc.w - 2 || over) {
        var lx = Math.max(sc.x + tw / 2 + 2, Math.min(sc.x + sc.w - tw / 2 - 2, g.cW / 2));
        if (sc.w < g.cW) lx = sc.x + sc.w / 2;
        g.fText(name, lx, y + 8);
      }
      if (over) {
        self.overContig = i;
        g.setCursor('pointer');
        if (g.MOUSE_UP_FAST) clickedLabel = i;
        if (!self.hover && g.mY > y + 6) self.hover = { lines: [name + '  ' + fmtBp(s.contig.length), 'click to open, double click to reset'] };
      }
      self.drawRuler(g, s, sc, y);
    });
    if (clickedLabel !== null) this.focus = this.focus === clickedLabel ? null : clickedLabel;
  };

  // Position ticks along a contig once it is wide enough to read them.
  View.prototype.drawRuler = function (g, s, sc, y) {
    if (sc.w < 260) return;
    var bpPerPx = s.contig.length / sc.w;
    var raw = bpPerPx * 110, p10 = Math.pow(10, Math.floor(Math.log10(raw)));
    var step = [1, 2, 5, 10].map(function (f) { return f * p10; }).find(function (v) { return v >= raw; });
    var first = Math.max(step, Math.ceil(((0 - sc.x) * bpPerPx) / step) * step);
    g.setText('rgba(255,255,255,0.4)', 9, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.context.fillStyle = 'rgba(255,255,255,0.3)';
    for (var p = first; p < s.contig.length; p += step) {
      var x = sc.x + (p - 0.5) / s.contig.length * sc.w;
      if (x > g.cW) break;
      if (x < 0) continue;
      g.context.fillRect(x, y + 24, 1, 4);
      g.fText(step >= 1e6 ? fmtBp(p) : p.toLocaleString(), x, y + 29); // full numbers once steps are finer than 1 Mb
    }
  };

  // Rows below the axis, one per measure, drawn as spaced ticks growing up
  // from a baseline. VCF: SNV, indel and het fraction. BAM: depth.
  // Right under the axis: the gVCF callable strip and the ClinVar strip.
  View.prototype.bandRows = function () {
    var d = this.data, L = this.layers;
    var rows = d.format === 'vcf'
      ? [{ key: 'snv', label: 'SNV/MNV', color: COLORS.snv, series: 'snv' },
         { key: 'indel', label: 'indel', color: COLORS.indel, series: 'indel' },
         { key: 'het', label: 'het fraction', color: 'rgb(190,140,255)', ratio: true }]
      : [{ key: 'depth', label: 'mean depth', color: COLORS.depth, series: 'depth', mean: true }];
    if (d.methyl) rows.push({ key: 'methyl', label: 'methylation', methyl: true });
    if (d.prs) rows.push({ key: 'prs', label: 'score ' + d.prs.info.id, prs: true });
    return rows.filter(function (r) { return L[r.key] !== false; });
  };

  View.prototype.drawBand = function (g, axisY) {
    var d = this.data, ctx = g.context, self = this, isVcf = d.format === 'vcf';
    var rows = this.bandRows(), STEP = 2;
    var top = axisY + (this.genesOn() ? 66 : 46), bandH = Math.max(40, g.cH - top - 62), gapY = 14;
    var rowH = rows.length ? (bandH - gapY * (rows.length - 1)) / rows.length : 0;
    var columns = [];

    this.segs.forEach(function (s) {
      var tr = d.tracks[s.contig.key];
      if (!tr) return;
      var sc = self.segScreen(s);
      var x0 = Math.max(0, Math.floor(sc.x)), x1 = Math.min(g.cW, Math.ceil(sc.x + sc.w));
      if (x1 <= x0) return;
      var bpPerPx = s.contig.length / sc.w, span = bpPerPx * STEP;
      var lv = {};
      rows.forEach(function (r) {
        if (r.ratio) { lv.het = tr.het.levelFor(span); lv.hom = tr.hom.levelFor(span); }
        else if (r.methyl) { var mt = d.methyl.tracks(s.contig.key); if (mt) { lv.mMeth = mt.meth.levelFor(span); lv.mCov = mt.cov.levelFor(span); } }
        else if (r.prs) { var pt = d.prs.tracks[s.contig.key]; if (pt) { lv.pPos = pt.pos.levelFor(span); lv.pNeg = pt.neg.levelFor(span); } }
        else lv[r.key] = tr[r.series].levelFor(span);
      });
      var cal = isVcf && d.isGvcf && self.layers.callable ? [tr.callable.levelFor(bpPerPx), tr.lowdp.levelFor(bpPerPx)] : null;
      var cvTrack = isVcf && self.layers.clinvar && self.clinvar && d.clinvarStatus === 'matched' ? self.clinvar.trackFor(s.contig.key, s.contig.length, d.binSize) : null;
      var cvLevel = cvTrack ? cvTrack.levelFor(bpPerPx) : null;
      var gwTrack = isVcf && self.layers.gwas && self.gwas && self.grch38() ? self.gwasTrack(s.contig.key, s.contig.length, d.binSize) : null;
      var gwLevel = gwTrack ? gwTrack.levelFor(bpPerPx) : null;
      for (var px = x0 - ((x0 - Math.floor(sc.x)) % STEP); px < x1; px += STEP) {
        var a = (px - sc.x) * bpPerPx, b = a + span;
        if (b < 0 || a > s.contig.length) continue;
        var v = {};
        rows.forEach(function (r) {
          if (r.ratio) {
            var he = sampleLevel(lv.het, a, b, 'sum'), ho = sampleLevel(lv.hom, a, b, 'sum');
            v.het = he + ho >= 3 ? he / (he + ho) : null; v.nGt = he + ho;
          } else if (r.prs) {
            v.prs = lv.pPos ? sampleLevel(lv.pPos, a, b, 'sum') - sampleLevel(lv.pNeg, a, b, 'sum') : null;
            if (v.prs === 0) v.prs = null;
          } else if (r.methyl) {
            var mc = lv.mCov ? sampleLevel(lv.mCov, a, b, 'sum') : 0;
            v.methyl = mc >= G.methylation.MIN_COV ? sampleLevel(lv.mMeth, a, b, 'sum') / mc : null; v.mCov = mc;
          } else v[r.key] = sampleLevel(lv[r.key], a, b, r.mean ? 'mean' : 'sum');
        });
        var calv = cal ? [sampleLevel(cal[0], a, b, 'mean'), sampleLevel(cal[1], a, b, 'mean')] : null;
        columns.push({ px: px, v: v, cal: calv, cv: cvLevel ? sampleLevel(cvLevel, a, b, 'sum') : null, gw: gwLevel ? sampleLevel(gwLevel, a, b, 'sum') : null, a: a, b: b, contig: s.contig });
      }
    });

    // Per-row gain: the 98th percentile of the visible columns, eased.
    rows.forEach(function (r) {
      if (r.ratio || r.methyl) { r.gain = 1; return; }
      var vals = columns.map(function (c) { return r.prs ? Math.abs(c.v[r.key] || 0) : c.v[r.key]; }).sort(function (x, y) { return x - y; });
      var target = vals.length ? vals[Math.floor(vals.length * 0.98)] || vals[vals.length - 1] : 0;
      self.gain[r.key] = !self.gain[r.key] ? target : self.gain[r.key] * 0.85 + target * 0.15;
      r.gain = self.gain[r.key] || 1;
    });
    var cvs = columns.map(function (c) { return c.cv || 0; }).sort(function (x, y) { return x - y; });
    var cvT = cvs.length ? cvs[Math.floor(cvs.length * 0.98)] : 0;
    this.gain.cv = !this.gain.cv ? cvT : this.gain.cv * 0.85 + cvT * 0.15;
    var cvGain = this.gain.cv || 1;

    rows.forEach(function (r, ri) {
      var base = top + ri * (rowH + gapY) + rowH;
      ctx.fillStyle = 'rgba(255,255,255,0.12)';
      ctx.fillRect(0, base, g.cW, 1);
      if (r.ratio || r.methyl) { // reference line at 0.5
        ctx.fillStyle = 'rgba(255,255,255,0.08)';
        ctx.fillRect(0, base - rowH * 0.5, g.cW, 1);
      }
      ctx.fillStyle = r.color;
      columns.forEach(function (c) {
        var val = c.v[r.key];
        if (val == null) return;
        var h = Math.min(1, val / r.gain) * rowH;
        if (r.methyl) { ctx.fillStyle = G.methylation.Methylation.color(val, 0.85); h = Math.max(1, h); } // height and colour: the methylated share
        if (r.prs) { ctx.fillStyle = val > 0 ? 'rgba(255,110,90,0.9)' : 'rgba(110,170,255,0.9)'; h = Math.max(1, Math.min(1, Math.abs(val) / r.gain) * rowH); }
        if (h >= 0.5) ctx.fillRect(c.px, base - h, 1, h);
      });
      g.setText('rgba(255,255,255,0.45)', 10, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
      g.fText(r.label + (r.methyl ? '  (' + (G.methylation.Methylation.CODE_NAMES[d.methyl.code] || d.methyl.code) + ' share, 0 to 1, blue to red; line at 0.5)' : r.ratio ? '  (0 to 1, line at 0.5; low = homozygous run or hemizygous)' :
        '  (full height = ' + (r.mean ? r.gain.toFixed(1) + 'x' : r.gain.toFixed(r.gain < 10 ? 1 : 0) + ' per ' + STEP + ' px') + ')'), 12, base - rowH - 2);
    });

    columns.forEach(function (c) {
      if (c.cal) { // three states: confidently called, low depth, not called
        var called = Math.min(1, c.cal[0]), low = Math.min(1, c.cal[1]);
        ctx.fillStyle = low > 0.05 ? withAlpha('rgb(230,170,40)', 0.3 + 0.7 * low) : called > 0.02 ? withAlpha('rgb(255,255,255)', 0.1 + 0.35 * called) : COLORS.nocall;
        ctx.fillRect(c.px, axisY + 3, STEP, 3);
      }
      if (c.gw) { // GWAS Catalog SNP density
        ctx.fillStyle = 'rgba(190,140,255,' + Math.min(0.9, 0.15 + 0.2 * Math.sqrt(c.gw)).toFixed(2) + ')';
        ctx.fillRect(c.px, axisY + 10, STEP, 2);
      }
      if (c.cv) { // ClinVar P/LP site density
        ctx.fillStyle = withAlpha('rgb(255,80,80)', Math.min(0.9, 0.15 + 0.75 * Math.sqrt(c.cv / cvGain)));
        ctx.fillRect(c.px, axisY + 7, STEP, 2);
      }
    });

    // Hover: every row's value for the column under the cursor.
    if (g.mY > top - 12 && g.mY < top + bandH) {
      var mx = Math.floor(g.mX);
      var col = columns.find(function (c) { return mx >= c.px && mx < c.px + STEP; });
      if (col) {
        ctx.fillStyle = 'rgba(255,255,255,0.4)'; ctx.fillRect(col.px, top, 1, bandH);
        var lines = [col.contig.name + ':' + Math.max(1, Math.round(col.a + 1)).toLocaleString() + '-' + Math.round(col.b).toLocaleString()];
        rows.forEach(function (r) {
          var val = col.v[r.key];
          if (r.prs) { lines.push('score contribution: ' + (val == null ? 'none here' : (val > 0 ? '+' : '') + val.toFixed(4))); return; }
          if (r.methyl) { lines.push('methylation: ' + (val == null ? 'too few calls' : Math.round(100 * val) + '% (' + Math.round(col.v.mCov) + ' calls)')); return; }
          lines.push(r.label + ': ' + (val == null ? 'too few genotypes' : r.ratio ? val.toFixed(2) + ' (' + Math.round(col.v.nGt) + ' genotypes)' :
            r.mean ? val.toFixed(2) + 'x' : val.toFixed(val < 10 ? 1 : 0)));
        });
        if (col.cal) lines.push('callable ' + Math.round(Math.min(1, col.cal[0]) * 100) + '%, low depth ' + Math.round(Math.min(1, col.cal[1]) * 100) + '%');
        if (col.cv != null) lines.push(col.cv.toFixed(col.cv < 10 ? 1 : 0) + ' ClinVar P/LP sites');
        this.hover = { lines: lines };
      }
    }
  };

  function sampleLevel(lv, a, b, kind) {
    var bs = lv.binSize, arr = lv.data;
    var i0 = Math.max(0, Math.floor(a / bs)), i1 = Math.min(arr.length - 1, Math.floor((b - 1e-6) / bs));
    if (i1 < i0) return 0;
    if (kind === 'mean') {
      var s = 0; for (var i = i0; i <= i1; i++) s += arr[i];
      return s / (i1 - i0 + 1);
    }
    // sum of counts in [a, b): partial bins contribute their overlap share
    var t = 0;
    for (var j = i0; j <= i1; j++) {
      var lo = Math.max(a, j * bs), hi = Math.min(b, (j + 1) * bs);
      if (hi > lo) t += arr[j] * (hi - lo) / bs;
    }
    return t;
  }

  function contigColor(i, n, a) {
    return 'hsla(' + ((i / Math.max(1, n)) * 300).toFixed(0) + ',70%,62%,' + a.toFixed(3) + ')';
  }
  G.contigColor = contigColor;

  // The window under the cursor when it is on the axis, like pressing a word
  // in moebio.com/attention: only that window's similarity arcs stay lit.
  View.prototype.windowAtCursor = function (g, axisY) {
    var mdl = this.windowModel;
    if (!mdl || Math.abs(g.mY - axisY) > 10) return null;
    for (var i = 0; i < this.segs.length; i++) {
      var sc = this.segScreen(this.segs[i]);
      if (g.mX < sc.x || g.mX > sc.x + sc.w) continue;
      var c = this.segs[i].contig, bp = (g.mX - sc.x) / sc.w * c.length;
      var wi = this.winIndex[c.key + ':' + Math.floor(bp / mdl.win)];
      return wi === undefined ? null : wi;
    }
    return null;
  };

  // Arcs: ellipses above the axis, flattened to fit the screen when tall.
  // Similarity arcs are drawn first, as a background; SV and read-pair arcs on top.
  View.prototype.drawArcs = function (g, axisY) {
    var ctx = g.context, self = this, nn = G.genome.normName;
    var maxRy = axisY - 40, best = null, bestErr = 4;
    var overKey = this.overContig != null && this.segs[this.overContig] ? this.segs[this.overContig].contig.key : null;
    this.overContig = null;
    var focusWin = this.layers.similar ? this.windowAtCursor(g, axisY) : null;
    var nC = this.segs.length, drawn = 0;

    if (focusWin !== null) { // mark the pressed window on the axis
      var fw = this.windowModel.windows[focusWin];
      var fx0 = this.bpToX(fw.contig.key, fw.start), fx1 = this.bpToX(fw.contig.key, fw.end);
      ctx.fillStyle = 'rgba(255,255,255,0.9)'; ctx.fillRect(fx0, axisY - 3, Math.max(2, fx1 - fx0), 6);
    }

    function one(a, sim) {
      var k0 = nn(a.c0), k1 = nn(a.c1);
      var xa = self.bpToX(k0, a.p0), xb = self.bpToX(k1, a.p1);
      if (xa === null || xb === null) return;
      var lx = Math.min(xa, xb), rx = Math.max(xa, xb);
      if (rx < 0 || lx > g.cW) return;
      var r = Math.max(sim ? 0 : 1.5, (rx - lx) / 2); // keep short SVs visible as small caps
      if (r < 0.6 || r > 60000) return;
      var ry = Math.min(r, maxRy * (0.35 + 0.65 * Math.min(1, r / (g.cW * 0.5))));
      var cx = (lx + rx) / 2, alpha, width, color;
      if (sim) {
        var t = Math.max(0, Math.min(1, (a.support - 0.9) / 0.1));
        var touches = focusWin !== null && (a.wi === focusWin || a.wj === focusWin);
        if (focusWin !== null && !touches) return;
        alpha = touches ? 0.85 : 0.04 + 0.4 * t * t * t;
        if (self.regVisible && !touches) alpha *= 0.2;
        width = touches ? 1.4 : 0.4 + 1.6 * Math.pow(t, 4);
        if (overKey && k0 !== overKey && k1 !== overKey) alpha *= 0.1;
        var far = self.segByKey[k1];
        color = contigColor(far ? far.i : 0, nC, alpha);
      } else {
        var strength = Math.sqrt(a.support / self.maxSupport);
        alpha = (0.35 + 0.6 * strength) * (a.pass ? 1 : 0.35);
        if (overKey && k0 !== overKey && k1 !== overKey) alpha *= 0.15;
        width = 1 + Math.log2(1 + a.support) * 0.7;
        color = withAlpha(COLORS.arc[a.type] || 'rgb(200,200,200)', alpha);
      }
      ctx.strokeStyle = color; ctx.lineWidth = width;
      ctx.beginPath(); ctx.ellipse(cx, axisY - 1, r, ry, 0, Math.PI, 2 * Math.PI); ctx.stroke();
      drawn++;
      if (g.mY < axisY - 2) { // hover test on the outline
        var nx = (g.mX - cx) / r, ny = (g.mY - axisY) / ry;
        var err = Math.abs(Math.sqrt(nx * nx + ny * ny) - 1) * Math.min(r, ry);
        if (err < bestErr || (err < 4 && !sim && best && best[4])) { bestErr = err; best = [a, cx, r, ry, sim]; }
      }
    }

    if (this.layers.similar) for (var i = 0; i < this.simArcs.length; i++) one(this.simArcs[i], true);
    if (this.layers.arcs) for (var j = 0; j < this.arcs.length; j++) one(this.arcs[j], false);

    if (best) {
      ctx.strokeStyle = 'white'; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.ellipse(best[1], axisY - 1, best[2], best[3], 0, Math.PI, 2 * Math.PI); ctx.stroke();
      var a0 = best[0];
      this.hover = { lines: best[4] ? [a0.label].concat(a0.detail) : [a0.label, 'support ' + a0.support + (a0.pass ? '' : ', filtered')] };
    }
    if (focusWin !== null && !this.hover) {
      var w = this.windowModel.windows[focusWin];
      var n = this.simArcs.filter(function (a) { return a.wi === focusWin || a.wj === focusWin; }).length;
      this.hover = { lines: [w.contig.name + ':' + w.start.toLocaleString() + '-' + w.end.toLocaleString(),
        G.windows.describe(this.windowModel, w), n + ' similar window' + (n === 1 ? '' : 's') + ' linked'] };
    }
    this.arcsDrawn = drawn;
  };

  // Single variants as ticks once zoomed in far enough.
  View.prototype.drawTicks = function (g, axisY) {
    var d = this.data, ctx = g.context, self = this, Z = G.vcf.Z;
    var best = null, bestDist = 5;
    this.segs.forEach(function (s) {
      var cols = d.variants[s.contig.key];
      if (!cols || !cols.n) return;
      var sc = self.segScreen(s);
      if (sc.x > g.cW || sc.x + sc.w < 0) return;
      var bpPerPx = s.contig.length / sc.w;
      if (bpPerPx > TICK_MAX_BP_PER_PX) return; // the rows below cover wider views
      var start = Math.max(1, Math.floor((0 - sc.x) * bpPerPx)), end = Math.ceil((g.cW - sc.x) * bpPerPx);
      var i0 = lowerBound(cols.pos, cols.n, start), i1 = lowerBound(cols.pos, cols.n, end + 1);
      if (i1 - i0 > TICK_MAX_VISIBLE) return;
      var showText = bpPerPx < 0.2;
      var gn = self.gnomadOn() && end - start <= GNOMAD_MAX_VIEW ? self.gnomad : null;
      if (gn && Math.abs(self.zoomTarget - self.zoom) / self.zoom < 0.01 && !g.MOUSE_PRESSED) gn.want(s.contig.name, start, end);
      for (var i = i0; i < i1; i++) {
        var x = sc.x + (cols.pos[i] - 0.5) / s.contig.length * sc.w;
        var z = cols.zyg[i], t = cols.type[i];
        var h = z === Z.HOM ? 20 : z === Z.HET ? 12 : 6;
        var ci = t >= 5 ? 5 : Math.min(t, 4);
        ctx.fillStyle = TICK_FILL[ci][cols.pass[i]];
        if (z === Z.MISSING || z === Z.REF) {
          ctx.strokeStyle = ctx.fillStyle; ctx.lineWidth = 1;
          ctx.strokeRect(x - 1.5, axisY - 4 - h, 3, h);
        } else {
          ctx.fillRect(x - (showText ? 1.5 : 0.5), axisY - 3 - h, showText ? 3 : 1.2, h);
        }
        if (showText && cols.labels[i]) {
          g.setText(TICK_FILL[ci][1], 10, 'Menlo, monospace', 'left', 'middle');
          ctx.save(); ctx.translate(x, axisY - 8 - h); ctx.rotate(-Math.PI / 2);
          ctx.fillText(cols.labels[i].split(' GT=')[0], 0, 0); ctx.restore();
        }
        if (gn) { // red dot: not in gnomAD; orange: rare (AF < 1%)
          var fr = gn.forVariant(s.contig.name, cols, i);
          if (fr && (fr.absent || fr.af < 0.01)) {
            ctx.fillStyle = fr.absent || fr.af < 1e-4 ? 'rgb(255,70,70)' : 'rgb(255,160,60)';
            ctx.beginPath(); ctx.arc(x, axisY - 8 - h, 2.5, 0, 2 * Math.PI); ctx.fill();
          }
        }
        var dist = Math.abs(g.mX - x);
        if (dist < bestDist && g.mY > axisY - 40 && g.mY < axisY + 4) { bestDist = dist; best = [s, i, x]; }
      }
    });
    if (best) {
      var cols = d.variants[best[0].contig.key], i = best[1];
      ctx.fillStyle = 'white'; ctx.fillRect(best[2] - 1, axisY - 26, 2, 24);
      this.hover = { lines: [
        best[0].contig.name + ':' + cols.pos[i].toLocaleString() + '  ' + G.vcf.TYPE_NAMES[cols.type[i]],
        cols.labels[i] || '(details not kept: file has more than 1.5M variants)',
        G.vcf.ZYG_NAMES[cols.zyg[i]] + (cols.pass[i] ? '' : ', failed filter'),
        this.gnomadText(best[0].contig.name, cols, i)
      ].filter(Boolean) };
    }
  };

  // Findings (ClinVar P/LP alleles in the sample, and dropped findings files):
  // a stem with a head at every zoom level, so they are never lost.
  View.prototype.drawFindings = function (g, axisY) {
    var ctx = g.context, self = this, placed = [], best = null, bestD = 9;
    var hl = this.highlightFinding;
    this.findings.forEach(function (f) {
      var x = self.bpToX(G.genome.normName(f.chrom), f.pos);
      if (x === null || x < -20 || x > g.cW + 20) return;
      var srcs = Object.keys(f.sources), both = srcs.length > 1 && srcs.every(function (k) { return f.sources[k].status === 'reported'; });
      var col = f.classification === 'Pathogenic' ? 'rgb(255,70,70)' : 'rgb(255,160,60)';
      var h = f.severity === 'Critical' ? 110 : f.severity === 'High' ? 95 : 80;
      var r = f.reported ? 5 : 3.5, y = axisY - h;
      var on = hl === f;
      ctx.strokeStyle = withAlpha(col, on ? 1 : 0.6); ctx.lineWidth = on ? 2 : 1;
      ctx.beginPath(); ctx.moveTo(x, axisY - 2); ctx.lineTo(x, y + r); ctx.stroke();
      ctx.fillStyle = f.reported ? col : BG;
      ctx.beginPath(); ctx.arc(x, y, on ? r + 2 : r, 0, 2 * Math.PI); ctx.fill();
      ctx.strokeStyle = both ? 'white' : col; ctx.lineWidth = both ? 1.5 : 1; ctx.stroke();
      // gene label if it does not collide with one already placed
      g.setText(on ? 'white' : 'rgba(255,255,255,0.8)', 11, 'Helvetica, Arial, sans-serif', 'center', 'bottom');
      var tw = g.getTextW(f.gene);
      if (on || placed.every(function (p) { return Math.abs(p[0] - x) > (p[1] + tw) / 2 + 4 || p[2] !== h; })) {
        g.fText(f.gene, x, y - r - 3); placed.push([x, tw, h]);
      }
      var dd = Math.hypot(g.mX - x, g.mY - y);
      if (dd < bestD) { bestD = dd; best = [f, x, y]; }
    });
    if (best) {
      var f = best[0];
      g.setCursor('pointer');
      var lines = [f.gene + '  ' + (f.classification || '') + (f.severity ? '  severity ' + f.severity : ''),
        (f.variant_name || '').slice(0, 90), (f.phenotype || '').split('|').slice(0, 2).join('; ').slice(0, 90),
        f.chrom + ':' + f.pos.toLocaleString() + ' ' + (f.ref || '').slice(0, 12) + '>' + (f.alt || '').slice(0, 12) + '  GT ' + (f.gt || '?') + '  ' + (f.zygosity || '')];
      Object.keys(f.sources).forEach(function (k) {
        var x = f.sources[k];
        lines.push(k + ': ' + (x.status === 'reported' ? 'reported' : 'not reported, ' + (x.reason || '')));
      });
      if (this.gnomadOn()) {
        var fr = this.gnomad.lookup(f.chrom, f.pos, f.ref, f.alt);
        lines.push(fr === undefined ? 'gnomAD: loading' : 'gnomAD v4: ' + G.gnomad.rarity(fr).text +
          (!fr.absent && fr.af > 0.05 ? '  (above 5%: ACMG BA1 would call this benign unless it is a listed exception)' : ''));
      }
      if (f.action) lines.push('Asclepius action: ' + f.action.slice(0, 90));
      lines.push('click to zoom to base level');
      this.hover = { lines: lines.filter(Boolean) };
      if (g.MOUSE_UP_FAST) this.goTo(f.chrom, f.pos - 60, f.pos + 60);
    }
  };

  // ----- genes and regulatory links (src/regulatory.js), GRCh38 only

  var REG_COLORS = { intergenic: 'rgb(255,190,70)', genic: 'rgb(90,200,255)', promoter: 'rgb(255,100,170)' };
  var REG_MAX_BP_PER_PX = 60000, GENE_MAX_BP_PER_PX = 400000, REG_MAX_DRAWN = 3000;
  var GENE_COLORS = { coding: 'rgba(255,255,255,0.6)', lncRNA: 'rgba(80,210,190,0.75)', smallRNA: 'rgba(255,110,190,0.9)', pseudogene: 'rgba(255,255,255,0.14)', other: 'rgba(255,255,255,0.22)' };
  G.GENE_COLORS = GENE_COLORS;
  G.REG_COLORS = REG_COLORS;

  View.prototype.grch38 = function () { return this.data && this.data.build === 'GRCh38'; };
  View.prototype.genesOn = function () { return this.layers.genes && this.genes && this.grch38(); };
  View.prototype.regOn = function () { return this.layers.regulatory && this.reg && this.reg.links && this.reg.links.length && this.grch38(); };

  // The gene whose links are lit: the clicked one, else the hovered one.
  View.prototype.litGene = function () { return this.focusGene || this.hoverGene || null; };

  View.prototype.drawGenes = function (g, axisY) {
    var ctx = g.context, self = this, y = axisY + 44, placed = [], over = null;
    var lit = this.litGene();
    this.hoverGene = null;
    this.segs.forEach(function (s) {
      var sc = self.segScreen(s);
      if (sc.x > g.cW || sc.x + sc.w < 0) return;
      var bpPerPx = s.contig.length / sc.w;
      if (bpPerPx > GENE_MAX_BP_PER_PX) return;
      var a = Math.max(1, (0 - sc.x) * bpPerPx), b = (g.cW - sc.x) * bpPerPx;
      var list = self.genes.inRange(s.contig.key, a, b);
      if (list.length > 6000) list = list.filter(function (x) { return x.type === 'protein_coding'; });
      list.forEach(function (gn) {
        if ((gn.cls === 'lncRNA' || gn.cls === 'smallRNA' || gn.cls === 'pseudogene') && !self.layers[gn.cls] && lit !== gn.name) return;
        var x0 = sc.x + (gn.start - 0.5) / s.contig.length * sc.w, x1 = sc.x + (gn.end + 0.5) / s.contig.length * sc.w;
        var coding = gn.type === 'protein_coding', on = lit === gn.name;
        ctx.fillStyle = on ? 'white' : GENE_COLORS[gn.cls] || GENE_COLORS.other;
        ctx.fillRect(x0, y, Math.max(1, x1 - x0), on ? 4 : 3);
        var xt = gn.strand > 0 ? x0 : x1; // TSS tick, with a strand arrow when zoomed in
        ctx.fillRect(xt - 0.5, y - 4, 1, 4);
        if (bpPerPx < 3000) { ctx.beginPath(); ctx.moveTo(xt, y - 4); ctx.lineTo(xt + 5 * gn.strand, y - 2); ctx.lineTo(xt, y); ctx.fill(); }
        if (g.mY > y - 6 && g.mY < y + 16 && g.mX >= Math.min(x0, x1) - 2 && g.mX <= Math.max(x0, x1) + 2) {
          if (!over || (x1 - x0) < over[1]) over = [gn, x1 - x0];
        }
        if ((coding || on || (gn.cls === 'lncRNA' && bpPerPx < 3000) || (gn.cls === 'smallRNA' && bpPerPx < 300)) && bpPerPx < 30000) {
          g.setText(on ? 'white' : 'rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'left', 'top');
          var tw = g.getTextW(gn.name), lx = Math.max(x0, Math.min(x1 - tw, 4));
          if (on) { // the lit gene's name sits on a dark plate so neighbours do not run into it
            ctx.fillStyle = BG; ctx.fillRect(lx - 3, y + 4, tw + 6, 13);
            g.setText('white', 10, 'Helvetica, Arial, sans-serif', 'left', 'top');
            g.fText(gn.name, lx, y + 5); placed.push([lx, lx + tw]);
          } else if (placed.every(function (p) { return lx > p[1] + 6 || lx + tw < p[0] - 6; })) { g.fText(gn.name, lx, y + 5); placed.push([lx, lx + tw]); }
        }
      });
    });
    if (over) {
      var gn = over[0], links = this.reg && this.reg.byGene ? this.reg.byGene.get(gn.name) || [] : [];
      this.hoverGene = gn.name;
      var withVar = links.filter(function (l) { return l.variants; }).length;
      var asPairs = (this.genes.pairsByKey && this.genes.pairsByKey[gn.chrom ? G.genome.normName(gn.chrom) : ''] || []).filter(function (p) { return p.lnc === gn || p.gene === gn; });
      this.hover = { lines: [gn.name + '  ' + gn.type.replace(/_/g, ' ') + '  ' + (gn.strand > 0 ? '+' : '-') + ' strand' +
        (asPairs.length ? '  |  antisense pair with ' + asPairs.map(function (p) { return (p.lnc === gn ? p.gene.name : p.lnc.name) + (p.coord ? ' (GeneChords rho ' + (+p.coord.rho).toFixed(2) + ', ' + p.coord.marks + ')' : ''); }).join(', ') : ''),
        gn.chrom + ':' + gn.start.toLocaleString() + '-' + gn.end.toLocaleString() + '  (' + fmtBp(gn.end - gn.start + 1) + ')',
        links.length ? links.length + ' regulatory links in the picked tissues' + (withVar ? ', ' + withVar + ' with a sample variant' : '') : 'no regulatory links in the picked tissues',
        this.focusGene === gn.name ? 'click to release' : 'click to keep its links lit'] };
      g.setCursor('pointer');
      if (g.MOUSE_UP_FAST) { this.focusGene = this.focusGene === gn.name ? null : gn.name; if (G.app.onFocusGene) G.app.onFocusGene(this.focusGene); }
    }
  };

  View.prototype.drawRegulatory = function (g, axisY) {
    var ctx = g.context, self = this, reg = this.reg, maxRy = axisY - 40;
    var lit = this.litGene(), drawn = [], elemsOver = null;
    this.segs.forEach(function (s) {
      var sc = self.segScreen(s);
      if (sc.x > g.cW || sc.x + sc.w < 0) return;
      var bpPerPx = s.contig.length / sc.w;
      if (bpPerPx > REG_MAX_BP_PER_PX && !lit) return;
      var a = Math.max(1, (0 - sc.x) * bpPerPx), b = (g.cW - sc.x) * bpPerPx;
      var list = reg.linksInRange(s.contig.key, a, b);
      if (lit) list = list.filter(function (l) { return l.gene === lit; }).concat(bpPerPx <= REG_MAX_BP_PER_PX ? list.filter(function (l) { return l.gene !== lit; }) : []);
      list.forEach(function (l) { drawn.push([l, s, sc]); });
    });
    // strongest first when there are too many; the lit gene always makes it
    if (drawn.length > REG_MAX_DRAWN) {
      drawn.sort(function (p, q) { return (q[0].gene === lit) - (p[0].gene === lit) || q[0].score - p[0].score; });
      drawn = drawn.slice(0, REG_MAX_DRAWN);
    }
    this.regVisible = drawn.length;
    var best = null, bestErr = 4;
    drawn.forEach(function (d) {
      var l = d[0], s = d[1], sc = d[2], L = s.contig.length;
      var xe = sc.x + (l.mid - 0.5) / L * sc.w, xt = sc.x + (l.tss - 0.5) / L * sc.w;
      var on = !lit || l.gene === lit, col = REG_COLORS[l.cls] || REG_COLORS.intergenic;
      // element box on the axis; white outline when the sample has a variant in it
      var ex0 = sc.x + (l.start - 0.5) / L * sc.w, ew = Math.max(1.5, (l.end - l.start + 1) / L * sc.w);
      ctx.fillStyle = withAlpha(col, on ? 0.9 : 0.15);
      ctx.fillRect(ex0, axisY - 5, ew, 4);
      if (l.variants && on) {
        var rare = self.linkRare(l);
        ctx.strokeStyle = rare ? 'rgb(255,70,70)' : 'white'; ctx.lineWidth = rare ? 2 : 1;
        ctx.strokeRect(ex0 - 1, axisY - 6, ew + 2, 6);
      }
      if (g.mY > axisY - 8 && g.mY < axisY && g.mX >= ex0 - 2 && g.mX <= ex0 + ew + 2) elemsOver = l;
      if (l.self) return; // a gene's own promoter: box only, no arc
      var lx = Math.min(xe, xt), rx = Math.max(xe, xt), r = Math.max(1.5, (rx - lx) / 2);
      if (r > 60000) return;
      var ry = Math.min(r, maxRy * (0.25 + 0.5 * Math.min(1, r / (g.cW * 0.5))));
      var nT = Object.keys(l.tissues).length, abcOnly = l.scores.e2g === undefined;
      ctx.strokeStyle = withAlpha(col, on ? 0.25 + 0.7 * l.score : 0.05);
      ctx.lineWidth = (on && lit ? 1.6 : 0.8) + 0.6 * (nT - 1) + (l.agree ? 0.9 : 0);
      if (abcOnly) ctx.setLineDash([4, 3]); // ABC only: dashed; ENCODE-rE2G: solid; both: solid and thicker
      ctx.beginPath(); ctx.ellipse((lx + rx) / 2, axisY - 1, r, ry, 0, Math.PI, 2 * Math.PI); ctx.stroke();
      if (abcOnly) ctx.setLineDash([]);
      if (on && g.mY < axisY - 2) {
        var nx = (g.mX - (lx + rx) / 2) / r, ny = (g.mY - axisY) / ry;
        var err = Math.abs(Math.sqrt(nx * nx + ny * ny) - 1) * Math.min(r, ry);
        if (err < bestErr) { bestErr = err; best = [l, (lx + rx) / 2, r, ry]; }
      }
    });
    var show = best ? best[0] : elemsOver;
    if (best) {
      ctx.strokeStyle = 'white'; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.ellipse(best[1], axisY - 1, best[2], best[3], 0, Math.PI, 2 * Math.PI); ctx.stroke();
    }
    if (show) this.hover = { lines: this.linkLines(show) };
    if (drawn.length === 0 && reg.links.length && !lit) {
      g.setText('rgba(255,255,255,0.35)', 10, 'Helvetica, Arial, sans-serif', 'right', 'bottom');
      g.fText('zoom in below ' + fmtBp(REG_MAX_BP_PER_PX) + ' per pixel for enhancer arcs, or type a gene name', g.cW - 12, axisY - 8);
    }
  };

  View.prototype.linkLines = function (l) {
    var self = this, d = this.data, cols = d && d.variants ? d.variants[l.key] : null;
    var tissues = Object.keys(l.tissues).map(function (t) {
      var ts = l.tissues[t];
      return t + ' ' + ['e2g', 'abc'].filter(function (k) { return ts[k] !== undefined; }).map(function (k) { return (k === 'e2g' ? 'rE2G ' : 'ABC ') + ts[k].toFixed(2); }).join('/');
    });
    var models = [];
    if (l.scores.e2g !== undefined) models.push('ENCODE-rE2G ' + l.scores.e2g.toFixed(2) + ' (' + l.sets.e2g + ' set' + (l.sets.e2g > 1 ? 's' : '') + ')');
    if (l.scores.abc !== undefined) models.push('ABC ' + l.scores.abc.toFixed(3) + ' (' + l.sets.abc + ' set' + (l.sets.abc > 1 ? 's' : '') + ')');
    var lines = [(l.self ? 'own promoter of ' : l.cls === 'promoter' ? 'promoter element regulating ' : l.cls === 'genic' ? 'enhancer in a gene body, regulating ' : 'distal enhancer regulating ') + l.gene,
      l.chrom + ':' + l.start.toLocaleString() + '-' + l.end.toLocaleString() + (l.self ? '' : '  ' + fmtBp(Math.abs(l.tss - l.mid)) + ' from the TSS'),
      models.join(', ') + (l.agree ? ': both models agree' : ''),
      'tissues: ' + tissues.slice(0, 4).join(', ') + (tissues.length > 4 ? ' +' + (tissues.length - 4) : '')];
    if (l.variants && cols) {
      lines.push(l.variants.length + ' sample variant' + (l.variants.length > 1 ? 's' : '') + ' in this element:');
      l.variants.slice(0, 4).forEach(function (i) {
        var al = cols.alleles(i);
        lines.push('  ' + cols.pos[i].toLocaleString() + '  ' + (al ? al.ref.slice(0, 8) + '>' + al.alts.join(',').slice(0, 12) : G.vcf.TYPE_NAMES[cols.type[i]]) + '  ' +
          G.vcf.ZYG_NAMES[cols.zyg[i]] + (self.gnomadOn() ? '  ' + self.gnomadText(l.chrom, cols, i).replace('gnomAD v4: ', '') : ''));
      });
    }
    return lines;
  };

  // ----- gnomAD (src/gnomad.js)

  var GNOMAD_MAX_VIEW = 100000; // ask for the visible region only below this span (4 tiles)
  View.prototype.gnomadOn = function () { return this.gnomad && this.gnomad.enabled && this.grch38(); };

  View.prototype.gnomadText = function (chrom, cols, i) {
    if (!this.gnomadOn()) return '';
    var r = this.gnomad.forVariant(chrom, cols, i);
    return r === undefined ? 'gnomAD: not loaded (zoom below ' + fmtBp(GNOMAD_MAX_VIEW) + ' or light its gene)' : 'gnomAD v4: ' + G.gnomad.rarity(r).text;
  };

  // True when a sample variant in the link's element is rare (< 1%) or absent
  // from gnomAD. Cached per link until new tiles arrive.
  View.prototype.linkRare = function (l) {
    if (!this.gnomadOn() || !l.variants) return false;
    var gn = this.gnomad;
    if (l.rareVersion === gn.version) return l.rare;
    var cols = this.data.variants[l.key], rare = false;
    for (var k = 0; k < l.variants.length && !rare; k++) {
      var r = gn.forVariant(l.chrom, cols, l.variants[k]);
      if (r && (r.absent || r.af < 0.01)) rare = true;
    }
    l.rare = rare; l.rareVersion = gn.version;
    return rare;
  };

  // ----- GWAS markers and antisense chords

  View.prototype.gwasTrack = function (key, length, binSize) {
    this.gwasTracks = this.gwasTracks || {};
    var id = key + ':' + binSize;
    if (this.gwasTracks[id]) return this.gwasTracks[id];
    var g = this.gwas.byContig[key];
    if (!g) return null;
    var t = new G.genome.Track(length, binSize);
    for (var i = 0; i < g.n; i++) t.add(g.pos[i], 1);
    t.buildPyramid();
    return (this.gwasTracks[id] = t);
  };

  var DOSE_COLORS = ['rgba(190,140,255,0.35)', 'rgb(190,140,255)', 'rgb(235,120,255)'];

  // SNPs as small diamonds below the ticks when zoomed in; the chosen trait's
  // loci stand on stems at any zoom, like findings.
  View.prototype.drawGwas = function (g, axisY) {
    var ctx = g.context, self = this, gw = this.gwas, d = this.data, best = null, bestD = 6;
    var trait = G.app.gwasTraitLoci || null, traitSet = trait ? trait.set : null;
    var diamond = function (x, y, r, fill, stroke) {
      ctx.beginPath(); ctx.moveTo(x, y - r); ctx.lineTo(x + r, y); ctx.lineTo(x, y + r); ctx.lineTo(x - r, y); ctx.closePath();
      if (fill) { ctx.fillStyle = fill; ctx.fill(); } if (stroke) { ctx.strokeStyle = stroke; ctx.lineWidth = 1; ctx.stroke(); }
    };
    this.segs.forEach(function (s) {
      var sc = self.segScreen(s);
      if (sc.x > g.cW || sc.x + sc.w < 0) return;
      var bpPerPx = s.contig.length / sc.w, k = s.contig.key, gc = gw.byContig[k];
      if (!gc) return;
      var a = Math.max(1, (0 - sc.x) * bpPerPx), b = (g.cW - sc.x) * bpPerPx;
      var idx = bpPerPx <= 5000 ? gw.inRange(k, a, b) : [];
      if (idx.length > 3000) idx = [];
      idx.forEach(function (i) {
        var x = sc.x + (gc.pos[i] - 0.5) / s.contig.length * sc.w, dz = gw.dosage(d, k, i).dosage;
        var y = axisY - 34;
        diamond(x, y, 3.5, dz === null ? null : DOSE_COLORS[dz], dz === null ? 'rgba(190,140,255,0.6)' : null);
        var dd = Math.hypot(g.mX - x, g.mY - y);
        if (dd < bestD) { bestD = dd; best = [k, i, x, y]; }
      });
      if (traitSet) trait.loci.forEach(function (l) {
        if (l.k !== k) return;
        var x = sc.x + (gc.pos[l.i] - 0.5) / s.contig.length * sc.w;
        if (x < -10 || x > g.cW + 10) return;
        var dz = gw.dosage(d, k, l.i).dosage, y = axisY - 130;
        ctx.strokeStyle = 'rgba(190,140,255,0.6)'; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(x, axisY - 2); ctx.lineTo(x, y); ctx.stroke();
        diamond(x, y, 6, dz === null ? '#141414' : DOSE_COLORS[dz], 'rgb(190,140,255)');
        var dd = Math.hypot(g.mX - x, g.mY - y);
        if (dd < bestD + 3) { bestD = dd; best = [k, l.i, x, y]; }
      });
    });
    if (best) {
      var k = best[0], i = best[1], gc = gw.byContig[k], dz = gw.dosage(d, k, i), as = gw.associations(k, i);
      ctx.strokeStyle = 'white'; ctx.lineWidth = 2; diamond(best[2], best[3], 7, null, 'white');
      var lines = [gc.rsid[i] + '  risk allele ' + gc.risk[i] + (gc.freq[i] != null ? ' (freq ' + gc.freq[i] + ')' : '') + '  near ' + (gc.gene[i] || '?'),
        'this sample: ' + (dz.dosage === null ? 'unknown (' + dz.reason + ')' : dz.dosage + ' cop' + (dz.dosage === 1 ? 'y' : 'ies') + ' of the risk allele' + (dz.reason ? ' (' + dz.reason + ')' : ''))];
      as.slice(0, 5).forEach(function (x) { lines.push('p ' + (x.p > 0 ? x.p.toExponential(0) : '< 1e-300') + '  ' + x.trait.slice(0, 70) + (x.effect ? '  (OR/beta ' + x.effect + ')' : '')); });
      if (as.length > 5) lines.push('+' + (as.length - 5) + ' more associations');
      var reg = G.app.reg;
      if (reg && reg.links) {
        var inEl = reg.linksInRange(k, gc.pos[i], gc.pos[i]).filter(function (l) { return l.start <= gc.pos[i] && l.end >= gc.pos[i] && !l.self; });
        if (inEl.length) lines.push('inside a regulatory element linked to ' + inEl.map(function (l) { return l.gene; }).filter(function (x, j, a2) { return a2.indexOf(x) === j; }).join(', '));
      }
      this.hover = { lines: lines };
    }
  };

  // lncRNA / protein-coding antisense pairs as small teal arcs between their
  // start sites; pairs in the GeneChords results are gold and thicker.
  View.prototype.drawAntisense = function (g, axisY) {
    var ctx = g.context, self = this, genes = this.genes, best = null, bestErr = 4;
    if (!genes.pairsByKey) return;
    this.segs.forEach(function (s) {
      var sc = self.segScreen(s), bpPerPx = s.contig.length / sc.w;
      if (bpPerPx > 20000 || sc.x > g.cW || sc.x + sc.w < 0) return;
      (genes.pairsByKey[s.contig.key] || []).forEach(function (p) {
        var x0 = sc.x + (p.lnc.tss - 0.5) / s.contig.length * sc.w, x1 = sc.x + (p.gene.tss - 0.5) / s.contig.length * sc.w;
        if (Math.max(x0, x1) < 0 || Math.min(x0, x1) > g.cW) return;
        var r = Math.max(2, Math.abs(x1 - x0) / 2), cx = (x0 + x1) / 2, ry = Math.min(r, 40);
        ctx.strokeStyle = p.coord ? 'rgba(255,205,80,0.95)' : 'rgba(80,210,190,0.55)'; ctx.lineWidth = p.coord ? 2.5 : 1;
        ctx.beginPath(); ctx.ellipse(cx, axisY + 44, r, ry, 0, Math.PI, 2 * Math.PI); ctx.stroke();
        if (g.mY < axisY + 44 && g.mY > axisY + 44 - ry - 6) {
          var nx = (g.mX - cx) / r, ny = (g.mY - axisY - 44) / ry, err = Math.abs(Math.sqrt(nx * nx + ny * ny) - 1) * Math.min(r, ry);
          if (err < bestErr) { bestErr = err; best = p; }
        }
      });
    });
    if (best) this.hover = { lines: ['antisense pair: ' + best.lnc.name + ' (lncRNA, ' + (best.lnc.strand > 0 ? '+' : '-') + ') and ' + best.gene.name + ' (' + (best.gene.strand > 0 ? '+' : '-') + ')',
      'overlap ' + fmtBp(best.overlap) + ', start sites ' + fmtBp(Math.abs(best.lnc.tss - best.gene.tss)) + ' apart',
      best.coord ? 'GeneChords: coordinated in ' + best.coord.marks + ', rho ' + (+best.coord.rho).toFixed(3) + ', AUC ' + best.coord.auc : 'not among the GeneChords coordination results'] };
  };

  View.prototype.drawReadout = function (g, axisY) {
    // Cursor genome position and scale, bottom left.
    var s = this.segs.find(function (s) { return true; }), self = this, txt = '';
    for (var i = 0; i < this.segs.length; i++) {
      var sc = this.segScreen(this.segs[i]);
      if (g.mX >= sc.x && g.mX <= sc.x + sc.w) {
        var c = this.segs[i].contig, bp = Math.max(1, Math.ceil((g.mX - sc.x) / sc.w * c.length));
        txt = c.name + ':' + bp.toLocaleString() + '   1 px = ' + fmtBp(c.length / sc.w);
        break;
      }
    }
    g.setText('rgba(255,255,255,0.45)', 11, 'Menlo, monospace', 'left', 'bottom');
    g.fText(txt, 12, g.cH - 38);
  };

  View.prototype.drawTooltip = function (g, lines) {
    var ctx = g.context;
    g.setText('white', 12, 'Helvetica, Arial, sans-serif', 'left', 'top');
    var w = Math.max.apply(null, lines.map(function (l) { return g.getTextW(l); })) + 16, h = lines.length * 16 + 10;
    var x = Math.min(g.mX + 14, g.cW - w - 6), y = g.mY + 18 + h > g.cH ? g.mY - h - 10 : g.mY + 18;
    ctx.fillStyle = 'rgba(0,0,0,0.85)'; ctx.fillRect(x, y, w, h);
    ctx.strokeStyle = 'rgba(255,255,255,0.25)'; ctx.lineWidth = 1; ctx.strokeRect(x + 0.5, y + 0.5, w - 1, h - 1);
    lines.forEach(function (l, i) {
      g.setText(i === 0 ? 'white' : 'rgba(255,255,255,0.75)', 12, 'Helvetica, Arial, sans-serif', 'left', 'top');
      g.fText(l, x + 8, y + 6 + i * 16);
    });
  };

  View.prototype.drawMessage = function (g, text) {
    g.setText('rgba(255,255,255,0.6)', 15, 'Helvetica, Arial, sans-serif', 'center', 'middle');
    g.fText(text, g.cX, g.cY);
  };

  // Idle screen: a slow arc animation with a hint.
  View.prototype.drawEmpty = function (g) {
    var ctx = g.context, t = Date.now() / 1000, y = g.cH * 0.6;
    ctx.strokeStyle = 'rgba(255,255,255,0.4)'; ctx.fillStyle = 'rgba(255,255,255,0.4)';
    ctx.fillRect(g.cW * 0.1, y, g.cW * 0.8, 1);
    for (var i = 0; i < 14; i++) {
      var a = g.cW * (0.1 + 0.8 * ((Math.sin(i * 12.9898) * 43758.5453) % 1 + 1) % 1);
      var b = g.cW * (0.1 + 0.8 * ((Math.sin(i * 78.233) * 12345.678) % 1 + 1) % 1);
      var r = Math.abs(b - a) / 2 * (0.9 + 0.1 * Math.sin(t + i));
      ctx.strokeStyle = withAlpha(['rgb(120,180,255)', 'rgb(235,90,200)', 'rgb(255,170,80)'][i % 3], 0.35);
      ctx.lineWidth = 1 + (i % 3);
      ctx.beginPath(); ctx.ellipse((a + b) / 2, y, r, Math.min(r, y - 60), 0, Math.PI, 2 * Math.PI); ctx.stroke();
    }
    g.setText('rgba(255,255,255,0.75)', 16, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText('Drop a VCF, gVCF (.vcf / .vcf.gz) or BAM file anywhere', g.cX, y + 24);
    g.setText('rgba(255,255,255,0.4)', 12, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText('everything is read in this browser tab; nothing is uploaded', g.cX, y + 48);
  };

  // Unaligned BAM (e.g. basecaller output): read length and quality only.
  View.prototype.drawUnaligned = function (g) {
    var d = this.data, ctx = g.context;
    var hist = d.lenHist, n = hist.length, max = Math.max.apply(null, Array.prototype.slice.call(hist)) || 1;
    var left = 60, right = g.cW - 60, base = g.cH * 0.6, H = g.cH * 0.35;
    var used = [];
    for (var i = 0; i < n; i++) if (hist[i] > 0) used.push(i);
    var b0 = used[0] || 0, b1 = used[used.length - 1] || n - 1, bw = (right - left) / (b1 - b0 + 1);
    ctx.fillStyle = 'rgba(255,255,255,0.4)'; ctx.fillRect(left, base, right - left, 1);
    var over = null;
    for (var b = b0; b <= b1; b++) {
      var h = Math.sqrt(hist[b] / max) * H, x = left + (b - b0) * bw;
      var isOver = g.mX >= x && g.mX < x + bw && g.mY < base && g.mY > base - H - 20;
      ctx.fillStyle = isOver ? 'white' : 'rgba(90,210,190,0.85)';
      ctx.fillRect(x + 1, base - h, Math.max(1, bw - 2), h);
      if (isOver) over = b;
    }
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'center', 'top');
    for (b = b0; b <= b1; b += Math.max(1, Math.round((b1 - b0) / 8))) {
      g.fText(fmtBp(d.lenBinToLength(b)), left + (b - b0 + 0.5) * bw, base + 8);
    }
    g.setText('rgba(255,255,255,0.8)', 14, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText('Unaligned reads: no reference coordinates, so no genome layout.', g.cX, base - H - 70);
    g.setText('rgba(255,255,255,0.5)', 12, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText(d.stats.reads.toLocaleString() + ' reads, ' + fmtBp(d.stats.bases).replace('bp', 'bases') +
      '. Read length distribution (log scale, sqrt height). Align the reads (e.g. minimap2) to see them on the genome.', g.cX, base - H - 48);

    // quality histogram, small, bottom right
    var q = d.qualHist, qmax = Math.max.apply(null, Array.prototype.slice.call(q)) || 1, qx = left, qy = Math.min(g.cH - 40, base + 130);
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
    g.fText('mean base quality per read (Phred)', qx, qy - 64);
    for (var k = 0; k <= 60; k++) {
      var qh = Math.sqrt(q[k] / qmax) * 50;
      ctx.fillStyle = 'rgba(255,170,80,0.85)'; ctx.fillRect(qx + k * 4, qy - qh, 3, qh);
    }
    g.setText('rgba(255,255,255,0.4)', 9, 'Helvetica, Arial, sans-serif', 'center', 'top');
    [0, 10, 20, 30, 40].forEach(function (v) { g.fText(String(v), qx + v * 4 + 1.5, qy + 4); });
    if (over !== null) {
      this.drawTooltip(g, [fmtBp(d.lenBinToLength(over)) + ' to ' + fmtBp(d.lenBinToLength(over + 1)), hist[over].toLocaleString() + ' reads']);
    }
  };

  View.prototype.drawStatus = function (g) {
    var s = this.status, w = Math.min(420, g.cW - 40), x = g.cX - w / 2, y = 24;
    g.context.fillStyle = 'rgba(255,255,255,0.12)'; g.context.fillRect(x, y, w, 4);
    g.context.fillStyle = 'rgb(90,210,190)'; g.context.fillRect(x, y, w * Math.min(1, s.fraction || 0), 4);
    g.setText('rgba(255,255,255,0.7)', 11, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText(s.text, g.cX, y + 10);
  };

  G.View = View;
  G.COLORS = COLORS;
  G.fmtBp = fmtBp;
})(globalThis.G = globalThis.G || {});
