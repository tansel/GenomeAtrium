/*
 * Circular genome (Circos style, after Krzywinski et al. 2009): chromosomes
 * around a ring with their cytogenetic bands, a density ring (variants, or
 * depth for BAM), a het fraction ring, SV and read-pair chords across the
 * middle, and findings outside. The familiar picture for translocations:
 * a chord between two chromosomes stands out at once.
 *
 * Drag rotates, hover a chord or chromosome for details, click a
 * chromosome to open it in Arcs.
 */
(function (G) {
  var TAU = Math.PI * 2, BINS = 720;
  var STAIN = { gneg: '#3a3a3a', gpos25: '#5a5a5a', gpos50: '#7a7a7a', gpos75: '#9a9a9a', gpos100: '#c0c0c0', acen: '#c05050', gvar: '#6d6da8', stalk: '#4a4a4a' };

  function Circos() { this.rot = -Math.PI / 2; }

  Circos.prototype.setData = function (d) {
    this.data = d; this.rings = null;
    if (!d || !d.genome || !d.genome.contigs.length) return;
    var total = d.genome.totalLength(), gap = 0.004 * TAU, n = d.genome.contigs.length;
    var avail = TAU - gap * n, a = 0, lay = {};
    d.genome.contigs.forEach(function (c, i) { var w = c.length / total * avail; lay[c.key] = { a0: a, a1: a + w, c: c, i: i }; a += w + gap; });
    this.lay = lay;
    this.buildRings();
  };

  Circos.prototype.angle = function (chrom, pos) {
    var l = this.lay[G.genome.normName(chrom)];
    return l ? this.rot + l.a0 + (pos / l.c.length) * (l.a1 - l.a0) : null;
  };

  // Per-angle-bin summaries for the rings.
  Circos.prototype.buildRings = function () {
    var d = this.data, dens = new Float32Array(BINS), het = new Float32Array(BINS).fill(-1), owner = new Int16Array(BINS).fill(-1);
    var self = this, samp = G.Hilbert.prototype.sample;
    Object.keys(this.lay).forEach(function (k) {
      var l = self.lay[k], tr = d.tracks[k];
      if (!tr) return;
      var b0 = Math.floor(l.a0 / TAU * BINS), b1 = Math.floor(l.a1 / TAU * BINS);
      for (var b = b0; b <= b1 && b < BINS; b++) {
        var p0 = Math.max(0, (b / BINS * TAU - l.a0) / (l.a1 - l.a0) * l.c.length), p1 = Math.min(l.c.length, ((b + 1) / BINS * TAU - l.a0) / (l.a1 - l.a0) * l.c.length);
        if (p1 <= p0) continue;
        var bpp = p1 - p0;
        owner[b] = l.i;
        if (d.format === 'vcf') {
          dens[b] = samp(tr.snv.levelFor(bpp), p0, p1, 'sum') + samp(tr.indel.levelFor(bpp), p0, p1, 'sum');
          var he = samp(tr.het.levelFor(bpp), p0, p1, 'sum'), ho = samp(tr.hom.levelFor(bpp), p0, p1, 'sum');
          if (he + ho >= 5) het[b] = he / (he + ho);
        } else dens[b] = samp(tr.depth.levelFor(bpp), p0, p1, 'mean');
      }
    });
    var sorted = Array.from(dens).filter(function (v) { return v > 0; }).sort(function (a, b) { return a - b; });
    this.rings = { dens: dens, het: het, owner: owner, hi: sorted.length ? sorted[Math.floor(sorted.length * 0.98)] : 1 };
  };

  Circos.prototype.draw = function (g) {
    var ctx = g.context, d = this.data, self = this;
    if (!this.rings) { g.setText('rgba(255,255,255,0.6)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText('Open a file first.', g.cX, g.cY); return; }
    if (g.MOUSE_PRESSED) this.rot += g.DX_MOUSE * 0.005;
    // In the Atrium's window there is no page chrome to leave room for: fill the height.
    var inWin = G.app.view.mode === 'atrium';
    var cx = g.cX, cy = inWin ? g.cY - 10 : g.cY + 20, R = inWin ? (g.cH - 70) * 0.43 : Math.min(g.cW, g.cH - 120) * 0.38;
    var hover = null, dist = Math.hypot(g.mX - cx, g.mY - cy), mang = Math.atan2(g.mY - cy, g.mX - cx);
    var bands = G.app.view.cytobands;

    // chromosomes with bands, and labels
    Object.keys(this.lay).forEach(function (k) {
      var l = self.lay[k], a0 = self.rot + l.a0, a1 = self.rot + l.a1, list = bands && bands[k];
      if (list) list.forEach(function (b) {
        var s0 = a0 + b[0] / l.c.length * (a1 - a0), s1 = a0 + b[1] / l.c.length * (a1 - a0);
        ctx.strokeStyle = STAIN[b[3]] || STAIN.gneg; ctx.lineWidth = 10;
        ctx.beginPath(); ctx.arc(cx, cy, R + 6, s0, s1 + 0.0005); ctx.stroke();
      });
      else { ctx.strokeStyle = l.i % 2 ? '#6a6a6a' : '#9a9a9a'; ctx.lineWidth = 10; ctx.beginPath(); ctx.arc(cx, cy, R + 6, a0, a1); ctx.stroke(); }
      var am = (a0 + a1) / 2;
      if (a1 - a0 > 0.03) {
        g.setText('rgba(255,255,255,0.75)', 11, 'Helvetica, Arial, sans-serif', 'center', 'middle');
        g.fText(l.c.name.replace(/^chr/i, ''), cx + Math.cos(am) * (R + 24), cy + Math.sin(am) * (R + 24));
      }
      var rel = ((mang - a0) % TAU + TAU) % TAU;
      if (dist > R - 4 && dist < R + 16 && rel <= a1 - a0) hover = { lines: [l.c.name + '  ' + G.fmtBp(l.c.length), 'click to open in Arcs'], contig: l.c };
    });

    // rings: density (outer, bars inward) and het fraction (inner, colour)
    var r = this.rings;
    for (var b = 0; b < BINS; b++) {
      if (r.owner[b] < 0) continue;
      var a = this.rot + (b + 0.5) / BINS * TAU, h = Math.sqrt(Math.min(1, r.dens[b] / r.hi)) * 34;
      ctx.strokeStyle = d.format === 'vcf' ? 'rgba(120,180,255,0.85)' : 'rgba(90,210,190,0.85)'; ctx.lineWidth = TAU * R / BINS * 0.8;
      ctx.beginPath(); ctx.moveTo(cx + Math.cos(a) * (R - 2), cy + Math.sin(a) * (R - 2)); ctx.lineTo(cx + Math.cos(a) * (R - 2 - h), cy + Math.sin(a) * (R - 2 - h)); ctx.stroke();
      if (r.het[b] >= 0) {
        var v = r.het[b];
        ctx.strokeStyle = 'rgb(' + Math.round(40 + 150 * v) + ',' + Math.round(30 + 110 * v) + ',' + Math.round(80 + 175 * v) + ')'; ctx.lineWidth = TAU * (R - 44) / BINS * 1.05;
        ctx.beginPath(); ctx.moveTo(cx + Math.cos(a) * (R - 40), cy + Math.sin(a) * (R - 40)); ctx.lineTo(cx + Math.cos(a) * (R - 52), cy + Math.sin(a) * (R - 52)); ctx.stroke();
      }
    }

    // chords: SV arcs and discordant pairs (not splice junctions, too local)
    var rc = R - 58, best = null, bestD = 5;
    (d.arcs || []).forEach(function (arc) {
      if (arc.type === 'junction') return;
      var a0 = self.angle(arc.c0, arc.p0), a1 = self.angle(arc.c1, arc.p1);
      if (a0 === null || a1 === null) return;
      var x0 = cx + Math.cos(a0) * rc, y0 = cy + Math.sin(a0) * rc, x1 = cx + Math.cos(a1) * rc, y1 = cy + Math.sin(a1) * rc;
      var span = Math.abs(a1 - a0) % TAU, pull = span < 0.02 ? 0.85 : 0.15; // near ends bulge out a little, far ends pass near the centre
      var qx = cx + ((x0 + x1) / 2 - cx) * pull, qy = cy + ((y0 + y1) / 2 - cy) * pull;
      var col = (G.COLORS.arc[arc.type] || 'rgb(200,200,200)');
      ctx.strokeStyle = col.replace('rgb(', 'rgba(').replace(')', ',0.7)'); ctx.lineWidth = 0.8 + Math.log2(1 + arc.support) * 0.6;
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.quadraticCurveTo(qx, qy, x1, y1); ctx.stroke();
      // hover: distance to the curve at a few points
      for (var t = 0; t <= 1; t += 0.05) {
        var bx = (1 - t) * (1 - t) * x0 + 2 * (1 - t) * t * qx + t * t * x1, by = (1 - t) * (1 - t) * y0 + 2 * (1 - t) * t * qy + t * t * y1;
        var dd = Math.hypot(g.mX - bx, g.mY - by);
        if (dd < bestD) { bestD = dd; best = [arc, x0, y0, qx, qy, x1, y1]; }
      }
    });
    if (best) {
      ctx.strokeStyle = 'white'; ctx.lineWidth = 2.5;
      ctx.beginPath(); ctx.moveTo(best[1], best[2]); ctx.quadraticCurveTo(best[3], best[4], best[5], best[6]); ctx.stroke();
      hover = { lines: [best[0].label, 'support ' + best[0].support] };
    }

    // findings outside
    (G.app.view.findings || []).forEach(function (f) {
      var a = self.angle(f.chrom, f.pos);
      if (a === null) return;
      var x = cx + Math.cos(a) * (R + 40), y = cy + Math.sin(a) * (R + 40);
      ctx.fillStyle = f.classification === 'Pathogenic' ? 'rgb(255,70,70)' : 'rgb(255,160,60)';
      ctx.beginPath(); ctx.arc(x, y, 4, 0, TAU); ctx.fill();
      g.setText('rgba(255,255,255,0.8)', 10, 'Helvetica, Arial, sans-serif', Math.cos(a) >= 0 ? 'left' : 'right', 'middle');
      g.fText(f.gene, x + Math.cos(a) * 8, y + Math.sin(a) * 8);
    });

    g.setText('rgba(255,255,255,0.5)', 11, 'Helvetica, Arial, sans-serif', 'center', 'bottom');
    g.fText('outer: bands; blue: ' + (d.format === 'vcf' ? 'variant density' : 'depth') + (d.format === 'vcf' ? '; purple: het fraction (dark = homozygous run or hemizygous)' : '') +
      '; chords: SVs and discordant pairs; dots: findings. Drag rotates.', g.cX, g.cH - 40);

    if (hover) {
      G.app.view.drawTooltip(g, hover.lines);
      if (hover.contig) { g.setCursor('pointer'); if (g.MOUSE_UP_FAST) { G.app.setMode('arcs'); G.app.view.goTo(hover.contig.name, 1, hover.contig.length); } }
    }
  };

  G.Circos = Circos;
})(globalThis.G = globalThis.G || {});
