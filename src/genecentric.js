/*
 * Gene view: one gene's regulatory neighbourhood, the /attention "pressed
 * word" turned radial. The gene sits in the middle; each regulatory element
 * linked to it (picked tissues and models) is a node around it:
 *   left half upstream, right half downstream of the TSS (strand aware),
 *   distance from the centre on a log scale of bp from the TSS,
 *   colour by element class, size by score, white ring when the sample has
 *   a variant in it, red ring when that variant is rare (gnomAD on).
 * Edges follow the Arcs convention: solid ENCODE-rE2G, dashed ABC only,
 * thicker when both models agree. On the right, a tissue x element heatmap
 * of scores shows which tissues use which element.
 */
(function (G) {
  var TAU = Math.PI * 2;
  var COLORS = { intergenic: 'rgb(255,190,70)', genic: 'rgb(90,200,255)', promoter: 'rgb(255,100,170)' };

  function GeneView() {}

  GeneView.prototype.draw = function (g) {
    var ctx = g.context, view = G.app.view, reg = G.app.reg, gene = view.focusGene;
    var msg = function (t) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(t, g.cX, g.cY); };
    if (!view.data || view.data.build !== 'GRCh38') return msg('The gene view needs a GRCh38 file.');
    if (!reg || !reg.byGene) return msg('Pick tissues (Regulatory, tissues) to load enhancer links.');
    if (!gene) return msg('Type a gene in search (Ctrl+K), or click one in the Arcs gene lane.');
    var links = (reg.byGene.get(gene) || []).slice();
    var gn = view.genes && view.genes.get(gene);
    if (!links.length) return msg(gene + ' has no regulatory links in the picked tissues.');
    var strand = gn ? gn.strand : 1, tss = links[0].tss;
    var heatW = Math.min(420, g.cW * 0.32), cx = (g.cW - heatW) / 2, cy = g.cY + 20, R = Math.min(cx - 60, g.cH / 2 - 90);
    var maxD = links.reduce(function (m, l) { return Math.max(m, Math.abs(l.mid - tss)); }, 1000);
    var rOf = function (dist) { return 40 + (R - 40) * Math.log10(1 + Math.abs(dist)) / Math.log10(1 + maxD); };

    // place nodes: upstream on the left half, downstream on the right, spread by rank
    var up = [], down = [];
    links.forEach(function (l) { ((l.mid - tss) * strand < 0 ? up : down).push(l); });
    var place = function (list, a0, a1) {
      list.sort(function (a, b) { return Math.abs(a.mid - tss) - Math.abs(b.mid - tss); });
      list.forEach(function (l, i) {
        var t = (i + 0.5) / list.length, a = a0 + (a1 - a0) * t;
        var r = l.self ? 0 : rOf(l.mid - tss);
        l._x = cx + Math.cos(a) * r; l._y = cy + Math.sin(a) * r;
      });
    };
    place(up, Math.PI * 0.6, Math.PI * 1.4);
    place(down, -Math.PI * 0.4, Math.PI * 0.4);

    // distance rings
    [1e3, 1e4, 1e5, 1e6].forEach(function (d) {
      if (d > maxD * 1.2) return;
      ctx.strokeStyle = 'rgba(255,255,255,0.08)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(cx, cy, rOf(d), 0, TAU); ctx.stroke();
      g.setText('rgba(255,255,255,0.35)', 9, 'Helvetica, Arial, sans-serif', 'center', 'bottom');
      g.fText(G.fmtBp(d), cx, cy - rOf(d) - 2);
    });
    g.setText('rgba(255,255,255,0.4)', 11, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText('upstream', cx - R * 0.8, cy + R + 8); g.fText('downstream', cx + R * 0.8, cy + R + 8);

    var over = null, bestD = 10;
    links.forEach(function (l) {
      if (l.self) return;
      var abcOnly = l.scores.e2g === undefined;
      ctx.strokeStyle = COLORS[l.cls].replace('rgb(', 'rgba(').replace(')', ',' + (0.25 + 0.6 * l.score).toFixed(2) + ')');
      ctx.lineWidth = 1 + (l.agree ? 1.2 : 0) + 0.5 * (Object.keys(l.tissues).length - 1);
      if (abcOnly) ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(l._x, l._y); ctx.stroke();
      if (abcOnly) ctx.setLineDash([]);
    });
    links.forEach(function (l) {
      var rad = l.self ? 22 : 4 + 9 * l.score, rare = view.linkRare && view.linkRare(l);
      ctx.fillStyle = COLORS[l.cls];
      ctx.beginPath(); ctx.arc(l._x, l._y, rad, 0, TAU); ctx.fill();
      if (l.variants) { ctx.strokeStyle = rare ? 'rgb(255,60,60)' : 'white'; ctx.lineWidth = rare ? 3 : 1.5; ctx.beginPath(); ctx.arc(l._x, l._y, rad + 3, 0, TAU); ctx.stroke(); }
      var dd = Math.hypot(g.mX - l._x, g.mY - l._y) - rad;
      if (dd < bestD) { bestD = dd; over = l; }
    });
    g.setText('white', 13, 'Helvetica, Arial, sans-serif', 'center', 'middle');
    g.fText(gene, cx, cy);
    g.setText('rgba(255,255,255,0.5)', 10, 'Helvetica, Arial, sans-serif', 'center', 'top');
    g.fText((gn ? gn.type.replace(/_/g, ' ') + ', ' + (strand > 0 ? '+' : '-') + ' strand, ' : '') + links.length + ' links', cx, cy + 26);

    // heatmap: rows = elements by distance, columns = tissues
    var tissues = [];
    links.forEach(function (l) { Object.keys(l.tissues).forEach(function (t) { if (tissues.indexOf(t) < 0) tissues.push(t); }); });
    var rows = links.slice().sort(function (a, b) { return (a.mid - tss) * strand - (b.mid - tss) * strand; });
    var hx = g.cW - heatW - 20, hy = 170, cw = Math.max(14, Math.min(40, (heatW - 150) / Math.max(1, tissues.length))), rh = Math.max(6, Math.min(16, (g.cH - hy - 90) / rows.length));
    g.setText('rgba(255,255,255,0.7)', 10, 'Helvetica, Arial, sans-serif', 'left', 'bottom');
    tissues.forEach(function (t, j) {
      ctx.save(); ctx.translate(hx + 150 + j * cw + cw / 2, hy - 4); ctx.rotate(-Math.PI / 4); ctx.fillText(t.length > 22 ? t.slice(0, 20) + '..' : t, 0, 0); ctx.restore();
    });
    rows.forEach(function (l, i) {
      var y = hy + i * rh, dist = (l.mid - tss) * strand;
      g.setText(l === over ? 'white' : 'rgba(255,255,255,0.6)', Math.min(10, rh), 'Helvetica, Arial, sans-serif', 'right', 'middle');
      g.fText((l.self ? 'own promoter' : (dist < 0 ? '-' : '+') + G.fmtBp(Math.abs(dist))) + (l.variants ? ' *' : ''), hx + 144, y + rh / 2);
      tissues.forEach(function (t, j) {
        var ts = l.tissues[t], v = ts ? (ts.e2g !== undefined ? ts.e2g : Math.min(1, Math.sqrt(ts.abc / 0.25))) : 0;
        ctx.fillStyle = ts ? COLORS[l.cls].replace('rgb(', 'rgba(').replace(')', ',' + (0.2 + 0.8 * v).toFixed(2) + ')') : 'rgba(255,255,255,0.04)';
        ctx.fillRect(hx + 150 + j * cw + 1, y + 1, cw - 2, rh - 2);
      });
      if (g.mX > hx && g.mX < hx + 150 + tissues.length * cw && g.mY >= y && g.mY < y + rh) over = l;
    });
    g.setText('rgba(255,255,255,0.45)', 10, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('rows: elements by distance from the TSS (* = sample variant); cells: score per tissue', hx, hy + rows.length * rh + 6);

    if (over) {
      if (over._x) { ctx.strokeStyle = 'white'; ctx.lineWidth = 2; ctx.beginPath(); ctx.arc(over._x, over._y, (over.self ? 22 : 4 + 9 * over.score) + 6, 0, TAU); ctx.stroke(); }
      view.drawTooltip(g, view.linkLines(over).concat(['click to open in Arcs']));
      g.setCursor('pointer');
      if (g.MOUSE_UP_FAST) { G.app.setMode('arcs'); view.goTo(over.chrom, Math.min(over.start, over.tss) - 1000, Math.max(over.end, over.tss) + 1000); }
    }
  };

  G.GeneView = GeneView;
})(globalThis.G = globalThis.G || {});
