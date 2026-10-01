/*
 * Pathways: Reactome (human) as a treemap you can drill into, coloured by
 * how a chosen gene set falls on it.
 *
 * Gene sets: the findings' genes, the genes of a GWAS trait where this
 * genome carries a risk allele, and the gene in focus. Enrichment per
 * pathway is a one-sided hypergeometric test over the genes Reactome
 * covers (P(X >= k)); with thousands of pathways tested, treat the p-values
 * as a ranking, not as findings.
 *
 * Top level: Reactome's top pathways (no parent), area = gene count. Click
 * a box to open its sub-pathways; the breadcrumb goes back. Genes of a leaf
 * pathway are listed with the set's genes marked; click one to focus it.
 */
(function (G) {
  var LF = [0];
  function lf(n) { while (LF.length <= n) LF.push(LF[LF.length - 1] + Math.log(LF.length)); return LF[n]; }
  function logChoose(n, k) { return k < 0 || k > n ? -Infinity : lf(n) - lf(k) - lf(n - k); }
  // P(X >= k) for X ~ Hypergeometric(N, K, n).
  function hyperTail(N, K, n, k) {
    if (k <= 0) return 1;
    var s = 0, denom = logChoose(N, n);
    for (var x = k; x <= Math.min(K, n); x++) s += Math.exp(logChoose(K, x) + logChoose(N - K, n - x) - denom);
    return Math.min(1, s);
  }

  // Squarified treemap (Bruls et al.) of items with .value into rect r.
  function squarify(items, r) {
    var out = [], rest = items.filter(function (x) { return x.value > 0; }).sort(function (a, b) { return b.value - a.value; });
    var total = rest.reduce(function (s, x) { return s + x.value; }, 0), x = r.x, y = r.y, w = r.w, h = r.h;
    var scale = w * h / (total || 1);
    while (rest.length) {
      var side = Math.min(w, h), row = [], best = Infinity;
      while (rest.length) {
        var cand = row.concat([rest[0]]), sum = cand.reduce(function (s, c) { return s + c.value * scale; }, 0);
        var worst = Math.max.apply(null, cand.map(function (c) { var a = c.value * scale; return Math.max(side * side * a / (sum * sum), sum * sum / (side * side * a)); }));
        if (worst > best && row.length) break;
        best = worst; row = cand; rest.shift();
      }
      var rs = row.reduce(function (s, c) { return s + c.value * scale; }, 0), thick = rs / side, off = 0;
      row.forEach(function (c) {
        var len = c.value * scale / thick;
        if (w >= h) out.push({ item: c, x: x, y: y + off, w: thick, h: len }); else out.push({ item: c, x: x + off, y: y, w: len, h: thick });
        off += len;
      });
      if (w >= h) { x += thick; w -= thick; } else { y += thick; h -= thick; }
    }
    return out;
  }

  function Pathways() { this.path = []; this.sets = { findings: true, gwas: true, focus: true }; }

  Pathways.prototype.load = function (doc, genes) {
    var P = doc.pathways, self = this, children = {}, universe = new Set();
    Object.keys(P).forEach(function (id) {
      P[id].id = id;
      P[id].symbols = P[id].genes.map(function (e) { var g = genes && genes.byId.get(e); return g ? g.name : null; }).filter(Boolean);
      P[id].symbols.forEach(function (s) { universe.add(s); });
      P[id].parents.forEach(function (pa) { (children[pa] = children[pa] || []).push(id); });
    });
    // a pathway's gene count includes its descendants' genes (Reactome annotates at the lowest level)
    var memo = {};
    var all = function (id) {
      if (memo[id]) return memo[id];
      var s = new Set(P[id].symbols);
      memo[id] = s;
      (children[id] || []).forEach(function (c) { all(c).forEach(function (g) { s.add(g); }); });
      return s;
    };
    Object.keys(P).forEach(function (id) { P[id].all = all(id); });
    this.P = P; this.children = children; this.universe = universe;
    this.top = Object.keys(P).filter(function (id) { return !P[id].parents.length; });
    this.source = doc.source;
    return this;
  };

  // The gene set, from the switches.
  Pathways.prototype.geneSet = function () {
    var view = G.app.view, set = new Set(), s = this.sets;
    if (s.findings) (view.findings || []).forEach(function (f) { String(f.gene).split(/[;,]/).forEach(function (g) { if (g) set.add(g); }); });
    if (s.gwas && G.app.gwasGenes) G.app.gwasGenes.forEach(function (g) { set.add(g); });
    if (s.focus && view.focusGene) set.add(view.focusGene);
    return set;
  };

  Pathways.prototype.score = function () {
    var set = this.geneSet(), U = this.universe, P = this.P, inU = [];
    set.forEach(function (g) { if (U.has(g)) inU.push(g); });
    var N = U.size, n = inU.length, key = inU.slice().sort().join(',');
    if (key === this.scoredFor) return;
    Object.keys(P).forEach(function (id) {
      var p = P[id], k = 0;
      inU.forEach(function (g) { if (p.all.has(g)) k++; });
      p.hits = k; p.pval = k ? hyperTail(N, p.all.size, n, k) : 1;
    });
    this.scoredFor = key; this.nSet = set.size; this.nInU = n;
    this.ranked = Object.keys(P).map(function (id) { return P[id]; }).filter(function (p) { return p.hits; })
      .sort(function (a, b) { return a.pval - b.pval || b.hits - a.hits; });
  };

  Pathways.prototype.draw = function (g) {
    var ctx = g.context, view = G.app.view, self = this;
    if (!this.P) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText('Reactome not loaded. Run python3 tools/fetch_annotations.py', g.cX, g.cY); return; }
    this.score();
    var P = this.P, cur = this.path.length ? this.path[this.path.length - 1] : null;
    var ids = cur ? (this.children[cur] || []) : this.top;
    var listW = Math.min(380, g.cW * 0.3), R = { x: 20, y: 180, w: g.cW - listW - 50, h: g.cH - 250 };
    var boxes = squarify(ids.map(function (id) { return { id: id, value: P[id].all.size }; }), R);
    var over = null;
    boxes.forEach(function (b) {
      var p = P[b.item.id], f = p.hits ? Math.min(1, -Math.log10(p.pval) / 6) : 0;
      ctx.fillStyle = p.hits ? 'rgb(' + Math.round(60 + 190 * f) + ',' + Math.round(60 + 60 * f) + ',' + Math.round(90 - 40 * f) + ')' : '#2b2f3a';
      ctx.fillRect(b.x + 1, b.y + 1, Math.max(0, b.w - 2), Math.max(0, b.h - 2));
      var inside = g.mX > b.x && g.mX < b.x + b.w && g.mY > b.y && g.mY < b.y + b.h;
      if (inside) { over = b; ctx.strokeStyle = 'white'; ctx.lineWidth = 2; ctx.strokeRect(b.x + 1, b.y + 1, b.w - 2, b.h - 2); }
      if (b.w > 50 && b.h > 16) {
        g.setText('rgba(255,255,255,0.9)', Math.min(13, Math.max(9, b.h / 5)), 'Helvetica, Arial, sans-serif', 'left', 'top');
        var name = p.name, maxW = b.w - 8;
        while (name.length > 4 && g.getTextW(name) > maxW) name = name.slice(0, -2);
        if (name !== p.name) name = name.slice(0, -1) + '..';
        g.fText(name, b.x + 5, b.y + 4);
        if (b.h > 34) { g.setText('rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'left', 'top'); g.fText(p.all.size + ' genes' + (p.hits ? ', ' + p.hits + ' in set' : ''), b.x + 5, b.y + 20); }
      }
    });
    // breadcrumb and title
    g.setText('white', 14, 'Helvetica, Arial, sans-serif', 'left', 'top');
    var crumb = ['Reactome'].concat(this.path.map(function (id) { return P[id].name; })).join('  >  ');
    g.fText(crumb.length > 140 ? '...' + crumb.slice(-137) : crumb, 20, 156);
    g.setText('rgba(255,255,255,0.5)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText('gene set: ' + this.nSet + ' genes (' + this.nInU + ' in Reactome). Colour: enrichment, -log10 p. Click to open, click the title to go back.', 20, R.y + R.h + 8);
    if (g.MOUSE_UP_FAST && g.mY > 152 && g.mY < 174 && g.mX < R.x + R.w) this.path.pop();

    // ranked list on the right
    var lx = g.cW - listW - 14, ly = 200, genes = this.geneSet();
    g.setText('white', 12, 'Helvetica, Arial, sans-serif', 'left', 'top'); g.fText(cur && !(this.children[cur] || []).length ? 'Genes of this pathway' : 'Most enriched pathways', lx, ly - 22);
    var rows = cur && !(this.children[cur] || []).length ? Array.from(P[cur].all).sort(function (a, b) { return genes.has(b) - genes.has(a) || (a < b ? -1 : 1); }).slice(0, 60)
      : (this.ranked || []).slice(0, 30);
    rows.forEach(function (r, i) {
      var y = ly + i * 17, isGene = typeof r === 'string';
      var hot = g.mX > lx && g.mX < lx + listW && g.mY >= y && g.mY < y + 17;
      g.setText(isGene ? (genes.has(r) ? 'rgb(255,120,100)' : 'rgba(255,255,255,0.6)') : 'rgba(255,255,255,0.8)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
      var txt = isGene ? r : (r.pval < 1e-3 ? r.pval.toExponential(0) : r.pval.toFixed(3)) + '  ' + r.hits + '/' + r.all.size + '  ' + r.name;
      g.fText(txt.length > 62 ? txt.slice(0, 60) + '..' : txt, lx, y);
      if (hot) {
        g.setCursor('pointer');
        if (g.MOUSE_UP_FAST) {
          if (isGene) { G.app.focusGene(r, false); G.app.setMode('gene'); }
          else { // open the pathway: path from a top-level ancestor down to it
            var chain = [r.id], up = r.parents[0];
            while (up) { chain.unshift(up); up = P[up].parents[0]; }
            self.path = chain;
          }
        }
      }
    });
    if (over) {
      var p = P[over.item.id], inSet = Array.from(p.all).filter(function (x) { return genes.has(x); });
      view.drawTooltip(g, [p.name + '  (' + p.id + ')', p.all.size + ' genes; ' + p.hits + ' from the set' + (p.hits ? ', p = ' + p.pval.toExponential(1) : ''),
        inSet.length ? 'set genes: ' + inSet.slice(0, 12).join(', ') : '', (this.children[p.id] || []).length ? 'click to open ' + this.children[p.id].length + ' sub-pathways' : 'click to list its genes'].filter(Boolean));
      g.setCursor('pointer');
      if (g.MOUSE_UP_FAST && !(g.mY > 152 && g.mY < 174)) this.path.push(p.id);
    }
  };

  G.Pathways = Pathways; G.hyperTail = hyperTail; G.squarify = squarify;
})(globalThis.G = globalThis.G || {});
