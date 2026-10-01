/*
 * Gene panels from PanelApp (Genomics England, and PanelApp Australia):
 * expert-curated gene lists per disease, each gene rated green (diagnostic
 * grade), amber or red, with its mode of inheritance.
 *
 * For a chosen panel the page asks, per gene, how much of the gene's span is
 * callable in this sample (gVCF reference blocks and calls passing filters),
 * from the binned callable track (bins of about 1 kb, so exons shorter than
 * a bin are judged with their neighbourhood). It also flags biallelic
 * (recessive) genes that fall in a run of homozygosity, and genes holding a
 * finding. Only panel IDs are requested; no sample data leaves the page.
 */
(function (G) {
  var SOURCES = [
    { id: 'ge', name: 'Genomics England', base: 'https://panelapp.genomicsengland.co.uk/api/v1/' },
    { id: 'au', name: 'PanelApp Australia', base: 'https://panelapp-aus.org/api/v1/' }
  ];
  var LEVEL = { 3: 'green', 2: 'amber', 1: 'red', 0: 'red' };

  async function getJson(url) {
    var r = await fetch(url);
    if (!r.ok) throw new Error('PanelApp: HTTP ' + r.status);
    return r.json();
  }

  var listCache = null;
  // Every panel of both sources: [{source, id, name, version, nGenes}].
  function listPanels() {
    if (listCache) return listCache;
    listCache = Promise.all(SOURCES.map(async function (src) {
      var out = [], url = src.base + 'panels/?page=1';
      while (url) {
        var d = await getJson(url);
        d.results.forEach(function (p) { out.push({ source: src.id, sourceName: src.name, id: p.id, name: p.name, version: p.version, nGenes: p.stats ? p.stats.number_of_genes : 0 }); });
        url = d.next;
      }
      return out;
    })).then(function (lists) { return [].concat.apply([], lists); });
    listCache.catch(function () { listCache = null; });
    return listCache;
  }

  function searchPanels(all, q) {
    q = q.trim().toLowerCase();
    if (q.length < 2) return [];
    return all.filter(function (p) { return p.name.toLowerCase().indexOf(q) >= 0; })
      .sort(function (a, b) { return a.name.length - b.name.length; }).slice(0, 25);
  }

  var panelCache = {};
  // One panel's genes: {source, id, name, version, genes: [{symbol, level, moi, biallelic, loc}]}.
  function loadPanel(source, id) {
    var k = source + ':' + id;
    if (panelCache[k]) return panelCache[k];
    var src = SOURCES.find(function (s) { return s.id === source; });
    panelCache[k] = getJson(src.base + 'panels/' + encodeURIComponent(id) + '/').then(function (d) {
      return {
        source: source, sourceName: src.name, id: d.id, name: d.name, version: d.version,
        genes: (d.genes || []).map(function (g) {
          var gd = g.gene_data || {}, e = (gd.ensembl_genes || {}).GRch38 || {}, ver = Object.keys(e)[0], loc = ver ? e[ver].location : null;
          var m = loc && /^(\w+):(\d+)-(\d+)$/.exec(loc);
          return {
            symbol: gd.gene_symbol || g.entity_name, level: LEVEL[+g.confidence_level] || 'red', moi: g.mode_of_inheritance || '',
            biallelic: /BIALLELIC/i.test(g.mode_of_inheritance || '') && !/X-LINKED/i.test(g.mode_of_inheritance || ''), // X-linked text names biallelic for females loc: m ? { chrom: m[1], start: +m[2], end: +m[3] } : null
          };
        })
      };
    });
    panelCache[k].catch(function () { delete panelCache[k]; });
    return panelCache[k];
  }

  // Fraction of [start, end] covered by a 'mean' track holding per-bin fractions.
  function fraction(track, start, end) {
    var a = track.levels[0], bs = track.binSize, b0 = Math.floor((start - 1) / bs), b1 = Math.floor((end - 1) / bs), s = 0, w = 0;
    for (var b = b0; b <= b1; b++) {
      var lo = Math.max(start, b * bs + 1), hi = Math.min(end, (b + 1) * bs), part = (hi - lo + 1) / bs;
      if (part <= 0) continue;
      s += Math.min(1, b < a.length ? a[b] : 0) * part; w += part;
    }
    return w ? s / w : 0;
  }

  // Per gene of a panel: callable and low-depth fractions, ROH and findings.
  // genes: optional G.Genes table (GENCODE spans); else PanelApp's location.
  // Returns null coverage when the file has no callable track worth reading.
  function assess(data, panel, opts) {
    opts = opts || {};
    var table = opts.genes, roh = opts.roh, findings = opts.findings || [], canCover = !!(data && data.isGvcf);
    var hit = {};
    findings.forEach(function (f) { String(f.gene).split(/[;,]/).forEach(function (g) { (hit[g] = hit[g] || []).push(f); }); });
    var rows = panel.genes.map(function (g) {
      var span = null, t = table && table.get(g.symbol);
      if (t) span = { chrom: t.chrom, start: t.start, end: t.end }; else if (g.loc) span = g.loc;
      var row = { symbol: g.symbol, level: g.level, moi: g.moi, biallelic: g.biallelic, span: span, findings: hit[g.symbol] || [], callable: null, lowdp: null, roh: [] };
      var c = span && data && data.genome.get(span.chrom);
      if (!c) return row;
      row.key = c.key;
      var tr = data.tracks && data.tracks[c.key];
      if (canCover && tr) { row.callable = fraction(tr.callable, span.start, span.end); row.lowdp = fraction(tr.lowdp, span.start, span.end); }
      // autosomal runs only: a male's X is single copy and reads as one long run
      if (roh && G.roh.isAutosome(c.name)) row.roh = G.roh.overlaps(roh, c.key, span.start, span.end);
      return row;
    });
    var status = function (r) { return r.callable === null ? 'unknown' : r.callable >= 0.95 ? 'complete' : r.callable >= 0.8 ? 'gaps' : 'poor'; };
    rows.forEach(function (r) { r.status = status(r); });
    var green = rows.filter(function (r) { return r.level === 'green'; });
    return {
      rows: rows, canCover: canCover,
      summary: {
        genes: rows.length, green: green.length,
        complete: green.filter(function (r) { return r.status === 'complete'; }).length,
        gaps: green.filter(function (r) { return r.status === 'gaps' || r.status === 'poor'; }),
        recessiveInRoh: rows.filter(function (r) { return r.biallelic && r.roh.length; }),
        withFindings: rows.filter(function (r) { return r.findings.length; }),
        unplaced: rows.filter(function (r) { return !r.key; }).length
      }
    };
  }

  G.panels = { SOURCES: SOURCES, listPanels: listPanels, searchPanels: searchPanels, loadPanel: loadPanel, assess: assess, fraction: fraction };
})(globalThis.G = globalThis.G || {});
