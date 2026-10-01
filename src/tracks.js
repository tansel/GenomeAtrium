/*
 * Tracks view: a linked track browser drawn by GenomeSpy (genomespy.app,
 * MIT, WebGL), fed with the data this page has already parsed.
 *
 * GenomeSpy's own readers take URLs only, so rows are handed over through
 * named datasets and refreshed for what is on screen whenever the view
 * pans or zooms. Tracks, top to bottom:
 *   cytobands, genes, regulatory "squid" (element below, gene TSS above,
 *   after the PISA squid plot), findings, variants (colour: type, or gnomAD
 *   rarity when on), allele fraction from FORMAT AD (the B-allele frequency
 *   plot of copy-number work), and depth ratio log2(DP / median DP).
 * BAM files get depth, discordant read pairs and a Sashimi plot instead.
 *
 * The library (about 1.3 MB) loads from jsDelivr the first time the view
 * is opened. Positions are 1-based here and handed over with offset -1
 * (GenomeSpy's locus scale is 0-based, half-open).
 */
(function (G) {
  var LIB = 'https://cdn.jsdelivr.net/npm/@genome-spy/core@0.89.0';
  var MAX_POINTS = 60000, MAX_LINKS = 4000, MAX_GENES = 3000;
  var STAINS = ['gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100', 'acen', 'gvar', 'stalk'];
  var STAIN_COLORS = ['#2a2a2a', '#555', '#777', '#999', '#bbb', '#c05050', '#6d6da8', '#444'];

  function loadLib() {
    if (window.genomeSpyEmbed) return Promise.resolve();
    return new Promise(function (resolve, reject) {
      var s = document.createElement('script');
      s.src = LIB;
      s.onload = resolve;
      s.onerror = function () { reject(new Error('GenomeSpy could not be loaded from jsDelivr (offline?)')); };
      document.head.appendChild(s);
    });
  }

  function lowerBound(arr, n, v) {
    var lo = 0, hi = n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  // Every track names the shared x scale; a child's own x encoding replaces
  // the root one, so a name only on the root would be lost.
  var X = { chrom: 'chrom', pos: 'pos', type: 'locus', offset: -1, scale: { name: 'gx', zoom: true } };
  function x2(field) { return { chrom: 'chrom', pos: field, offset: -1 }; }

  function track(title, height, body) {
    // width: the 'dark' theme (after vega-themes) defaults views to 300 px wide
    var v = { title: { text: title, orient: 'left', anchor: 'middle', angle: 0, align: 'right', fontSize: 11, color: '#bbb' }, height: height, width: 'container' };
    for (var k in body) v[k] = body[k];
    return v;
  }

  // Spec for a parsed VCF or BAM. Row data comes later through datasets.
  function buildSpec(d, opts) {
    var contigs = d.genome.contigs.map(function (c) { return { name: c.name, size: c.length }; });
    var isVcf = d.format === 'vcf', grch38 = d.build === 'GRCh38';
    var datasets = { bands: [], genes: [], links: [], findings: [], variants: [], sampleGts: [], depthBins: [], pairs: [], junctions: [], gwas: [] };
    var views = [];

    views.push(track('bands', 22, {
      data: { name: 'bands' },
      layer: [
        { mark: { type: 'rect', minOpacity: 1, tooltip: null }, encoding: { x: X, x2: x2('end'),
          color: { field: 'stain', type: 'nominal', scale: { domain: STAINS, range: STAIN_COLORS }, legend: null } } },
        { mark: { type: 'text', size: 9, squeeze: true, tooltip: null }, encoding: { x: X, x2: x2('end'), text: { field: 'name' },
          color: { value: '#ddd' } } }
      ]
    }));
    if (grch38) views.push(track('genes', 42, {
      data: { name: 'genes' },
      layer: [
        { mark: { type: 'rect', y: 0.62, y2: 0.8, minWidth: 1, minOpacity: 1 }, encoding: { x: X, x2: x2('end'),
          color: { field: 'cls', type: 'nominal', scale: { domain: ['coding', 'lncRNA', 'smallRNA', 'pseudogene', 'other'], range: ['#d8d8d8', '#50d2be', '#ff6ebe', '#3a3a3a', '#6a6a6a'] }, legend: null },
          tooltip: [{ field: 'name', title: 'Gene' }, { field: 'type', title: 'Type' }, { field: 'strand', title: 'Strand' }] } },
        { mark: { type: 'text', y: 0.3, size: 10, squeeze: true, tooltip: null }, encoding: { x: X, x2: x2('end'), text: { field: 'name' },
          color: { value: '#ccc' } } }
      ]
    }));
    if (grch38 && opts.hasLinks) views.push(track('enhancers', 130, {
      data: { name: 'links' },
      layer: [
        { mark: { type: 'link', linkShape: 'diagonal', orient: 'vertical', minPickingSize: 3 },
          encoding: { x: X, x2: x2('tss'), y: { datum: 0, type: 'quantitative', scale: { domain: [0, 1] }, axis: null }, y2: { datum: 1 },
            color: { field: 'cls', type: 'nominal', scale: { domain: ['intergenic', 'genic', 'promoter'], range: ['#ffbe46', '#5ac8ff', '#ff64aa'] },
              legend: null },
            opacity: { field: 'strength', type: 'quantitative', scale: { domain: [0, 1], range: [0.15, 0.95] }, legend: null },
            size: { field: 'width', type: 'quantitative', scale: { type: 'identity' }, legend: null },
            tooltip: [{ field: 'gene', title: 'Target gene' }, { field: 'cls', title: 'Element' }, { field: 'models', title: 'Models' },
              { field: 'scores', title: 'Scores' }, { field: 'tissues', title: 'Tissues' }, { field: 'variants', title: 'Sample variants in element' }] } },
        { mark: { type: 'rect', y: 0, y2: 0.07, minWidth: 2, minOpacity: 1, tooltip: null }, encoding: { x: X, x2: x2('end'),
          color: { field: 'hit', type: 'nominal', scale: { domain: ['rare', 'yes', 'no'], range: ['#ff4646', '#ffffff', '#777777'] }, legend: null } } }
      ]
    }));
    if (grch38 && opts.hasGwas) views.push(track('GWAS', 90, {
      data: { name: 'gwas' },
      mark: { type: 'point', size: { expr: 'min(6 * pow(zoomLevel(), 1.1), 90)' }, filled: true, strokeWidth: 0 },
      encoding: { x: X, y: { field: 'lp', type: 'quantitative', scale: { domain: [7, 60], clamp: true }, axis: { title: null, tickCount: 3, labelColor: '#999' } },
        color: { field: 'dose', type: 'nominal', scale: { domain: ['2', '1', '0', 'unknown', 'trait'], range: ['#eb78ff', '#be8cff', '#5a4a70', '#555', '#ffd050'] }, legend: null },
        tooltip: [{ field: 'rsid', title: 'SNP' }, { field: 'trait', title: 'Top trait' }, { field: 'p', title: 'p' }, { field: 'dose', title: 'Risk allele copies here' }, { field: 'gene', title: 'Mapped gene' }] }
    }));
    if (isVcf) views.push(track('findings', 34, {
      data: { name: 'findings' },
      layer: [
        { mark: { type: 'point', size: 90, filled: true, y: 0.35 }, encoding: { x: X,
          color: { field: 'classification', type: 'nominal', scale: { domain: ['Pathogenic', 'Likely pathogenic'], range: ['#ff4646', '#ffa03c'] }, legend: null },
          tooltip: [{ field: 'gene', title: 'Gene' }, { field: 'classification', title: 'Class' }, { field: 'zygosity', title: 'Zygosity' },
            { field: 'af', title: 'gnomAD' }, { field: 'name', title: 'Variant' }] } },
        { mark: { type: 'text', y: 0.8, size: 10, tooltip: null }, encoding: { x: X, text: { field: 'gene' }, color: { value: '#eee' } } }
      ]
    }));
    if (isVcf) {
      views.push(track('variants', 56, {
        data: { name: 'variants' },
        mark: { type: 'point', size: { expr: 'min(4 * pow(zoomLevel(), 1.2), 90)' }, filled: true, strokeWidth: 0 },
        encoding: { x: X,
          y: { field: 'zygosity', type: 'nominal', scale: { domain: ['hom', 'het', 'other'] }, axis: { labelColor: '#999', title: null } },
          color: { field: 'colorBy', type: 'nominal',
            scale: { domain: ['SNV', 'indel', 'SV', 'not in gnomAD', 'ultra-rare', 'rare', 'low frequency', 'common'],
              range: ['#78b4ff', '#ffaa50', '#eb5ac8', '#ff3c3c', '#ff6e3c', '#ffa03c', '#c8c878', '#5a6e8c'] },
            legend: null },
          opacity: { field: 'passAlpha', type: 'quantitative', scale: { type: 'identity' }, legend: null },
          tooltip: [{ field: 'label', title: 'Variant' }, { field: 'zygosity', title: 'Zygosity' }, { field: 'baf', title: 'Allele fraction', format: '.2f' },
            { field: 'dp', title: 'Depth' }, { field: 'gnomad', title: 'gnomAD' }] }
      }));
      if (d.samples.length > 1) views.push(track('samples', Math.min(160, 16 * Math.min(d.samples.length, 16) + 10), {
        data: { name: 'sampleGts' },
        mark: { type: 'point', size: { expr: 'min(4 * pow(zoomLevel(), 1.2), 60)' }, filled: true, strokeWidth: 0 },
        encoding: { x: X, y: { field: 'sample', type: 'nominal', scale: { domain: d.samples.slice(0, 16) }, axis: { labelColor: '#999', title: null } },
          color: { field: 'gt', type: 'nominal', scale: { domain: ['het', 'hom', '0/0', 'missing'], range: ['#be8cff', '#78b4ff', '#555', '#a33'] }, legend: null },
          tooltip: [{ field: 'sample', title: 'Sample' }, { field: 'gt', title: 'Genotype' }, { field: 'label', title: 'Variant' }] }
      }));
      if (opts.hasBaf) {
        views.push(track('allele fraction', 120, {
          data: { name: 'variants' },
          transform: [{ type: 'filter', expr: 'datum.baf !== null' }],
          layer: [
            { mark: { type: 'rule', color: '#555', strokeDash: [3, 3], tooltip: null }, data: { values: [{ y: 0.5 }] }, encoding: { y: { field: 'y', type: 'quantitative' }, x: null } },
            { mark: { type: 'point', size: { expr: 'min(3 * pow(zoomLevel(), 1.2), 60)' }, filled: true, strokeWidth: 0 },
              encoding: { x: X, y: { field: 'baf', type: 'quantitative', scale: { domain: [0, 1] }, axis: { title: null, tickCount: 3, labelColor: '#999' } },
                color: { field: 'zygosity', type: 'nominal', scale: { domain: ['het', 'hom', 'other'], range: ['#be8cff', '#78b4ff', '#888'] }, legend: null },
                opacity: { value: 0.55 } } }
          ]
        }));
        views.push(track('depth ratio', 90, {
          data: { name: 'variants' },
          transform: [{ type: 'filter', expr: 'datum.logr !== null' }],
          layer: [
            { mark: { type: 'rule', color: '#555', strokeDash: [3, 3], tooltip: null }, data: { values: [{ y: 0 }] }, encoding: { y: { field: 'y', type: 'quantitative' }, x: null } },
            { mark: { type: 'point', size: { expr: 'min(3 * pow(zoomLevel(), 1.2), 60)' }, filled: true, strokeWidth: 0 },
              encoding: { x: X, y: { field: 'logr', type: 'quantitative', scale: { domain: [-2, 2], clamp: true }, axis: { title: null, tickCount: 3, labelColor: '#999' } },
                color: { value: '#5ad2be' }, opacity: { value: 0.45 } } }
          ]
        }));
      }
    } else {
      views.push(track('depth', 110, {
        data: { name: 'depthBins' },
        mark: { type: 'rect', minWidth: 0.5, minOpacity: 1 },
        encoding: { x: X, x2: x2('end'), y: { field: 'depth', type: 'quantitative', axis: { title: null, labelColor: '#999' } }, y2: { datum: 0 },
          color: { value: '#5ad2be' }, tooltip: [{ field: 'depth', title: 'Mean depth', format: '.2f' }] }
      }));
      views.push(track('splice junctions (Sashimi)', 110, {
        data: { name: 'junctions' },
        layer: [
          { mark: { type: 'link', linkShape: 'dome', maxChordLength: 100000000 },
            encoding: { x: X, x2: x2('end'), y: { field: 'span', type: 'quantitative', axis: null },
              size: { field: 'reads', type: 'quantitative', scale: { type: 'sqrt', range: [0.5, 4] }, legend: null }, color: { value: '#78dc78' },
              tooltip: [{ field: 'label', title: 'Junction' }, { field: 'reads', title: 'Reads' }] } },
          { mark: { type: 'text', dy: -8, size: 10, tooltip: null }, encoding: { x: { chrom: 'chrom', pos: 'center', type: 'locus', offset: -1, scale: { name: 'gx', zoom: true } },
            y: { field: 'span', type: 'quantitative' }, text: { field: 'reads' }, color: { value: '#cfc' } } }
        ]
      }));
      views.push(track('discordant pairs', 90, {
        data: { name: 'pairs' },
        mark: { type: 'link', linkShape: 'arc', arcHeightFactor: 0.8 },
        encoding: { x: X, x2: { chrom: 'chrom2', pos: 'pos2', offset: -1 }, color: { field: 'kind', type: 'nominal',
          scale: { domain: ['across chromosomes', 'long insert'], range: ['#eb5ac8', '#ff5f5f'] }, legend: null },
          size: { field: 'support', type: 'quantitative', scale: { type: 'sqrt', range: [0.5, 4] }, legend: null },
          tooltip: [{ field: 'label', title: 'Pairs' }] }
      }));
    }

    return {
      theme: 'dark',
      background: '#141414',
      width: 'container',
      genomes: { sample: { contigs: contigs } },
      assembly: 'sample',
      datasets: datasets,
      resolve: { scale: { x: 'shared' }, axis: { x: 'shared' } },
      encoding: { x: { chrom: 'chrom', pos: 'pos', type: 'locus', offset: -1, scale: { name: 'gx', zoom: true }, axis: { orient: 'top' } } },
      spacing: 6,
      vconcat: views
    };
  }

  function Tracks(el) { this.el = el; this.api = null; this.data = null; this.pending = null; }

  Tracks.prototype.open = async function (data, range) {
    await loadLib();
    if (this.data !== data || !this.api) {
      if (this.api) { this.api.finalize(); this.api = null; }
      this.data = data;
      this.prepare(data);
      this.hasLinks = !!(G.app.reg && G.app.reg.links && G.app.reg.links.length);
      this.hasGwas = !!G.app.gwas;
      var spec = buildSpec(data, { hasLinks: this.hasLinks, hasBaf: this.hasBaf, hasGwas: this.hasGwas });
      this.api = await genomeSpyEmbed.embed(this.el, spec, { renderer: 'auto' });
      var self = this;
      this.scale = this.api.getScaleResolutionByName('gx');
      this.scale.addEventListener('domain', function () { self.schedule(); });
      this.pushStatic();
    }
    if (range) await this.zoomTo(range);
    this.refresh();
  };

  Tracks.prototype.close = function () { /* the embed stays alive; it is re-used when the view reopens */ };

  // Median depth and whether the file has allele fractions, for the ratio track.
  Tracks.prototype.prepare = function (d) {
    this.hasBaf = false; this.medianDp = 0;
    if (d.format !== 'vcf') return;
    var sample = [];
    Object.keys(d.variants).forEach(function (k) {
      var c = d.variants[k], step = Math.max(1, Math.floor(c.n / 20000));
      for (var i = 0; i < c.n; i += step) if (c.dp[i]) sample.push(c.dp[i]);
      for (var j = 0; j < Math.min(c.n, 2000); j++) if (c.baf[j]) { this.hasBaf = true; break; }
    }, this);
    sample.sort(function (a, b) { return a - b; });
    this.medianDp = sample.length ? sample[sample.length >> 1] : 0;
  };

  Tracks.prototype.zoomTo = async function (range) {
    if (!this.scale || !range) return;
    try {
      // animated only when the page is visible: hidden tabs pause animation frames
      var anim = document.hidden ? false : { duration: 600 };
      await this.scale.zoomTo([{ chrom: range.chrom, pos: Math.max(0, range.start - 1) }, { chrom: range.chrom, pos: range.end }], anim);
      this.refresh();
    } catch (e) { /* out of range: stay put */ }
  };

  // Visible range as [{chrom, pos}, {chrom, pos}] and span in bp.
  Tracks.prototype.visible = function () {
    var dom = this.scale ? this.scale.getComplexDomain() : null;
    return dom && dom.length === 2 ? dom : null;
  };

  Tracks.prototype.schedule = function () {
    var self = this;
    clearTimeout(this.pending);
    this.pending = setTimeout(function () { self.refresh(); }, 120);
  };

  Tracks.prototype.setData = function (name, rows) { this.api.datasets.set(name, rows); };

  Tracks.prototype.pushStatic = function () {
    var d = this.data, view = G.app.view;
    var bands = [], cyto = view.cytobands || {};
    d.genome.contigs.forEach(function (c) {
      (cyto[c.key] || []).forEach(function (b) { bands.push({ chrom: c.name, pos: b[0] + 1, end: b[1], name: b[2], stain: b[3] }); });
    });
    this.setData('bands', bands);
    if (d.format === 'vcf') {
      var gn = G.app.gnomad;
      this.setData('findings', (view.findings || []).filter(function (f) { return d.genome.get(f.chrom); }).map(function (f) {
        var r = gn && gn.enabled ? gn.lookup(f.chrom, f.pos, f.ref, f.alt) : undefined;
        return { chrom: d.genome.get(f.chrom).name, pos: f.pos, gene: f.gene, classification: f.classification, zygosity: f.zygosity,
          name: f.variant_name, af: r === undefined ? '' : G.gnomad.rarity(r).text };
      }));
    }
    if (d.format === 'bam') {
      var pairs = (d.arcs || []).filter(function (a) { return a.type !== 'junction' && d.genome.get(a.c0) && d.genome.get(a.c1); });
      this.setData('pairs', pairs.slice(0, MAX_LINKS).map(function (a) {
        return { chrom: d.genome.get(a.c0).name, pos: a.p0, chrom2: d.genome.get(a.c1).name, pos2: a.p1, support: a.support, label: a.label,
          kind: a.type === 'pair_inter' ? 'across chromosomes' : 'long insert' };
      }));
      this.setData('junctions', (d.arcs || []).filter(function (a) { return a.type === 'junction'; }).map(function (a) {
        return { chrom: d.genome.get(a.c0).name, pos: a.p0, end: a.p1, center: (a.p0 + a.p1) >> 1, span: a.p1 - a.p0, reads: a.support, label: a.label };
      }));
    }
  };

  // Rows for what is on screen.
  Tracks.prototype.refresh = function () {
    if (!this.api || !this.data) return;
    var d = this.data, dom = this.visible();
    if (!dom) return;
    var parts = [], ga = d.genome.get(dom[0].chrom), gb = d.genome.get(dom[1].chrom);
    if (!ga || !gb) return;
    for (var i = ga.index; i <= gb.index; i++) {
      var c = d.genome.contigs[i];
      parts.push({ c: c, a: i === ga.index ? Math.max(1, dom[0].pos) : 1, b: i === gb.index ? Math.min(c.length, dom[1].pos) : c.length });
    }
    var span = parts.reduce(function (s, p) { return s + (p.b - p.a); }, 0);
    this.span = span;
    if (d.format === 'vcf') this.pushVariants(parts, span);
    else this.pushDepth(parts, span);
    // genes and enhancer links only make sense zoomed in: below 20 Mb
    if (d.build === 'GRCh38') {
      if (this.hasGwas) this.pushGwas(parts);
      if (span <= 20e6) { this.pushGenes(parts, span); this.pushLinks(parts); }
      else { this.setData('genes', []); if (this.hasLinks) this.setData('links', []); }
    }
  };

  Tracks.prototype.pushVariants = function (parts, span) {
    var d = this.data, rows = [], total = 0, gn = G.app.gnomad && G.app.gnomad.enabled ? G.app.gnomad : null, med = this.medianDp;
    var Z = G.vcf.Z, T = G.vcf.T;
    parts.forEach(function (p) { var c = d.variants[p.c.key]; if (c) total += lowerBound(c.pos, c.n, p.b + 1) - lowerBound(c.pos, c.n, p.a); });
    var step = Math.max(1, Math.ceil(total / MAX_POINTS));
    if (gn && span <= 100000) parts.forEach(function (p) { gn.want(p.c.name, p.a, p.b); });
    parts.forEach(function (p) {
      var c = d.variants[p.c.key];
      if (!c) return;
      var i0 = lowerBound(c.pos, c.n, p.a), i1 = lowerBound(c.pos, c.n, p.b + 1);
      for (var i = i0; i < i1; i += step) {
        var z = c.zyg[i], t = c.type[i], baf = c.bafAt(i), al = c.alleles(i);
        var r = gn ? gn.forVariant(p.c.name, c, i) : undefined, rr = r ? G.gnomad.rarity(r) : null;
        rows.push({ chrom: p.c.name, pos: c.pos[i],
          zygosity: z === Z.HET ? 'het' : z === Z.HOM ? 'hom' : 'other',
          colorBy: rr ? rr.text.split(',')[0].replace('not in gnomAD', 'not in gnomAD') : t <= T.MNV ? 'SNV' : t <= T.COMPLEX ? 'indel' : 'SV',
          passAlpha: c.pass[i] ? 0.9 : 0.3,
          baf: baf, dp: c.dp[i] || null, logr: c.dp[i] && med ? Math.log2(c.dp[i] / med) : null,
          label: al ? al.ref.slice(0, 10) + '>' + al.alts.join(',').slice(0, 14) : G.vcf.TYPE_NAMES[t],
          gnomad: rr ? rr.text : '' });
      }
    });
    this.sampled = step > 1 ? step : 0;
    this.setData('variants', rows);
    var ns = Math.min(d.samples.length, 16);
    if (ns > 1) { // per-sample genotypes, thinned further so the track stays light
      var srows = [], names = ['no genotype', 'het', 'hom', '0/0', 'missing'], sstep = step * Math.max(1, Math.ceil(ns / 4));
      parts.forEach(function (p) {
        var c = d.variants[p.c.key];
        if (!c) return;
        for (var i = lowerBound(c.pos, c.n, p.a); i < c.n && c.pos[i] <= p.b; i += sstep) {
          var al = c.alleles(i), lab = al ? al.ref.slice(0, 8) + '>' + al.alts.join(',').slice(0, 10) : '';
          for (var s2 = 0; s2 < ns; s2++) { var z = c.zygOf(i, s2); if (z) srows.push({ chrom: p.c.name, pos: c.pos[i], sample: d.samples[s2], gt: names[z], label: lab }); }
        }
      });
      this.setData('sampleGts', srows);
    }
    if (G.app.onTracksInfo) G.app.onTracksInfo(this);
  };

  Tracks.prototype.pushDepth = function (parts, span) {
    var d = this.data, rows = [], perPart = Math.max(200, Math.floor(3000 / parts.length));
    parts.forEach(function (p) {
      var tr = d.tracks[p.c.key];
      if (!tr) return;
      var lv = tr.depth.levelFor((p.b - p.a) / perPart), bs = lv.binSize;
      for (var b = Math.floor((p.a - 1) / bs); b <= Math.floor((p.b - 1) / bs) && b < lv.data.length; b++) {
        if (lv.data[b] > 0) rows.push({ chrom: p.c.name, pos: b * bs + 1, end: (b + 1) * bs, depth: lv.data[b] });
      }
    });
    this.setData('depthBins', rows);
  };

  Tracks.prototype.pushGenes = function (parts, span) {
    var genes = G.app.view.genes, rows = [];
    if (!genes) return;
    parts.forEach(function (p) {
      genes.inRange(p.c.key, p.a, p.b).forEach(function (g) {
        if (span > 5e6 && g.type !== 'protein_coding') return;
        rows.push({ chrom: p.c.name, pos: g.start, end: g.end, name: g.name, type: g.type, cls: g.cls, strand: g.strand > 0 ? '+' : '-' });
      });
    });
    this.setData('genes', rows.slice(0, MAX_GENES));
  };

  // GWAS SNPs in view, strongest first, capped; the chosen trait's loci in gold.
  Tracks.prototype.pushGwas = function (parts) {
    var gw = G.app.gwas, d = this.data, rows = [], sel = G.app.gwasTraitLoci;
    parts.forEach(function (p) {
      var gc = gw.byContig[p.c.key];
      if (!gc) return;
      gw.inRange(p.c.key, p.a, p.b).forEach(function (i) { rows.push([p, i, gc.lp[i]]); });
    });
    rows.sort(function (a, b) { return b[2] - a[2]; });
    this.setData('gwas', rows.slice(0, 20000).map(function (r) {
      var p = r[0], i = r[1], gc = gw.byContig[p.c.key], dz = gw.dosage(d, p.c.key, i).dosage, top = gw.associations(p.c.key, i)[0];
      return { chrom: p.c.name, pos: gc.pos[i], lp: gc.lp[i], rsid: gc.rsid[i], gene: gc.gene[i], p: top ? (top.p > 0 ? top.p.toExponential(1) : '< 1e-300') : '', trait: top ? top.trait : '',
        dose: sel && sel.set.has(p.c.key + ':' + i) ? 'trait' : dz === null ? 'unknown' : String(dz) };
    }));
  };

  Tracks.prototype.pushLinks = function (parts) {
    var reg = G.app.reg, view = G.app.view;
    if (!reg || !reg.links) return;
    var rows = [];
    parts.forEach(function (p) {
      reg.linksInRange(p.c.key, p.a, p.b).forEach(function (l) {
        if (l.self) return;
        var rare = l.variants && view.linkRare ? view.linkRare(l) : false;
        rows.push({ chrom: p.c.name, pos: l.start, end: l.end, tss: l.tss, gene: l.gene, cls: l.cls, strength: l.score, width: l.agree ? 2.2 : 1,
          models: l.agree ? 'ENCODE-rE2G and ABC' : l.scores.e2g !== undefined ? 'ENCODE-rE2G' : 'ABC',
          scores: [l.scores.e2g !== undefined ? 'rE2G ' + l.scores.e2g.toFixed(2) : '', l.scores.abc !== undefined ? 'ABC ' + l.scores.abc.toFixed(3) : ''].filter(Boolean).join(', '),
          tissues: Object.keys(l.tissues).join(', '), variants: l.variants ? l.variants.length : 0,
          hit: rare ? 'rare' : l.variants ? 'yes' : 'no' });
      });
    });
    rows.sort(function (a, b) { return b.strength - a.strength; });
    this.setData('links', rows.slice(0, MAX_LINKS));
  };

  G.Tracks = Tracks;
  G.tracksSpec = buildSpec;
})(globalThis.G = globalThis.G || {});
