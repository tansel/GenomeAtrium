/*
 * Regulatory layer: GENCODE genes plus enhancer-gene links for the tissues
 * picked in the page, from two models ENCODE publishes side by side:
 * ENCODE-rE2G (a classifier) and ABC (activity by contact). Where both
 * predict the same element and gene, the link is marked as agreed.
 *
 * ENCODE-rE2G predicts, per biosample, which candidate element (a DNase
 * peak) regulates which gene. Each thresholded link has an element class:
 *   promoter    the element is a promoter (of this gene, or of another one)
 *   genic       an enhancer inside a gene body
 *   intergenic  a distal enhancer between genes
 * Links from several biosamples (and several donor sets of one biosample)
 * merge on element + gene; a link keeps every set that predicts it.
 *
 * Coordinates: element BED starts and TSS positions are 0-based in the
 * files and converted to 1-based here, like everything else in the page.
 * GRCh38 only (Asclepius D3).
 */
(function (G) {
  var SETS_PER_BIOSAMPLE = 5; // donor sets merged per picked biosample
  var ENCODE = 'https://www.encodeproject.org';
  var CLASSES = ['intergenic', 'genic', 'promoter'];

  function lowerBound(arr, n, v) {
    var lo = 0, hi = n;
    while (lo < hi) { var mid = (lo + hi) >> 1; if (arr[mid] < v) lo = mid + 1; else hi = mid; }
    return lo;
  }

  // ---------------------------------------------------------------- genes

  function Genes() { this.byContig = {}; this.byName = new Map(); this.byId = new Map(); this.meta = {}; this.n = 0; }

  // GENCODE gene_type -> display class.
  var SMALL_RNA = { miRNA: 1, snoRNA: 1, snRNA: 1, misc_RNA: 1, scaRNA: 1, vault_RNA: 1, Y_RNA: 1, rRNA: 1, ribozyme: 1, sRNA: 1, scRNA: 1, Mt_rRNA: 1, Mt_tRNA: 1 };
  function geneClass(type) {
    if (type === 'protein_coding') return 'coding';
    if (type === 'lncRNA') return 'lncRNA';
    if (SMALL_RNA[type]) return 'smallRNA';
    if (/pseudogene/.test(type)) return 'pseudogene';
    return 'other';
  }

  Genes.prototype.load = async function (blob) {
    var own = G.bgzf.own, tmp = {}, self = this;
    for await (var line of G.bgzf.lines(blob)) {
      if (line[0] === '#') { var eq = line.indexOf('='); if (eq > 0) this.meta[line.slice(1, eq)] = own(line.slice(eq + 1)); continue; }
      var f = line.split('\t');
      if (f.length < 7) continue;
      var g = { chrom: own(f[0]), start: +f[1], end: +f[2], strand: f[3] === '-' ? -1 : 1, name: own(f[4]), id: own(f[5]), type: own(f[6]) };
      g.tss = g.strand > 0 ? g.start : g.end;
      g.cls = geneClass(g.type);
      this.byId.set(g.id, g);
      var key = G.genome.normName(g.chrom);
      (tmp[key] = tmp[key] || []).push(g);
      var up = g.name.toUpperCase(), prev = this.byName.get(up);
      if (!prev || (prev.type !== 'protein_coding' && g.type === 'protein_coding')) this.byName.set(up, g);
      this.n++;
    }
    Object.keys(tmp).forEach(function (k) {
      var list = tmp[k].sort(function (a, b) { return a.start - b.start; });
      var maxLen = list.reduce(function (m, g) { return Math.max(m, g.end - g.start); }, 0);
      self.byContig[k] = { list: list, starts: Int32Array.from(list.map(function (g) { return g.start; })), maxLen: maxLen };
    });
    return this;
  };

  // Genes overlapping [a, b] on a contig (keys as normName gives).
  Genes.prototype.inRange = function (key, a, b) {
    var c = this.byContig[key];
    if (!c) return [];
    var i = lowerBound(c.starts, c.list.length, a - c.maxLen), out = [];
    for (; i < c.list.length && c.list[i].start <= b; i++) if (c.list[i].end >= a) out.push(c.list[i]);
    return out;
  };

  Genes.prototype.get = function (name) { return this.byName.get(String(name).toUpperCase()) || null; };

  // lncRNA / protein-coding pairs that overlap on opposite strands (antisense
  // pairs). Computed once. coordination: optional rows from GeneChords
  // (antisense_lncRNA, sense_gene, marks, auc, rho, p) attached by name.
  Genes.prototype.antisensePairs = function (coordination) {
    if (this.pairs && !coordination) return this.pairs;
    var out = [], self = this, coord = {};
    (coordination || []).forEach(function (r) { coord[r.antisense_lncRNA + '|' + r.sense_gene] = r; });
    Object.keys(this.byContig).forEach(function (k) {
      var list = self.byContig[k].list;
      list.forEach(function (lnc) {
        if (lnc.cls !== 'lncRNA') return;
        self.inRange(k, lnc.start, lnc.end).forEach(function (pc) {
          if (pc.cls !== 'coding' || pc.strand === lnc.strand) return;
          var ov = Math.min(lnc.end, pc.end) - Math.max(lnc.start, pc.start) + 1;
          if (ov <= 0) return;
          out.push({ key: k, lnc: lnc, gene: pc, overlap: ov, coord: coord[lnc.name + '|' + pc.name] || null });
        });
      });
    });
    this.pairs = out;
    this.pairsByKey = {};
    out.forEach(function (p) { (self.pairsByKey[p.key] = self.pairsByKey[p.key] || []).push(p); });
    return out;
  };

  // ---------------------------------------------------------------- links

  // Parses one ENCODE-rE2G thresholded file (full ENCODE format or the slim
  // copy from tools/fetch_regulatory.py; columns are found by header name).
  async function parseLinks(blob, label) {
    var cols = null, out = [], intern = new Map();
    function share(s) { var v = intern.get(s); if (v === undefined) { v = G.bgzf.own(s); intern.set(v, v); } return v; }
    for await (var line of G.bgzf.lines(blob)) {
      var f = line.split('\t');
      if (!cols) {
        var h = f.map(function (x) { return x.replace(/^#/, ''); });
        cols = { chr: h.indexOf('chr'), start: h.indexOf('start'), end: h.indexOf('end'), cls: h.indexOf('class'),
          gene: h.indexOf('TargetGene'), tss: h.indexOf('TargetGeneTSS'), self: h.indexOf('isSelfPromoter'), score: h.indexOf('Score') };
        if (cols.chr < 0 || cols.gene < 0 || cols.tss < 0 || cols.score < 0) throw new Error(label + ': not an ENCODE-rE2G links file');
        continue;
      }
      if (f.length <= cols.score) continue;
      out.push({ chrom: share(f[cols.chr]), start: +f[cols.start] + 1, end: +f[cols.end], cls: share(f[cols.cls] || 'intergenic'),
        gene: share(f[cols.gene]), tss: +f[cols.tss] + 1, self: f[cols.self] === 'TRUE', score: +f[cols.score] });
    }
    return out;
  }

  function Regulatory() {
    this.catalog = null; this.genes = null; this.sets = {}; this.picked = []; this.links = null; this.byContig = {};
  }

  Regulatory.prototype.loadCatalog = async function (blob) {
    this.catalog = JSON.parse(await blob.text());
    var groups = new Map();
    this.catalog.sets.forEach(function (s) {
      var g = groups.get(s.biosample);
      if (!g) { g = { biosample: s.biosample, classification: s.classification, organs: s.organs, sets: [] }; groups.set(s.biosample, g); }
      g.sets.push(s);
    });
    this.biosamples = Array.from(groups.values());
    // Group under the most specific organ: the organ term used by the fewest
    // biosamples (ENCODE lists broad terms such as "exocrine gland" first).
    var freq = {};
    this.biosamples.forEach(function (b) { b.organs.forEach(function (o) { freq[o] = (freq[o] || 0) + 1; }); });
    this.biosamples.forEach(function (b) {
      b.organ = b.organs.slice().sort(function (x, y) { return freq[x] - freq[y] || (x < y ? -1 : 1); })[0] || 'other';
    });
    return this;
  };

  // Local copy first (tools/fetch_regulatory.py --prefetch), then ENCODE.
  async function fetchSet(acc) {
    try {
      var r = await fetch('data/regulatory/' + acc + '.tsv.gz');
      if (r.ok) return { blob: await r.blob(), from: 'local' };
    } catch (e) { /* not served locally: fall through */ }
    var r2 = await fetch(ENCODE + '/files/' + acc + '/@@download/' + acc + '.bed.gz');
    if (!r2.ok) throw new Error(acc + ': ENCODE returned ' + r2.status);
    return { blob: await r2.blob(), from: 'ENCODE' };
  }

  var MODELS = { e2g: 'ENCODE-rE2G', abc: 'ABC' };

  // Display strength 0..1 of a merged link. ENCODE-rE2G scores are already
  // probabilities. ABC scores are fractions of regulatory input, mostly
  // 0.02 to 0.3, so they are scaled (sqrt of score / 0.25) for drawing only;
  // the tooltips show the raw scores.
  function strength(m) {
    var e = m.scores.e2g, a = m.scores.abc;
    var sa = a === undefined ? 0 : Math.min(1, Math.sqrt(a / 0.25));
    return Math.max(e === undefined ? 0 : e, sa);
  }

  // Loads the picked biosamples for the picked models ('e2g', 'abc') and
  // merges their links. onProgress(text).
  Regulatory.prototype.pick = async function (biosampleNames, models, onProgress) {
    var self = this, todo = [], token = this.pickSeq = (this.pickSeq || 0) + 1;
    models = models && models.length ? models : ['e2g'];
    this.picked = biosampleNames.slice(); this.models = models.slice();
    biosampleNames.forEach(function (name) {
      var g = self.biosamples.find(function (b) { return b.biosample === name; });
      if (!g) return;
      g.sets.slice(0, SETS_PER_BIOSAMPLE).forEach(function (s) {
        models.forEach(function (mod) {
          var acc = mod === 'abc' ? s.abc : s.file;
          if (acc) todo.push({ file: acc, model: mod, biosample: name, of: g.sets.length });
        });
      });
    });
    for (var i = 0; i < todo.length; i++) {
      var t = todo[i];
      if (this.sets[t.file]) continue;
      if (onProgress) onProgress('loading ' + t.biosample + ' ' + MODELS[t.model] + ' (' + t.file + ', ' + (i + 1) + ' of ' + todo.length + ')');
      var got = await fetchSet(t.file);
      this.sets[t.file] = { links: await parseLinks(got.blob, t.file), biosample: t.biosample, from: got.from };
    }
    if (token !== this.pickSeq) return null; // a newer pick started meanwhile: it merges, this one does not
    this.merge(todo);
    return this;
  };

  // Merges links of the loaded sets on element + gene, across tissues and
  // models. A link keeps each model's best score and each tissue's scores.
  Regulatory.prototype.merge = function (todo) {
    var map = new Map(), self = this, nn = G.genome.normName;
    this.used = todo.map(function (t) { return { file: t.file, model: t.model, biosample: t.biosample, of: t.of, from: self.sets[t.file] && self.sets[t.file].from }; });
    todo.forEach(function (t) {
      var s = self.sets[t.file];
      if (!s) return;
      s.links.forEach(function (l) {
        var key = l.chrom + ':' + l.start + ':' + l.gene;
        var m = map.get(key);
        if (!m) {
          m = { key: nn(l.chrom), chrom: l.chrom, start: l.start, end: l.end, cls: l.cls, gene: l.gene, tss: l.tss, self: l.self,
            scores: {}, sets: { e2g: 0, abc: 0 }, tissues: {} };
          map.set(key, m);
        }
        m.sets[t.model]++;
        var ts = m.tissues[t.biosample] || (m.tissues[t.biosample] = {});
        ts[t.model] = Math.max(ts[t.model] || 0, l.score);
        m.scores[t.model] = Math.max(m.scores[t.model] || 0, l.score);
      });
    });
    map.forEach(function (m) {
      m.score = strength(m);
      m.agree = m.scores.e2g !== undefined && m.scores.abc !== undefined;
    });
    var byContig = {};
    map.forEach(function (m) {
      m.mid = (m.start + m.end) >> 1;
      m.left = Math.min(m.mid, m.tss); m.right = Math.max(m.mid, m.tss);
      (byContig[m.key] = byContig[m.key] || []).push(m);
    });
    Object.keys(byContig).forEach(function (k) {
      var list = byContig[k].sort(function (a, b) { return a.left - b.left; });
      var maxSpan = list.reduce(function (x, m) { return Math.max(x, m.right - m.left); }, 0);
      var elems = list.slice().sort(function (a, b) { return a.start - b.start; });
      var maxEl = elems.reduce(function (x, m) { return Math.max(x, m.end - m.start); }, 0);
      byContig[k] = { links: list, lefts: Int32Array.from(list.map(function (m) { return m.left; })), maxSpan: maxSpan,
        elems: elems, elemStarts: Int32Array.from(elems.map(function (m) { return m.start; })), maxEl: maxEl };
    });
    this.byContig = byContig;
    this.links = Array.from(map.values());
    this.byGene = new Map();
    for (var i = 0; i < this.links.length; i++) {
      var l = this.links[i], list = this.byGene.get(l.gene);
      if (list) list.push(l); else this.byGene.set(l.gene, [l]);
    }
    this.hits = null;
  };

  // Links with any part in [a, b].
  Regulatory.prototype.linksInRange = function (key, a, b) {
    var c = this.byContig[key];
    if (!c) return [];
    var i = lowerBound(c.lefts, c.links.length, a - c.maxSpan), out = [];
    for (; i < c.links.length && c.links[i].left <= b; i++) if (c.links[i].right >= a) out.push(c.links[i]);
    return out;
  };

  // Sample variants inside linked elements. Adds l.variants (indexes into
  // the VCF columns) to each link whose element carries at least one.
  Regulatory.prototype.intersect = function (data) {
    var total = 0, elems = 0, self = this;
    if (!data || data.format !== 'vcf') { this.hits = null; return null; }
    Object.keys(this.byContig).forEach(function (k) {
      var cols = data.variants[k], c = self.byContig[k];
      c.links.forEach(function (l) { l.variants = null; });
      if (!cols || !cols.n) return;
      c.links.forEach(function (l) {
        var i = lowerBound(cols.pos, cols.n, l.start), v = null;
        for (; i < cols.n && cols.pos[i] <= l.end; i++) (v = v || []).push(i);
        if (v) { l.variants = v; elems++; total += v.length; }
      });
    });
    this.hits = { links: elems, variants: total };
    return this.hits;
  };

  G.regulatory = { Genes: Genes, Regulatory: Regulatory, parseLinks: parseLinks, CLASSES: CLASSES, SETS_PER_BIOSAMPLE: SETS_PER_BIOSAMPLE, geneClass: geneClass,
    MODELS: MODELS, strength: strength };
})(globalThis.G = globalThis.G || {});
