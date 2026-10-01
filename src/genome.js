/*
 * Genome model shared by the parsers and the view: contigs, name
 * normalisation, build detection and binned tracks with a zoom pyramid.
 *
 * Coordinates: everything stored here is 1-based inclusive, like VCF POS.
 * The BAM parser converts from its 0-based positions before storing.
 */
(function (G) {
  // chr1 lengths that identify the common human builds.
  var BUILDS = { 248956422: 'GRCh38', 249250621: 'GRCh37' };

  // "chr1" and "1" are the same contig; so are "chrM", "M" and "MT".
  function normName(name) {
    var n = String(name).replace(/^chr/i, '');
    if (n === 'M') n = 'MT';
    return n.toUpperCase();
  }

  function isPrimary(name) {
    return /^([0-9]{1,2}|X|Y|MT|W|Z)$/.test(normName(name));
  }

  // Natural order: numbers first by value, then X, Y, MT, then everything else.
  function contigRank(name) {
    var n = normName(name);
    if (/^[0-9]+$/.test(n)) return [0, +n, n];
    var named = { X: 1, Y: 2, W: 3, Z: 4, MT: 5 }[n];
    if (named) return [1, named, n];
    return [2, 0, n];
  }
  function compareContigs(a, b) {
    var ra = contigRank(a), rb = contigRank(b);
    return ra[0] - rb[0] || ra[1] - rb[1] || (ra[2] < rb[2] ? -1 : ra[2] > rb[2] ? 1 : 0);
  }

  function detectBuild(contigs) {
    for (var i = 0; i < contigs.length; i++) {
      if (normName(contigs[i].name) === '1' && BUILDS[contigs[i].length]) return BUILDS[contigs[i].length];
    }
    return null;
  }

  // Bin width so a whole genome fits in about four million bins.
  function chooseBinSize(totalLength) {
    return Math.max(10, Math.ceil(totalLength / 4e6));
  }

  // A per-contig Float32Array track plus coarser copies for fast drawing.
  // kind 'sum' adds bins together (counts); 'mean' averages them (depth).
  function Track(length, binSize, kind) {
    this.binSize = binSize;
    this.kind = kind || 'sum';
    this.levels = [new Float32Array(Math.ceil(length / binSize) + 1)];
  }
  // Grows the base level when data runs past the declared contig length
  // (VCFs without ##contig lines, or lengths that are wrong).
  Track.prototype._fit = function (b) {
    var a = this.levels[0];
    if (b < a.length) return a;
    var grown = new Float32Array(Math.max(b + 1, a.length * 2));
    grown.set(a);
    this.levels[0] = grown;
    return grown;
  };
  Track.prototype.add = function (pos, value) {
    var b = Math.floor((pos - 1) / this.binSize);
    if (b >= 0) this._fit(b)[b] += value;
  };
  // Adds value * (overlap fraction) to every bin touched by [start, end].
  Track.prototype.addSpan = function (start, end, value) {
    var bs = this.binSize, a = this.levels[0];
    var b0 = Math.max(0, Math.floor((start - 1) / bs));
    var b1 = Math.floor((end - 1) / bs);
    a = this._fit(b1);
    for (var b = b0; b <= b1; b++) {
      var lo = Math.max(start, b * bs + 1), hi = Math.min(end, (b + 1) * bs);
      if (hi >= lo) a[b] += value * (hi - lo + 1) / bs;
    }
  };
  Track.prototype.buildPyramid = function () {
    var cur = this.levels[0];
    this.levels.length = 1;
    while (cur.length > 64) {
      var next = new Float32Array(Math.ceil(cur.length / 4));
      for (var i = 0; i < cur.length; i++) next[i >> 2] += cur[i];
      if (this.kind === 'mean') for (var j = 0; j < next.length; j++) next[j] /= 4;
      this.levels.push(next);
      cur = next;
    }
  };
  // Returns { data, binSize } at the coarsest level still finer than bpPerPixel.
  Track.prototype.levelFor = function (bpPerPixel) {
    var li = 0, bs = this.binSize;
    while (li + 1 < this.levels.length && bs * 4 <= bpPerPixel) { li++; bs *= 4; }
    return { data: this.levels[li], binSize: bs };
  };
  Track.prototype.max = function () {
    var a = this.levels[0], m = 0;
    for (var i = 0; i < a.length; i++) if (a[i] > m) m = a[i];
    return m;
  };

  // Genome: ordered contigs with lengths and a lookup by normalised name.
  function Genome(contigs) {
    this.contigs = [];
    this.byName = {};
    var self = this;
    contigs.forEach(function (c) { self.addContig(c.name, c.length); });
  }
  Genome.prototype.addContig = function (name, length) {
    var key = normName(name);
    if (this.byName[key]) {
      if (length > this.byName[key].length) this.byName[key].length = length;
      return this.byName[key];
    }
    var c = { name: name, key: key, length: length || 0, index: this.contigs.length };
    this.contigs.push(c);
    this.byName[key] = c;
    return c;
  };
  Genome.prototype.get = function (name) { return this.byName[normName(name)]; };
  Genome.prototype.totalLength = function () {
    return this.contigs.reduce(function (s, c) { return s + c.length; }, 0);
  };
  // Chooses the contigs to lay out and sorts them naturally. When the file
  // has primary chromosomes, unplaced, alt and decoy contigs are hidden even
  // if they hold data: a GRCh38 gVCF has refblocks on ~2400 of them, which
  // would otherwise draw as slivers. Otherwise any contig with data is kept.
  Genome.prototype.finish = function (hasData) {
    var withData = this.contigs.filter(function (c) { return c.length > 0 && hasData(c); });
    var primary = this.contigs.filter(function (c) { return c.length > 0 && isPrimary(c.name); });
    var keep = primary.length ? primary : withData;
    this.hidden = withData.filter(function (c) { return keep.indexOf(c) < 0; });
    keep.sort(function (a, b) { return compareContigs(a.name, b.name); });
    keep.forEach(function (c, i) { c.index = i; });
    this.contigs = keep;
    this.build = detectBuild(keep);
  };

  G.genome = {
    Genome: Genome, Track: Track, normName: normName, isPrimary: isPrimary,
    compareContigs: compareContigs, detectBuild: detectBuild, chooseBinSize: chooseBinSize
  };
})(globalThis.G = globalThis.G || {});
