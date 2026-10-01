/*
 * gnomAD allele frequencies, looked up on demand from the gnomAD GraphQL API
 * (gnomad.broadinstitute.org/api, dataset gnomad_r4, GRCh38).
 *
 * Privacy: only region queries are used, over fixed 25 kb tiles. What leaves
 * the page is "chromosome, tile start, tile end". The sample's alleles and
 * genotypes are never sent; the tile's gnomAD variants come back and are
 * matched here. Lookups are off until the user switches them on.
 *
 * A variant absent from a loaded tile is "not in gnomAD". Frequencies
 * combine genomes and exomes: (AC genomes + AC exomes) / (AN genomes + AN exomes).
 * The API allows regions up to about 100 kb; one 25 kb tile takes 2 to 4 s.
 */
(function (G) {
  var API = 'https://gnomad.broadinstitute.org/api';
  var TILE = 25000;
  var GAP_MS = 250;      // pause between requests, to stay gentle with a public API
  var RETRY_MS = 60000;  // a failed tile is not asked for again within a minute

  function query(chrom, start, stop) {
    return '{ region(chrom: "' + chrom + '", start: ' + start + ', stop: ' + stop + ', reference_genome: GRCh38) {' +
      ' variants(dataset: gnomad_r4) { pos ref alt genome { ac an af } exome { ac an af } } } }';
  }

  function rarity(r) {
    if (!r) return null;
    if (r.absent) return { cls: 'novel', text: 'not in gnomAD' };
    var af = r.af;
    if (af < 1e-4) return { cls: 'ultra', text: 'ultra-rare, AF ' + fmtAf(af) };
    if (af < 0.01) return { cls: 'rare', text: 'rare, AF ' + fmtAf(af) };
    if (af < 0.05) return { cls: 'low', text: 'low frequency, AF ' + fmtAf(af) };
    return { cls: 'common', text: 'common, AF ' + fmtAf(af) };
  }
  function fmtAf(af) {
    if (af === 0) return '0';
    if (af >= 0.01) return (af * 100).toFixed(af >= 0.1 ? 0 : 1) + '%';
    return af.toExponential(1);
  }

  function Gnomad(fetchImpl) {
    this.enabled = false;
    this.tiles = new Map();    // "chrom:tile" -> 'queued' | 'loading' | 'done' | 'error'
    this.byKey = new Map();    // normalised variant key -> { af, ac, an, genomeAf, exomeAf }
    this.queue = [];           // tile ids, newest first
    this.errorAt = new Map();  // tile id -> time of its last failure
    this.busy = false;
    this.requests = 0; this.failures = 0;
    this.onUpdate = null;
    this.fetch = fetchImpl || (typeof fetch !== 'undefined' ? fetch.bind(globalThis) : null);
  }

  // gnomAD names chromosomes 1..22, X, Y; the mitochondrion is a separate dataset.
  function gchrom(chrom) {
    var n = G.genome.normName(chrom);
    return /^([0-9]{1,2}|X|Y)$/.test(n) ? n : null;
  }

  Gnomad.prototype.tileId = function (chrom, pos) { return gchrom(chrom) + ':' + Math.floor((pos - 1) / TILE); };

  // Frequency record for one allele: a record, { absent: true }, or
  // undefined while its tile has not been loaded.
  Gnomad.prototype.lookup = function (chrom, pos, ref, alt) {
    var c = gchrom(chrom);
    if (!c) return undefined;
    var r = this.byKey.get(G.clinvar.variantKey(c, pos, ref, alt));
    if (r) return r;
    return this.tiles.get(c + ':' + Math.floor((pos - 1) / TILE)) === 'done' ? { absent: true } : undefined;
  };

  // Rarest allele of a sample variant (columns index i), or undefined if unknown yet.
  Gnomad.prototype.forVariant = function (chrom, cols, i) {
    var al = cols.alleles(i);
    if (!al) return undefined;
    var best;
    for (var k = 0; k < al.alts.length; k++) {
      var r = this.lookup(chrom, cols.pos[i], al.ref, al.alts[k]);
      if (r === undefined) return undefined;
      if (!best || (r.absent ? 0 : r.af) < (best.absent ? 0 : best.af)) best = r;
    }
    return best;
  };

  // Asks for the tiles covering [a, b]. Newest requests are served first.
  Gnomad.prototype.want = function (chrom, a, b) {
    var c = gchrom(chrom);
    if (!this.enabled || !c) return;
    for (var t = Math.floor((a - 1) / TILE); t <= Math.floor((b - 1) / TILE); t++) this.wantTile(c + ':' + t);
  };
  Gnomad.prototype.wantTile = function (id) {
    var st = this.tiles.get(id);
    if (st === 'done' || st === 'loading') return;
    if (st === 'error' && Date.now() - (this.errorAt.get(id) || 0) < RETRY_MS) return; // no retry storm
    if (st === 'queued') this.queue.splice(this.queue.indexOf(id), 1);
    this.tiles.set(id, 'queued');
    this.queue.unshift(id);
    this.pump();
  };
  // Bumped whenever tiles arrive, so drawing code can cache per-item lookups.
  Gnomad.prototype.pending = function () { return this.queue.length + (this.busy ? 1 : 0); };

  Gnomad.prototype.pump = async function () {
    if (this.busy || !this.queue.length || !this.enabled) return;
    this.busy = true;
    var id = this.queue.shift(), parts = id.split(':'), chrom = parts[0], t = +parts[1];
    this.tiles.set(id, 'loading');
    try {
      var r = await this.fetch(API, { method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ query: query(chrom, t * TILE + 1, (t + 1) * TILE) }) });
      var body = await r.json();
      if (!r.ok || body.errors) throw new Error((body.errors && body.errors[0].message) || 'HTTP ' + r.status);
      this.store(chrom, body.data.region.variants);
      this.version = (this.version || 0) + 1;
      this.tiles.set(id, 'done');
      this.requests++;
    } catch (err) {
      this.tiles.set(id, 'error');
      this.errorAt.set(id, Date.now());
      this.failures++;
      this.lastError = err.message;
    }
    this.busy = false;
    if (this.onUpdate) this.onUpdate();
    var self = this;
    setTimeout(function () { self.pump(); }, GAP_MS);
  };

  Gnomad.prototype.store = function (chrom, variants) {
    for (var i = 0; i < variants.length; i++) {
      var v = variants[i], g = v.genome, e = v.exome;
      var ac = (g ? g.ac : 0) + (e ? e.ac : 0), an = (g ? g.an : 0) + (e ? e.an : 0);
      this.byKey.set(G.clinvar.variantKey(chrom, v.pos, v.ref, v.alt), {
        af: an ? ac / an : 0, ac: ac, an: an, genomeAf: g ? g.af : null, exomeAf: e ? e.af : null
      });
    }
  };

  G.gnomad = { Gnomad: Gnomad, rarity: rarity, fmtAf: fmtAf, TILE: TILE, query: query };
})(globalThis.G = globalThis.G || {});
