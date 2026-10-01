/*
 * Hi-C view: chromatin contacts around the region on screen, from an ENCODE
 * .hic contact map (GRCh38), read by range requests with hic-straw
 * (igvteam, MIT) so only the needed blocks of a many-GB file are fetched.
 * The contacts come from the chosen cell type, not from the sample: they
 * give the 3D context (contact domains) for the enhancer links and SVs.
 *
 * Drawn as the usual triangle: genome along the bottom, contact distance
 * upward. Enhancer-gene links sit at the cell of (element, TSS) as white
 * rings; SV and read-pair ends as cyan rings. Only the file accession and
 * byte ranges go to ENCODE.
 *
 * hic-straw 2.1.4 is used on purpose: 4.x imports a dependency straight from
 * GitHub, which jsDelivr cannot serve to browsers (404), while 2.1.4 ships a
 * self-contained browser bundle with the same reading API.
 */
(function (G) {
  var LIB = 'https://cdn.jsdelivr.net/npm/hic-straw@2.1.4/dist/hic-straw.min.js';
  var ENCODE = 'https://www.encodeproject.org';
  var MAX_BINS = 260, MAX_SPAN = 12e6;
  var Straw = null;

  function lib() {
    if (Straw) return Promise.resolve(Straw);
    return new Promise(function (resolve, reject) {
      var s = document.createElement('script');
      s.src = LIB;
      s.onload = function () { Straw = window.HicStraw; resolve(Straw); };
      s.onerror = function () { reject(new Error('hic-straw could not be loaded from jsDelivr (offline?)')); };
      document.head.appendChild(s);
    });
  }

  function HicView() { this.status = 'idle'; this.files = null; this.file = null; this.norm = null; }

  // GRCh38 in situ Hi-C contact maps on ENCODE, one per biosample (largest first).
  HicView.prototype.listFiles = async function () {
    if (this.files) return this.files;
    var r = await fetch(ENCODE + '/search/?type=File&file_format=hic&assembly=GRCh38&status=released' +
      '&output_type=mapping+quality+thresholded+contact+matrix&limit=all&format=json' +
      '&field=accession&field=file_size&field=biosample_ontology.term_name&field=dataset', { headers: { Accept: 'application/json' } });
    var d = await r.json(), per = {};
    (d['@graph'] || []).forEach(function (f) {
      var b = f.biosample_ontology && f.biosample_ontology.term_name;
      if (!b) return;
      if (!per[b] || f.file_size > per[b].size) per[b] = { acc: f.accession, size: f.file_size, biosample: b };
    });
    this.files = Object.keys(per).sort().map(function (k) { return per[k]; });
    return this.files;
  };

  HicView.prototype.open = async function (acc) {
    var S = await lib();
    this.file = acc;
    this.straw = new S({ url: ENCODE + '/files/' + acc + '/@@download/' + acc + '.hic' });
    this.meta = await this.straw.getMetaData();
    var norms = [];
    try { norms = await this.straw.getNormalizationOptions(); } catch (e) { norms = ['NONE']; }
    this.norms = norms;
    this.norm = norms.indexOf('KR') >= 0 ? 'KR' : norms.indexOf('SCALE') >= 0 ? 'SCALE' : norms.indexOf('VC_SQRT') >= 0 ? 'VC_SQRT' : 'NONE';
    this.chromNames = {};
    var self = this;
    (this.meta.chromosomes || []).forEach(function (c) { self.chromNames[G.genome.normName(c.name)] = c.name; });
    this.region = null; this.img = null;
  };

  // Fetches contacts for a region and renders them into an offscreen image.
  HicView.prototype.load = async function (region) {
    if (!this.straw) return;
    var chrom = this.chromNames[G.genome.normName(region.chrom)];
    if (!chrom) throw new Error(region.chrom + ' is not in this .hic file');
    var span = Math.min(MAX_SPAN, region.end - region.start), mid = (region.start + region.end) / 2;
    var a = Math.max(0, Math.round(mid - span / 2)), b = Math.round(mid + span / 2);
    var res = (this.meta.resolutions || []).slice().sort(function (x, y) { return x - y; });
    var bin = res.find(function (r) { return (b - a) / r <= MAX_BINS; }) || res[res.length - 1];
    a = Math.floor(a / bin) * bin; b = Math.ceil(b / bin) * bin;
    this.status = 'loading contacts at ' + G.fmtBp(bin) + ' resolution';
    var recs = await this.straw.getContactRecords(this.norm, { chr: chrom, start: a, end: b }, { chr: chrom, start: a, end: b }, 'BP', bin);
    var n = Math.round((b - a) / bin), m = new Float32Array(n * n), vals = [];
    recs.forEach(function (r) {
      var i = Math.round((r.bin1 * bin - a) / bin), j = Math.round((r.bin2 * bin - a) / bin);
      if (i < 0 || j < 0 || i >= n || j >= n || !isFinite(r.counts)) return;
      m[i * n + j] = m[j * n + i] = r.counts; vals.push(r.counts);
    });
    vals.sort(function (x, y) { return x - y; });
    var hi = vals.length ? vals[Math.floor(vals.length * 0.97)] : 1;
    // image: cell (i, j) with i <= j at row (j - i), column (i + j) / 2, in a 2n x n grid
    var cv = document.createElement('canvas'); cv.width = 2 * n; cv.height = n;
    var ctx = cv.getContext('2d'), im = ctx.createImageData(2 * n, n), px = im.data;
    for (var i = 0; i < n; i++) for (var j = i; j < n; j++) {
      var v = m[i * n + j]; if (!v) continue;
      var f = Math.min(1, Math.log1p(v) / Math.log1p(hi));
      var row = n - 1 - (j - i), col = i + j; // two image columns per diagonal step
      [col, col + 1].forEach(function (c) {
        if (c >= 2 * n) return;
        var o = (row * 2 * n + c) * 4;
        px[o] = 255; px[o + 1] = Math.round(255 * (1 - f * 0.85)); px[o + 2] = Math.round(255 * (1 - f)); px[o + 3] = Math.round(40 + 215 * f);
      });
    }
    ctx.putImageData(im, 0, 0);
    this.region = { chrom: region.chrom, a: a, b: b, bin: bin, n: n, m: m, hi: hi };
    this.img = cv;
    this.status = 'ready';
  };

  HicView.prototype.draw = function (g) {
    var ctx = g.context, view = G.app.view, R = this.region;
    var msg = function (t) { g.setText('rgba(255,255,255,0.65)', 14, 'Helvetica, Arial, sans-serif', 'center', 'middle'); g.fText(t, g.cX, g.cY); };
    if (!view.data || view.data.build !== 'GRCh38') return msg('The Hi-C view reads GRCh38 contact maps; open a GRCh38 file.');
    if (this.error) return msg(this.error);
    if (!R || !this.img) return msg(this.status === 'idle' ? 'Pick a contact map in the bar above.' : this.status + '...');
    var left = 60, right = g.cW - 60, W = right - left, H = Math.min(W / 2, g.cH - 260), base = 150 + H;
    var xOf = function (pos) { return left + (pos - R.a) / (R.b - R.a) * W; };
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(this.img, left, base - H, W, H);
    ctx.strokeStyle = 'rgba(255,255,255,0.3)'; ctx.beginPath(); ctx.moveTo(left, base); ctx.lineTo(right, base); ctx.stroke();
    var cellXY = function (p, q) { // position pair -> screen point in the triangle
      var x = (xOf(p) + xOf(q)) / 2, y = base - Math.abs(xOf(q) - xOf(p)) * H / W;
      return [x, y];
    };
    var key = G.genome.normName(R.chrom), over = null, bestD = 7;
    // enhancer links as white rings at (element, TSS)
    var reg = G.app.reg;
    if (reg && reg.links) reg.linksInRange(key, R.a, R.b).forEach(function (l) {
      if (l.self || l.tss < R.a || l.tss > R.b || l.mid < R.a || l.mid > R.b) return; // both ends inside the map
      var p = cellXY(l.mid, l.tss), lit = !view.focusGene || l.gene === view.focusGene;
      ctx.strokeStyle = lit ? 'rgba(255,255,255,0.9)' : 'rgba(255,255,255,0.2)'; ctx.lineWidth = l.agree ? 2 : 1;
      ctx.beginPath(); ctx.arc(p[0], p[1], 3, 0, Math.PI * 2); ctx.stroke();
      var dd = Math.hypot(g.mX - p[0], g.mY - p[1]);
      if (dd < bestD) { bestD = dd; over = view.linkLines(l); }
    });
    // SV and read-pair ends
    (view.data.arcs || []).forEach(function (a) {
      if (G.genome.normName(a.c0) !== key || G.genome.normName(a.c1) !== key || a.p0 < R.a || a.p1 > R.b) return;
      var p = cellXY(a.p0, a.p1);
      ctx.strokeStyle = 'rgba(120,255,220,0.9)'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(p[0], p[1], 5, 0, Math.PI * 2); ctx.stroke();
      if (Math.hypot(g.mX - p[0], g.mY - p[1]) < bestD) over = [a.label];
    });
    // genes and findings under the axis
    var genes = view.genes, placed = [];
    if (genes) genes.inRange(key, R.a, R.b).forEach(function (gn) {
      if (gn.type !== 'protein_coding' && gn.name !== view.focusGene) return;
      var x0 = xOf(gn.start), x1 = xOf(gn.end), on = gn.name === view.focusGene;
      ctx.fillStyle = on ? 'white' : 'rgba(255,255,255,0.5)'; ctx.fillRect(x0, base + 8, Math.max(1, x1 - x0), on ? 4 : 3);
      g.setText(on ? 'white' : 'rgba(255,255,255,0.6)', 10, 'Helvetica, Arial, sans-serif', 'left', 'top');
      var tw = g.getTextW(gn.name);
      if (on || placed.every(function (p) { return x0 > p + 4; })) { g.fText(gn.name, x0, base + 15); placed.push(x0 + tw); }
    });
    (view.findings || []).forEach(function (f) {
      if (G.genome.normName(f.chrom) !== key || f.pos < R.a || f.pos > R.b) return;
      ctx.fillStyle = 'rgb(255,70,70)'; ctx.beginPath(); ctx.arc(xOf(f.pos), base + 34, 4, 0, Math.PI * 2); ctx.fill();
    });
    // hover a contact cell
    if (!over && g.mY < base && g.mY > base - H && g.mX > left && g.mX < right) {
      var t = (g.mX - left) / W, dy = (base - g.mY) / H; // position of the cell midpoint and its distance
      var p = R.a + (t - dy / 2) * (R.b - R.a), q = R.a + (t + dy / 2) * (R.b - R.a);
      if (p >= R.a && q <= R.b) {
        var i = Math.floor((p - R.a) / R.bin), j = Math.floor((q - R.a) / R.bin), v = R.m[i * R.n + j];
        over = ['contact ' + R.chrom + ':' + G.fmtBp(p) + ' with ' + G.fmtBp(q) + ' (' + G.fmtBp(q - p) + ' apart)', 'score ' + (v || 0).toFixed(2) + ' (' + this.norm + ', ' + G.fmtBp(R.bin) + ' bins)'];
      }
    }
    g.setText('rgba(255,255,255,0.55)', 11, 'Helvetica, Arial, sans-serif', 'left', 'top');
    g.fText(R.chrom + ':' + R.a.toLocaleString() + '-' + R.b.toLocaleString() + '  |  ' + G.fmtBp(R.bin) + ' bins, ' + this.norm + ' normalised  |  white rings: enhancer to TSS links; cyan: SV ends; red: findings', left, base + 52);
    if (over) view.drawTooltip(g, over);
  };

  G.HicView = HicView;
})(globalThis.G = globalThis.G || {});
