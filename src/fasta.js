/*
 * Reference sequence: from a FASTA file you open, or else from Ensembl.
 *
 * FASTA: plain (not compressed) .fa/.fasta, read with File.slice so only the
 * bytes of the region asked for are read. Its .fai index (samtools faidx) is
 * used when opened with it; without one the page builds the index by one
 * pass over the file (seconds for a small FASTA, about a minute for a whole
 * genome). A compressed FASTA cannot be sliced this way and is refused with
 * a note.
 *
 * refseq.get(build, chrom, start, end): the bases of 1-based [start, end],
 * upper case. Order: a loaded FASTA whose build matches and that has the
 * contig; else Ensembl REST (rest.ensembl.org for GRCh38,
 * grch37.rest.ensembl.org for GRCh37), which receives the coordinates of the
 * region only. Results are cached in 10 kb chunks.
 */
(function (G) {
  function Fasta(file, index, build) { this.file = file; this.index = index; this.build = build; this.name = file.name; }

  // .fai: name, length, offset, line bases, line width.
  function parseFai(text) {
    var idx = {};
    text.split(/\r?\n/).forEach(function (l) {
      var f = l.split('\t');
      if (f.length >= 5) idx[G.genome.normName(f[0])] = { name: f[0], length: +f[1], offset: +f[2], lineBases: +f[3], lineWidth: +f[4] };
    });
    return idx;
  }

  // Build the index by one pass over the bytes (no .fai given).
  async function scanIndex(file, onProgress) {
    var reader = file.stream().getReader(), idx = {}, cur = null, pos = 0, inHeader = false, header = [], lineLen = 0, lineStart = true;
    var finish = function () { if (cur && !cur.lineBases) { cur.lineBases = cur.length; cur.lineWidth = cur.length + 1; } };
    for (;;) {
      var r = await reader.read();
      if (r.done) break;
      var b = r.value;
      for (var i = 0; i < b.length; i++, pos++) {
        var c = b[i];
        if (lineStart && c === 62) { finish(); inHeader = true; header = []; lineStart = false; continue; } // '>'
        if (inHeader) {
          if (c === 10) {
            var name = String.fromCharCode.apply(null, header).trim().split(/\s/)[0];
            cur = { name: name, length: 0, offset: pos + 1, lineBases: 0, lineWidth: 0 };
            idx[G.genome.normName(name)] = cur; inHeader = false; lineStart = true; lineLen = 0;
          } else if (header.length < 200) header.push(c);
          continue;
        }
        if (c === 10) { // end of a sequence line: its width sets the line layout (first line)
          if (cur && !cur.lineBases && lineLen) { cur.lineBases = lineLen - (b[i - 1] === 13 ? 1 : 0); cur.lineWidth = lineLen + 1; }
          lineStart = true; lineLen = 0; continue;
        }
        lineStart = false; lineLen++;
        if (cur && c !== 13) cur.length++;
      }
      if (onProgress) onProgress(pos / file.size);
    }
    finish();
    return idx;
  }

  async function openFasta(file, faiFile, onProgress) {
    if (/\.(gz|bgz|bgzf|zip)$/i.test(file.name)) throw new Error(file.name + ' is compressed: open the uncompressed FASTA (with its .fai from samtools faidx) instead.');
    var index = faiFile ? parseFai(await faiFile.text()) : await scanIndex(file, onProgress);
    var contigs = Object.keys(index).map(function (k) { return { name: index[k].name, length: index[k].length }; });
    if (!contigs.length) throw new Error(file.name + ' holds no FASTA records.');
    return new Fasta(file, index, G.genome.detectBuild(contigs));
  }

  // Bases of 1-based [start, end] on chrom, or null when the contig is absent.
  Fasta.prototype.seq = async function (chrom, start, end) {
    var e = this.index[G.genome.normName(chrom)];
    if (!e) return null;
    start = Math.max(1, start); end = Math.min(e.length, end);
    if (end < start) return '';
    var off = function (p) { return e.offset + Math.floor((p - 1) / e.lineBases) * e.lineWidth + (p - 1) % e.lineBases; };
    var text = await this.file.slice(off(start), off(end) + 1).text();
    return text.replace(/[\r\n]/g, '').toUpperCase();
  };

  // ----- the reference source used by the views
  var CHUNK = 10000, cache = new Map();
  var refseq = {
    fasta: null,
    source: function (build) { return this.fasta && (!build || !this.fasta.build || this.fasta.build === build) ? 'fasta' : build === 'GRCh38' || build === 'GRCh37' ? 'ensembl' : null; },
    get: async function (build, chrom, start, end) {
      var src = this.source(build);
      if (!src) return null;
      var c0 = Math.floor((start - 1) / CHUNK), c1 = Math.floor((end - 1) / CHUNK), parts = [];
      for (var c = c0; c <= c1; c++) parts.push(this.chunk(src, build, chrom, c));
      var s = (await Promise.all(parts)).join('');
      var from = start - 1 - c0 * CHUNK;
      return s.slice(from, from + end - start + 1);
    },
    chunk: function (src, build, chrom, c) {
      var key = src + ':' + build + ':' + G.genome.normName(chrom) + ':' + c, self = this;
      if (cache.has(key)) return cache.get(key);
      var a = c * CHUNK + 1, b = (c + 1) * CHUNK, p;
      if (src === 'fasta') p = this.fasta.seq(chrom, a, b).then(function (s) { return s || ''; });
      else {
        var host = build === 'GRCh37' ? 'https://grch37.rest.ensembl.org' : 'https://rest.ensembl.org';
        var name = G.genome.normName(chrom);
        p = fetch(host + '/sequence/region/human/' + name + ':' + a + '..' + b + ':1?content-type=text/plain').then(function (r) {
          if (!r.ok) throw new Error('Ensembl sequence: HTTP ' + r.status);
          return r.text();
        }).then(function (t) { return t.trim().toUpperCase(); });
      }
      p.catch(function () { cache.delete(key); });
      cache.set(key, p);
      return p;
    },
    clear: function () { cache.clear(); }
  };

  G.fasta = { openFasta: openFasta, parseFai: parseFai, scanIndex: scanIndex, Fasta: Fasta };
  G.refseq = refseq;
})(globalThis.G = globalThis.G || {});
