/*
 * Byte and line readers for plain, gzip and BGZF files.
 *
 * BGZF (used by BAM, .vcf.gz, .bcf) is a chain of small gzip members. The
 * browser's DecompressionStream does not reliably read past the first gzip
 * member, so each BGZF block is split out here and inflated as raw deflate.
 * Works in the browser and in Node 18+, which both provide
 * DecompressionStream, Blob and Response.
 */
(function (G) {
  var SLICE = 8 * 1024 * 1024; // bytes read from disk at a time
  var BATCH = 64;              // BGZF blocks inflated in parallel

  function inflateRaw(bytes) {
    var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
    return new Response(stream).arrayBuffer().then(function (b) { return new Uint8Array(b); });
  }

  function sniff(head) {
    if (head.length < 18 || head[0] !== 0x1f || head[1] !== 0x8b) return 'plain';
    var flg = head[3];
    if (!(flg & 4)) return 'gzip';
    var xlen = head[10] | (head[11] << 8);
    for (var i = 12; i + 4 <= 12 + xlen && i + 4 <= head.length;) {
      var slen = head[i + 2] | (head[i + 3] << 8);
      if (head[i] === 66 && head[i + 1] === 67) return 'bgzf';
      i += 4 + slen;
    }
    return 'gzip';
  }

  // Yields decompressed Uint8Array chunks, in file order.
  // opts.onProgress(bytesRead, totalBytes), opts.signal (AbortSignal)
  async function* chunks(blob, opts) {
    opts = opts || {};
    var total = blob.size;
    var head = new Uint8Array(await blob.slice(0, 64).arrayBuffer());
    var kind = sniff(head);
    var progress = opts.onProgress || function () {};

    if (kind === 'plain') {
      for (var p = 0; p < total; p += SLICE) {
        if (opts.signal && opts.signal.aborted) return;
        yield new Uint8Array(await blob.slice(p, p + SLICE).arrayBuffer());
        progress(Math.min(p + SLICE, total), total);
      }
      return;
    }

    if (kind === 'gzip') {
      var reader = blob.stream().pipeThrough(new DecompressionStream('gzip')).getReader();
      var seen = 0;
      for (;;) {
        if (opts.signal && opts.signal.aborted) { reader.cancel(); return; }
        var r = await reader.read();
        if (r.done) return;
        seen += r.value.length;
        progress(Math.min(seen / 4, total), total); // compressed position unknown; rough
        yield r.value;
      }
    }

    // BGZF: walk block headers, inflate batches of blocks in parallel.
    var buf = new Uint8Array(0), off = 0, filePos = 0;
    for (;;) {
      if (opts.signal && opts.signal.aborted) return;
      if (filePos < total) {
        var next = new Uint8Array(await blob.slice(filePos, filePos + SLICE).arrayBuffer());
        filePos += next.length;
        var merged = new Uint8Array(buf.length - off + next.length);
        merged.set(buf.subarray(off), 0);
        merged.set(next, buf.length - off);
        buf = merged; off = 0;
      }
      var jobs = [];
      while (off + 18 <= buf.length) {
        var xlen = buf[off + 10] | (buf[off + 11] << 8);
        var bsize = -1;
        for (var i = off + 12; i < off + 12 + xlen;) {
          var slen = buf[i + 2] | (buf[i + 3] << 8);
          if (buf[i] === 66 && buf[i + 1] === 67) bsize = buf[i + 4] | (buf[i + 5] << 8);
          i += 4 + slen;
        }
        if (bsize < 0) throw new Error('BGZF block without BC field at byte ' + (filePos - buf.length + off));
        var blockLen = bsize + 1;
        if (off + blockLen > buf.length) break;
        var cdata = buf.subarray(off + 12 + xlen, off + blockLen - 8);
        var isize = buf[off + blockLen - 4] | (buf[off + blockLen - 3] << 8) |
                    (buf[off + blockLen - 2] << 16) | (buf[off + blockLen - 1] << 24);
        if (isize !== 0) jobs.push(inflateRaw(cdata.slice()));
        off += blockLen;
        if (jobs.length >= BATCH) {
          var outs = await Promise.all(jobs); jobs = [];
          for (var j = 0; j < outs.length; j++) yield outs[j];
        }
      }
      if (jobs.length) {
        var rest = await Promise.all(jobs);
        for (var k = 0; k < rest.length; k++) yield rest[k];
      }
      progress(filePos - (buf.length - off), total);
      if (filePos >= total) return;
    }
  }

  // Yields text lines without the trailing newline.
  async function* lines(blob, opts) {
    var dec = new TextDecoder();
    var carry = '';
    for await (var chunk of chunks(blob, opts)) {
      var text = carry + dec.decode(chunk, { stream: true });
      var start = 0, nl;
      while ((nl = text.indexOf('\n', start)) !== -1) {
        var end = nl > start && text.charCodeAt(nl - 1) === 13 ? nl - 1 : nl;
        yield text.slice(start, end);
        start = nl + 1;
      }
      carry = text.slice(start);
    }
    carry += dec.decode();
    if (carry.length) yield carry;
  }

  // Pull-style byte reader over chunks(), for binary formats such as BAM.
  function ByteReader(blob, opts) {
    this.it = chunks(blob, opts);
    this.buf = new Uint8Array(0);
    this.off = 0;
    this.done = false;
  }
  // Ensures n bytes are available from this.off. Returns false at end of data.
  ByteReader.prototype.need = async function (n) {
    while (this.buf.length - this.off < n) {
      if (this.done) return false;
      var r = await this.it.next();
      if (r.done) { this.done = true; return this.buf.length - this.off >= n; }
      var merged = new Uint8Array(this.buf.length - this.off + r.value.length);
      merged.set(this.buf.subarray(this.off), 0);
      merged.set(r.value, this.buf.length - this.off);
      this.buf = merged; this.off = 0;
    }
    return true;
  };
  ByteReader.prototype.take = function (n) {
    var out = this.buf.subarray(this.off, this.off + n);
    this.off += n;
    return out;
  };

  // A string cut from a line (split, slice) can keep the whole decoded text
  // block it came from alive in V8, about 4 MB per block. Anything kept after
  // parsing goes through own() first, which makes a flat copy of just its
  // characters. Without this a 30x gVCF held most of its decompressed text in
  // memory and the browser tab ran out of heap.
  function own(s) { return s == null ? s : (' ' + s).slice(1); }

  G.bgzf = { chunks: chunks, lines: lines, ByteReader: ByteReader, sniff: sniff, own: own };
})(globalThis.G = globalThis.G || {});
