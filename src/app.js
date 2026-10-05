/*
 * Page wiring: file open and drop, format detection, progress, the info
 * panel and legend, go-to box and view mode buttons.
 */
(function (G) {
  // The project was called MoebioToBio until 2026-10-02: carry its saved settings
  // (layers, panels, tissues, bookmarks...) over to the new names once.
  try {
    var OLD = 'moebio' + 'tobio.', oldKeys = []; // list first: adding keys while walking storage by index reorders it
    for (var si = 0; si < localStorage.length; si++) if (localStorage.key(si).indexOf(OLD) === 0) oldKeys.push(localStorage.key(si));
    oldKeys.forEach(function (k) { var nk = 'genomeatrium.' + k.slice(OLD.length); if (localStorage.getItem(nk) === null) localStorage.setItem(nk, localStorage.getItem(k)); });
  } catch (e) { /* storage off */ }
  var view = new G.View('#maindiv');
  var $ = function (id) { return document.getElementById(id); };
  var current = null, controller = null;
  G.app = { view: view, clinvar: null };

  // ClinVar P/LP table exported from Asclepius (tools/export_clinvar.py).
  // Fetched when served over http; can also be dropped on the page.
  var clinvarReady = fetch('data/clinvar_grch38.tsv.gz')
    .then(function (r) { if (!r.ok) throw new Error(r.status); return r.blob(); })
    .then(function (b) { return loadClinvar(b); })
    .catch(function () { clinvarNote('ClinVar table not loaded. Run tools/export_clinvar.py, or drop data/clinvar_grch38.tsv.gz here.'); });

  async function loadClinvar(blob) {
    clinvarNote('loading ClinVar table...');
    var cv = await new G.clinvar.ClinVar().load(blob);
    G.app.clinvar = cv;
    view.clinvar = cv;
    clinvarNote('ClinVar P/LP: ' + n(cv.n) + ' GRCh38 sites (' + esc((cv.meta.source || '').replace(/\s*\(.*\)/, '').trim()) + ')');
    return cv;
  }
  function clinvarNote(html) { $('clinvarNote').innerHTML = html; }

  async function detectFormat(file) {
    var name = file.name.toLowerCase();
    if (/\.bam$/.test(name)) return 'bam';
    if (/\.vcf(\.b?gz)?$/.test(name)) return 'vcf';
    var first = (await G.bgzf.chunks(file).next()).value || new Uint8Array(0);
    if (first[0] === 66 && first[1] === 65 && first[2] === 77 && first[3] === 1) return 'bam';
    if (new TextDecoder().decode(first.subarray(0, 16)).startsWith('##fileformat=VCF')) return 'vcf';
    if (/\.(bcf|cram)$/.test(name)) throw new Error(name.split('.').pop().toUpperCase() + ' is not supported yet. Convert with bcftools view / samtools view -b.');
    throw new Error('Could not tell the format. Expected VCF, gVCF or BAM.');
  }

  async function load(file) {
    if (/\.json$/i.test(file.name)) return loadFindings(file);
    if (/bedmethyl|methyl|\.bed(\.b?gz)?$|\.cov(\.b?gz)?$/i.test(file.name)) return loadMethylation(file);
    // a structural-variant VCF (Sniffles, wf_sv) opened after a genome adds its SVs to it
    if (current && current.format === 'vcf' && /(^|[._-])(sv|svs|sniffles|wf_sv)([._-]|$)/i.test(file.name) && /\.vcf(\.b?gz)?$/i.test(file.name)) return loadSv(file);
    if (/clinvar.*\.tsv(\.gz)?$/i.test(file.name)) return loadClinvar(file).then(function () {
      clinvarNote($('clinvarNote').innerHTML + '<br><span class="dim">Reopen the VCF to match against it.</span>');
    });
    if (controller) controller.abort();
    controller = new AbortController();
    var ctl = controller;
    $('stop').style.display = 'inline-block';
    $('summary').innerHTML = '<span class="file">' + esc(file.name) + '</span><br><span class="dim">reading ' + G.fmtBp(file.size).replace('bp', 'B').replace('kb', 'kB').replace('Mb', 'MB') + '...</span>';
    $('legend').innerHTML = '';
    var t0 = performance.now();
    try {
      var fmt = await detectFormat(file);
      if (fmt === 'vcf') await clinvarReady;
      var opts = {
        clinvar: fmt === 'vcf' ? G.app.clinvar : null,
        signal: ctl.signal,
        onProgress: function (done, total) {
          view.setStatus({ fraction: done / total, text: 'reading ' + file.name + '  ' + Math.round(100 * done / total) + '%' });
        }
      };
      var data = await G[fmt].parse(file, opts);
      if (ctl !== controller) return; // a newer file replaced this one
      data.fileName = file.name;
      data.seconds = (performance.now() - t0) / 1000;
      data.overlays = {};
      if (data.clinvarHits) data.overlays.GenomeAtrium = G.clinvar.toFindings(data.clinvarHits,
        { name: 'GenomeAtrium', code: 'src/clinvar.js' }, { file: file.name }, data.build, 'ClinVar P/LP from Asclepius');
      current = data;
      setPeople(data);
      view.reads = null; view.pileupState = {}; $('readsNote').innerHTML = '';
      view.setStatus(null);
      view.history = []; view.histIdx = null; view.selection = null;
      if (view.mode === 'tracks' || view.mode === 'atrium') setMode('arcs');
      view.setData(data);
      loadCytobands(data.build);
      G.app.onHistory();
      refreshFindings();
      if (reg.links) reg.intersect(data);
      refreshRegulatory();
      refreshSamples();
      refreshNc(); gwasNote(); renderGwas();
      refreshRoh(); refreshPanels(); renderMethylation(); renderPrs(); renderPeople(); renderDepression(false);
      if (G.landscape) G.landscape.setData(data);
      if (G.matrix) G.matrix.setData(data);
      G.app.hilbert.setData(data); G.app.circos.setData(data); renderHilbertLayers();
      describe(data);
      setMode(landingView(data)); // the Atrium first, where it can be shown
    } catch (err) {
      console.error(err);
      view.setStatus(null);
      $('summary').innerHTML = '<span class="file">' + esc(file.name) + '</span><br><span class="warn">' + esc(err.message) + '</span>';
    } finally {
      if (ctl === controller) $('stop').style.display = 'none';
    }
  }

  // The first view for a newly loaded genome: ?view=<mode> if given, else the Atrium
  // (the central room every view opens from), else Arcs when there is no genome to
  // place (an unaligned BAM) or no WebGL.
  var VIEWS = ['arcs', 'tracks', 'circos', 'hilbert', 'gene', 'protein', 'hic', 'pathways', 'mito', 'atrium', 'matrix', '3d'];
  function hasWebGL() {
    try { var c = document.createElement('canvas'); return !!(c.getContext('webgl2') || c.getContext('webgl')); } catch (e) { return false; }
  }
  function landingView(d) {
    var asked = new URLSearchParams(location.search).get('view');
    if (asked === 'arena') asked = 'atrium';
    if (asked && VIEWS.indexOf(asked) >= 0) return asked;
    var placeable = d && d.genome && d.genome.contigs.length && !(d.format === 'bam' && !d.stats.aligned);
    return placeable && hasWebGL() ? 'atrium' : 'arcs';
  }

  function similarText() {
    var m = view.windowModel;
    return m ? 'similarity arcs: ' + n(view.simArcs.length) + ' links between ' + G.fmtBp(m.win) + ' windows with alike profiles (hover the line to pick one window)' : 'similarity arcs';
  }

  function pct(a, b) { return (100 * a / b).toFixed(a / b < 0.01 ? 2 : 1) + '%'; }

  // ----- findings: our ClinVar matches and any dropped findings files

  async function loadFindings(file) {
    try {
      var doc = G.clinvar.parseFindings(await file.text());
      if (!current || !current.genome) throw new Error('Open the VCF first, then drop its findings file.');
      if (current.build && doc.build && current.build !== doc.build) throw new Error('Findings are ' + doc.build + ', the open file is ' + current.build + '. Coordinates are not comparable (D3).');
      current.overlays[(doc.producer && doc.producer.name) || file.name] = doc;
      refreshFindings();
    } catch (err) {
      $('findings').innerHTML = '<span class="warn">' + esc(file.name) + ': ' + esc(err.message) + '</span>';
    }
  }

  // Merges every overlay into one list keyed by normalised variant.
  function mergeFindings(overlays) {
    var byKey = new Map();
    Object.keys(overlays).forEach(function (src) {
      overlays[src].findings.forEach(function (f) {
        var key = f.key || G.clinvar.variantKey(f.chrom, f.pos, f.ref, f.alt);
        var m = byKey.get(key);
        if (!m) { m = { key: key, chrom: f.chrom, pos: f.pos, ref: f.ref, alt: f.alt, gene: f.gene, sources: {} }; byKey.set(key, m); }
        m.sources[src] = f;
        ['classification', 'variant_name', 'phenotype', 'zygosity', 'gt', 'severity', 'action'].forEach(function (k) { if (f[k] && !m[k]) m[k] = f[k]; });
        if (f.status === 'reported') m.reported = true;
      });
    });
    var rank = { Critical: 0, High: 1, Medium: 2, Low: 3 };
    return Array.from(byKey.values()).sort(function (a, b) {
      return (rank[a.severity] == null ? 5 : rank[a.severity]) - (rank[b.severity] == null ? 5 : rank[b.severity]) ||
        (a.classification === 'Pathogenic' ? 0 : 1) - (b.classification === 'Pathogenic' ? 0 : 1) || G.genome.compareContigs(a.chrom, b.chrom) || a.pos - b.pos;
    });
  }

  function refreshFindings() {
    var d = current;
    if (!d || !d.overlays) { $('findings').innerHTML = ''; view.setFindings([]); return; }
    var list = mergeFindings(d.overlays), srcs = Object.keys(d.overlays);
    view.setFindings(list);
    if (G.app.refreshRegulatory) G.app.refreshRegulatory();
    if (G.app.gnomadWants) G.app.gnomadWants();
    if (!srcs.length) { $('findings').innerHTML = ''; return; }
    var head = '<div class="fhead">Findings <span class="dim">(' + srcs.map(function (s) { return esc(s) + ' ' + d.overlays[s].findings.filter(function (f) { return f.status === 'reported'; }).length; }).join(', ') +
      ')</span> <button id="exportFindings" title="Save GenomeAtrium matches as a findings file">export</button></div>';
    var body = list.length ? list.map(function (f, i) {
      var badges = srcs.map(function (s) {
        var x = f.sources[s];
        return x ? '<span class="src ' + (x.status === 'reported' ? '' : 'nr') + '" title="' + esc(s + (x.status === 'reported' ? ': reported' : ': seen, not reported: ' + (x.reason || ''))) + '">' + esc(s[0]) + '</span>' : '<span class="src none" title="' + esc(s) + ': absent">' + esc(s[0]) + '</span>';
      }).join('');
      return '<div class="frow" data-i="' + i + '"><span class="dot ' + (f.classification === 'Pathogenic' ? 'p' : 'lp') + '"></span><b>' + esc(f.gene) + '</b> ' +
        (f.severity ? '<span class="sev ' + esc(f.severity) + '">' + esc(f.severity) + '</span> ' : '') +
        '<span class="dim">' + esc(f.zygosity || '') + ' &middot; ' + esc(f.chrom) + ':' + n(f.pos) + '</span> ' + badges + afBadge(f) +
        '<div class="dim fsub">' + esc((f.phenotype || '').split('|')[0]) + '</div></div>';
    }).join('') : '<span class="dim">No ClinVar P/LP allele found in this file.</span>';
    $('findings').innerHTML = head + '<div class="flist">' + body + '</div>';
    Array.prototype.forEach.call(document.querySelectorAll('.frow'), function (el) {
      el.onclick = function () {
        var f = list[+el.dataset.i];
        if (view.mode === 'protein') { G.app.proteinView.gene = String(f.gene).split(/[;,]/)[0]; return; }
        setMode('arcs'); view.goTo(f.chrom, f.pos - 60, f.pos + 60);
      };
      el.onmouseenter = function () { view.highlightFinding = list[+el.dataset.i]; };
      el.onmouseleave = function () { view.highlightFinding = null; };
    });
    var ex = $('exportFindings');
    if (ex) ex.onclick = function () {
      var doc = d.overlays.GenomeAtrium;
      if (!doc) return;
      var a = document.createElement('a');
      a.href = URL.createObjectURL(new Blob([JSON.stringify(doc, null, 1)], { type: 'application/json' }));
      a.download = d.fileName.replace(/\.(vcf|g\.vcf)(\.gz)?$/i, '') + '.genomeatrium.findings.json';
      a.click();
    };
  }

  // ----- gnomAD frequencies (src/gnomad.js): opt-in, tile coordinates only

  var gnomad = new G.gnomad.Gnomad();
  G.app.gnomad = gnomad; view.gnomad = gnomad;
  gnomad.enabled = (function () { try { return localStorage.getItem('genomeatrium.gnomad') === 'on'; } catch (e) { return false; } })();
  $('gnomadOn').checked = gnomad.enabled;
  var gnomadTimer = null;
  gnomad.onUpdate = function () { // batch UI refreshes while tiles stream in
    gnomadNote();
    if (gnomadTimer) return;
    gnomadTimer = setTimeout(function () { gnomadTimer = null; refreshFindings(); refreshRegulatory(); }, 600);
  };
  function gnomadNote() {
    if (!gnomad.enabled) { $('gnomadNote').textContent = 'off. When on, the page asks gnomad.broadinstitute.org for 25 kb regions (coordinates only, never your alleles).'; return; }
    var done = 0; gnomad.tiles.forEach(function (st) { if (st === 'done') done++; });
    $('gnomadNote').textContent = 'on: ' + done + ' tiles loaded' + (gnomad.pending() ? ', ' + gnomad.pending() + ' pending' : '') +
      (gnomad.failures ? ', ' + gnomad.failures + ' failed (' + gnomad.lastError + ')' : '') + '. Sent: 25 kb tile coordinates only.';
  }
  // What to ask for: every finding, and the elements of the lit gene.
  function gnomadWants() {
    if (!gnomad.enabled || !current || current.build !== 'GRCh38') return;
    (view.findings || []).forEach(function (f) { gnomad.want(f.chrom, f.pos, f.pos); });
    var fg = view.focusGene, ls = fg && reg.byGene ? reg.byGene.get(fg) || [] : [];
    ls.forEach(function (l) { if (l.variants) gnomad.want(l.chrom, l.start, l.end); });
    gnomadNote();
  }
  G.app.gnomadWants = gnomadWants;
  $('gnomadOn').onchange = function () {
    gnomad.enabled = $('gnomadOn').checked;
    try { localStorage.setItem('genomeatrium.gnomad', gnomad.enabled ? 'on' : 'off'); } catch (e) { /* not kept */ }
    if (gnomad.enabled) { gnomadWants(); gnomad.pump(); }
    gnomadNote(); refreshFindings();
  };
  gnomadNote();

  // Frequency badge for a finding, with the ACMG BA1 check (AF above 5%).
  function afBadge(f) {
    if (!gnomad.enabled || !current || current.build !== 'GRCh38') return '';
    var r = gnomad.lookup(f.chrom, f.pos, f.ref, f.alt);
    if (r === undefined) return '<span class="af">AF ...</span>';
    var rr = G.gnomad.rarity(r);
    if (!r.absent && r.af > 0.05 && (f.classification === 'Pathogenic' || f.classification === 'Likely pathogenic'))
      return '<span class="af flag" title="Above 5% in gnomAD: ACMG BA1 (stand-alone benign) would apply unless this is a listed exception. Worth checking the classification.">AF ' + G.gnomad.fmtAf(r.af) + ' BA1?</span>';
    return '<span class="af ' + rr.cls + '" title="gnomAD v4 genomes+exomes: ' + rr.text + '">' + (r.absent ? 'not in gnomAD' : 'AF ' + G.gnomad.fmtAf(r.af)) + '</span>';
  }
  G.app.afBadge = afBadge;


  // ----- navigation: cytobands, history buttons, bookmarks, search palette, shared selection

  var bandCache = {};
  function loadCytobands(build) {
    view.cytobands = null;
    if (!build) return;
    (bandCache[build] || (bandCache[build] = fetch('data/reference/cytobands_' + build + '.tsv.gz').then(okBlob).then(async function (b) {
      var lines = [];
      for await (var line of G.bgzf.lines(b)) lines.push(line);
      return lines.join('\n');
    }).catch(function () { return null; }))).then(function (text) { if (text && current && current.build === build) view.setCytobands(text); });
  }

  G.app.onHistory = function () {
    $('navBack').disabled = !(view.histIdx > 0);
    $('navFwd').disabled = !(view.history && view.histIdx < view.history.length - 1);
  };
  $('navBack').onclick = function () { setMode('arcs'); view.historyGo(-1); };
  $('navFwd').onclick = function () { setMode('arcs'); view.historyGo(1); };

  var BOOK_KEY = 'genomeatrium.bookmarks';
  function bookmarks() { try { return JSON.parse(localStorage.getItem(BOOK_KEY) || '[]'); } catch (e) { return []; } }
  $('bookmark').onclick = function () {
    var b = view.visibleSpan();
    if (!b) return;
    b.name = (view.focusGene ? view.focusGene + ' ' : '') + b.chrom + ':' + b.start.toLocaleString() + '-' + b.end.toLocaleString();
    var list = bookmarks().filter(function (x) { return x.name !== b.name; });
    list.unshift(b);
    try { localStorage.setItem(BOOK_KEY, JSON.stringify(list.slice(0, 50))); } catch (e) { /* not kept */ }
    $('bookmark').innerHTML = '&#9733;';
    setTimeout(function () { $('bookmark').innerHTML = '&#9734;'; }, 900);
  };

  // Shared selection: Arcs (alt+drag), Matrix and Landscape (click a window).
  G.app.onSelection = function (sel) {
    if (G.matrix) G.matrix.selection = sel;
    if (G.landscape) G.landscape.selection = sel;
  };
  G.app.select = function (sel) { view.setSelection(sel); };

  // Search palette: one box for genes, regions, findings, bookmarks, tissues and views.
  var palItems = [], palOn = 0;
  function openPalette() { $('palette').hidden = false; $('paletteInput').value = ''; renderPalette(); $('paletteInput').focus(); }
  function closePalette() { $('palette').hidden = true; }
  function paletteCandidates(q) {
    var out = [], ql = q.trim().toLowerCase();
    var m = /^\s*([^:\s]+):([\d.,]+)([kKmM]?)(?:-([\d.,]+)([kKmM]?))?\s*$/.exec(q);
    var num = function (v, u) { var x = parseFloat(v.replace(/,/g, '')); return Math.round(x * ({ k: 1e3, m: 1e6 }[(u || '').toLowerCase()] || 1)); };
    if (m && current && current.genome && current.genome.get(m[1])) {
      var a = num(m[2], m[3]), b = m[4] ? num(m[4], m[5]) : a + 500;
      if (!m[4]) a = Math.max(1, a - 500);
      out.push({ kind: 'region', label: m[1] + ':' + a.toLocaleString() + '-' + b.toLocaleString(), go: function () { setMode('arcs'); view.goTo(m[1], a, b); } });
    }
    if (current && current.genome) current.genome.contigs.forEach(function (c) {
      if (ql && (c.name.toLowerCase() === ql || c.name.toLowerCase() === 'chr' + ql)) out.push({ kind: 'chromosome', label: c.name, go: function () { setMode('arcs'); view.goTo(c.name, 1, c.length); } });
    });
    (view.findings || []).forEach(function (f) {
      if (!ql || String(f.gene).toLowerCase().indexOf(ql) >= 0) out.push({ kind: 'finding', label: f.gene + '  ' + (f.classification || '') + '  ' + f.chrom + ':' + f.pos.toLocaleString(),
        go: function () { setMode('arcs'); view.goTo(f.chrom, f.pos - 60, f.pos + 60); } });
    });
    if (genes && ql.length >= 2) {
      var n = 0;
      genes.byName.forEach(function (g, name) {
        if (n >= 12 || name.toLowerCase().indexOf(ql) !== 0) return;
        n++;
        out.push({ kind: 'gene', label: g.name + '  ' + g.type.replace(/_/g, ' ') + '  ' + g.chrom + ':' + g.start.toLocaleString(), go: function () { focusGene(g.name, true); } });
      });
    }
    bookmarks().forEach(function (b) {
      if (!ql || b.name.toLowerCase().indexOf(ql) >= 0) out.push({ kind: 'bookmark', label: b.name, go: function () { setMode('arcs'); view.goTo(b.chrom, b.start, b.end); } });
    });
    [['arcs', 'Arcs view'], ['tracks', 'Tracks view'], ['circos', 'Circos view'], ['hilbert', 'Hilbert map'], ['gene', 'Gene view'], ['protein', 'Protein view'], ['hic', 'Hi-C view'], ['pathways', 'Pathways view'], ['mito', 'Mito view (mitochondrial genome)'], ['atrium', 'Atrium (VR)'], ['matrix', 'Matrix view'], ['3d', 'Landscape view']].forEach(function (v) {
      if (ql && v[1].toLowerCase().indexOf(ql) >= 0) out.push({ kind: 'view', label: v[1], go: function () { setMode(v[0]); } });
    });
    if (ql && 'tissues'.indexOf(ql) === 0) out.push({ kind: 'tissues', label: 'pick tissues and models', go: function () { $('openPicker').click(); } });
    if (ql && 'reset whole genome'.indexOf(ql) >= 0) out.push({ kind: 'view', label: 'whole genome', go: function () { setMode('arcs'); view.reset(); } });
    return out.slice(0, 40);
  }
  function renderPalette() {
    palItems = paletteCandidates($('paletteInput').value);
    palOn = Math.min(palOn, Math.max(0, palItems.length - 1));
    $('paletteList').innerHTML = palItems.map(function (it, i) {
      return '<div class="pitem' + (i === palOn ? ' on' : '') + '" data-i="' + i + '"><span class="pkind">' + it.kind + '</span>' + esc(it.label) + '</div>';
    }).join('') || '<div class="dim">type a gene, a region such as chr7:117.5M-117.7M, a finding, a bookmark or a view</div>';
    Array.prototype.forEach.call(document.querySelectorAll('.pitem'), function (el) {
      el.onclick = function () { var it = palItems[+el.dataset.i]; closePalette(); it.go(); };
    });
  }
  $('openPalette').onclick = openPalette;
  $('paletteInput').oninput = function () { palOn = 0; renderPalette(); };
  $('paletteInput').onkeydown = function (e) {
    if (e.key === 'ArrowDown') { palOn = Math.min(palItems.length - 1, palOn + 1); renderPalette(); e.preventDefault(); }
    else if (e.key === 'ArrowUp') { palOn = Math.max(0, palOn - 1); renderPalette(); e.preventDefault(); }
    else if (e.key === 'Enter' && palItems[palOn]) { var it = palItems[palOn]; closePalette(); it.go(); }
    else if (e.key === 'Escape') closePalette();
  };
  document.addEventListener('keydown', function (e) {
    var tag = (e.target && e.target.tagName) || '';
    if ((e.key === 'k' && (e.metaKey || e.ctrlKey)) || (e.key === '/' && tag !== 'INPUT')) { e.preventDefault(); openPalette(); }
    else if (e.key === 'Escape' && !$('palette').hidden) closePalette();
  });


  // ----- samples and ancestry

  var aimPanel = null;
  fetch('data/reference/aim_panel_grch38.json').then(function (r) { return r.ok ? r.json() : null; }).then(function (p) { aimPanel = p; if (current) refreshSamples(); }).catch(function () {});
  function refreshSamples() {
    var d = current, html = [];
    if (!d || d.format !== 'vcf') { $('samplesBox').innerHTML = ''; return; }
    if (aimPanel && d.build === 'GRCh38') {
      var a = G.samples.ancestry(d, aimPanel), pops = aimPanel.superpops.slice().sort(function (x, y) { return a.posterior[y] - a.posterior[x]; });
      html.push('<div class="fhead">Ancestry <span class="dim">(illustrative: ' + a.used.length + ' of ' + a.total + ' Asclepius markers read)</span></div>');
      pops.forEach(function (p) {
        var v = a.posterior[p];
        html.push('<div class="abar"><span style="width:34px">' + p + '</span><span class="b" style="width:' + Math.max(1, Math.round(v * 160)) + 'px"></span>' + (v * 100).toFixed(v < 0.01 ? 1 : 0) + '%</div>');
      });
      html.push('<div class="dim" style="font-size:11px">Hardy-Weinberg likelihood over a small marker panel, as in Asclepius; not admixture or PCA.</div>');
    }
    var n = Math.min(d.samples.length, 16);
    if (n > 1) {
      var cc = G.samples.concordance(d), short = function (x) { return esc(x.length > 14 ? x.slice(0, 12) + '..' : x); };
      html.push('<div class="fhead" style="margin-top:6px">Samples <span class="dim">genotype concordance where both are called</span></div><table><tr><td></td>' +
        d.samples.slice(0, n).map(function (x) { return '<td>' + short(x) + '</td>'; }).join('') + '</tr>' +
        d.samples.slice(0, n).map(function (x, i) {
          return '<tr><td>' + short(x) + '</td>' + d.samples.slice(0, n).map(function (y, j) { var r = i === j ? null : cc.rate(i, j); return '<td>' + (r === null ? '' : (r * 100).toFixed(1) + '%') + '</td>'; }).join('') + '</tr>';
        }).join('') + '</table>');
      if (n >= 3) {
        var opts = function (sel) { return d.samples.slice(0, n).map(function (x, i) { return '<option value="' + i + '"' + (i === sel ? ' selected' : '') + '>' + short(x) + '</option>'; }).join(''); };
        html.push('<div class="prow" style="margin-top:4px">trio: child <select id="trioC">' + opts(0) + '</select> father <select id="trioF">' + opts(1) + '</select> mother <select id="trioM">' + opts(2) + '</select> <button id="trioRun">check</button></div><div id="trioOut" class="dim"></div>');
      }
      html.push('<div class="dim" style="font-size:11px">Views follow the first sample; Tracks shows every sample\'s genotypes.</div>');
    }
    $('samplesBox').innerHTML = html.join('');
    if ($('trioRun')) $('trioRun').onclick = function () {
      var t = G.samples.trio(d, +$('trioC').value, +$('trioF').value, +$('trioM').value);
      var rate = t.checked ? (100 * t.errors / t.checked).toFixed(2) : '0';
      $('trioOut').innerHTML = t.checked.toLocaleString() + ' sites called in all three; ' + t.errors.toLocaleString() + ' Mendelian inconsistencies (' + rate + '%); ' +
        t.deNovo.length + ' de novo candidates (child het, both parents 0/0, PASS)' +
        (t.deNovo.length ? ': ' + t.deNovo.slice(0, 12).map(function (x, i) { return '<a href="#" data-dn="' + i + '">' + x.key + ':' + x.pos.toLocaleString() + '</a>'; }).join(', ') : '');
      Array.prototype.forEach.call(document.querySelectorAll('[data-dn]'), function (el) {
        el.onclick = function (e) { e.preventDefault(); var x = t.deNovo[+el.dataset.dn]; setMode('arcs'); view.goTo(d.genome.get(x.key).name, x.pos - 60, x.pos + 60); };
      });
    };
  }


  // ----- GWAS, ncRNA, pathways

  var gwas = null;
  G.app.gwasGenes = new Set();
  Promise.all([fetch('data/annotations/gwas_snps_grch38.tsv.gz').then(okBlob), fetch('data/annotations/gwas_traits.json').then(function (r) { return r.json(); })])
    .then(function (x) { $('gwasNote').textContent = 'loading...'; return new G.Gwas().load(x[0], x[1]); })
    .then(function (gw) { gwas = G.app.gwas = view.gwas = gw; gwasNote(); renderGwas(); })
    .catch(function () { $('gwasNote').innerHTML = 'not loaded: run <code>python3 tools/fetch_annotations.py</code>'; });
  function gwasNote() {
    if (!gwas) return;
    $('gwasNote').textContent = gwas.n.toLocaleString() + ' SNPs at p <= 5e-8, ' + gwas.traits.length.toLocaleString() + ' traits' + (current && current.build && current.build !== 'GRCh38' ? ' (GRCh38 only: not drawn for ' + current.build + ')' : '');
  }
  var gwasTimer = null;
  $('gwasSearch').oninput = function () { clearTimeout(gwasTimer); gwasTimer = setTimeout(renderGwas, 250); };
  // ----- people (people.js): several persons on one genome, each in a colour
  G.app.people = view.people = [];
  function setPeople(d) {
    var ps = [];
    if (d.format === 'vcf') {
      var names = d.samples && d.samples.length ? d.samples : [d.fileName];
      for (var s = 0; s < Math.min(names.length, G.people.MAX); s++) ps.push(new G.people.Person(d, s, names[s] || d.fileName, s));
    }
    G.app.people = view.people = ps;
  }
  // A genome VCF that should join as a person (not an SV, methylation or findings file).
  function isPersonFile(f) {
    return /\.vcf(\.b?gz)?$/i.test(f.name) && !/(^|[._-])(sv|svs|sniffles|wf_sv)([._-]|$)/i.test(f.name);
  }
  // Several files: the first genome is primary, later genome VCFs join as people, the rest load as usual.
  // Index files pair with their data: a .bai with its BAM (reads for the pileup), a .fai with its FASTA.
  function loadMany(fs) {
    if (!fs.length) return;
    var idx = {};
    fs = fs.filter(function (f) { if (/\.(bai|fai)$/i.test(f.name)) { idx[f.name.replace(/\.(bai|fai)$/i, '').toLowerCase()] = f; return false; } return true; });
    var indexFor = function (f) { var k = f.name.toLowerCase(); return idx[k] || idx[k.replace(/\.(bam|fa|fasta|fna)$/, '')] || null; };
    var fasta = fs.filter(function (f) { return /\.(fa|fasta|fna)(\.gz)?$/i.test(f.name); });
    fs = fs.filter(function (f) { return fasta.indexOf(f) < 0; });
    var first = true;
    var after = function () { return fasta.reduce(function (p, f) { return p.then(function () { return loadFasta(f, indexFor(f)); }); }, Promise.resolve()); };
    var bams = fs.filter(function (f) { return /\.bam$/i.test(f.name) && indexFor(f); });
    if (bams.length && current && current.format === 'vcf') { // reads for the loaded genome
      fs = fs.filter(function (f) { return bams.indexOf(f) < 0; });
      bams.forEach(function (f) { attachReads(G.reads.fileSource(f), G.reads.fileSource(indexFor(f)), f.name); });
    }
    fs.sort(function (a, b) { return isPersonFile(b) - isPersonFile(a); }); // genomes before overlays
    fs.reduce(function (p, f) {
      return p.then(function () {
        if (isPersonFile(f) && !first) return addPerson(f);
        if (isPersonFile(f)) first = false;
        var bai = /\.bam$/i.test(f.name) && indexFor(f);
        return load(f).then(function () { if (bai) return attachReads(G.reads.fileSource(f), G.reads.fileSource(bai), f.name); });
      });
    }, Promise.resolve()).then(after);
  }

  // Reads for the pileup (reads.js): a BAM with its index, local or a URL; same build only.
  async function attachReads(bamSrc, baiSrc, name) {
    try {
      var src = await G.reads.open(bamSrc, baiSrc);
      if (current && src.build && current.build && src.build !== current.build) throw new Error(name + ' is ' + src.build + ', the genome is ' + current.build + ': coordinates are never compared across builds.');
      view.reads = src; view.pileupState = {};
      $('readsNote').innerHTML = 'Reads: ' + esc(name) + ' <span class="dim">(' + (src.build || 'build unknown') + '). Zoom to 4 kb or less in Arcs for the pileup.</span>';
    } catch (err) { console.error(err); $('readsNote').innerHTML = '<span class="warn">Reads: ' + esc(err.message) + '</span>'; }
  }
  G.app.attachReads = attachReads;

  // A reference FASTA (fasta.js): bases at base-level zoom and for the pileup's mismatches.
  async function loadFasta(file, fai) {
    try {
      $('readsNote').innerHTML = '<span class="dim">indexing ' + esc(file.name) + '...</span>';
      var fa = await G.fasta.openFasta(file, fai, function (f) { $('readsNote').innerHTML = '<span class="dim">indexing ' + esc(file.name) + ': ' + Math.round(100 * f) + '%</span>'; });
      if (current && fa.build && current.build && fa.build !== current.build) throw new Error(file.name + ' is ' + fa.build + ', the genome is ' + current.build + '.');
      G.refseq.fasta = fa; G.refseq.clear(); view.pileupState = {};
      $('readsNote').innerHTML = 'Reference: ' + esc(file.name) + ' <span class="dim">(' + Object.keys(fa.index).length + ' sequences, ' + (fa.build || 'build unknown') + (fai ? '' : ', indexed in the page') + ')</span>' + (view.reads ? '<br>Reads: ' + esc(view.reads.name) : '');
    } catch (err) { console.error(err); $('readsNote').innerHTML = '<span class="warn">' + esc(err.message) + '</span>'; }
  }
  async function addPerson(file) {
    var base = current;
    if (!base || base.format !== 'vcf') return load(file); // nothing to add to: it becomes the genome
    if (G.app.people.length >= G.people.MAX) { $('peopleBox').innerHTML += '<div class="warn">Up to ' + G.people.MAX + ' people at a time.</div>'; return; }
    view.setStatus({ fraction: 0, text: 'reading ' + file.name + ' (adding a person)' });
    try {
      await clinvarReady;
      var d = await G.vcf.parse(file, { clinvar: G.app.clinvar, onProgress: function (done, total) { view.setStatus({ fraction: done / total, text: 'reading ' + file.name + ' (adding a person) ' + Math.round(100 * done / total) + '%' }); } });
      if (current !== base) return;
      if (d.build !== base.build) throw new Error(file.name + ' is ' + d.build + ', the genome is ' + base.build + ': coordinates are never compared across builds.');
      d.fileName = file.name;
      var names = d.samples && d.samples.length ? d.samples : [file.name];
      for (var s = 0; s < names.length && G.app.people.length < G.people.MAX; s++) G.app.people.push(new G.people.Person(d, s, names[s] || file.name, G.app.people.length));
      view.people = G.app.people;
      view.setStatus(null);
      renderPeople();
    } catch (err) {
      console.error(err); view.setStatus(null);
      $('peopleBox').innerHTML += '<div class="warn">' + esc(err.message) + '</div>';
    }
  }
  G.app.addPerson = addPerson;
  function renderPeople() {
    var ps = G.app.people, box = $('peopleBox');
    if (!ps || ps.length < 2) { box.innerHTML = ps && ps.length === 1 && current && current.format === 'vcf' ? '<div class="dim">Add person: open more VCFs (Add person, or several files at once) to see a family together.</div>' : ''; if (G.app.onExtrasChanged) G.app.onExtrasChanged(); return; }
    var ROLES = ['', 'child', 'mother', 'father', 'sibling', 'other'];
    var html = '<div class="fhead">People <span class="dim">first is primary; each in its colour</span></div>' + ps.map(function (p, i) {
      var fs = p.findings(), r = p.roh();
      return '<div><input type="checkbox" data-pvis="' + i + '"' + (p.visible ? ' checked' : '') + (i === 0 ? ' disabled' : '') + '>' +
        '<span class="sw" style="background:' + p.color + '"></span><b>' + esc(p.name) + '</b> ' +
        '<select data-prole="' + i + '">' + ROLES.map(function (x) { return '<option' + (x === p.role ? ' selected' : '') + ' value="' + x + '">' + (x || 'role') + '</option>'; }).join('') + '</select>' +
        ' <span class="dim">' + fs.length + ' findings' + (r ? ', F<sub>ROH</sub> ' + (100 * r.fRoh).toFixed(2) + '%' : '') + (p.data !== current ? '' : p.sample ? ', sample ' + (p.sample + 1) : '') + '</span>' +
        (i ? ' <a href="#" data-prm="' + i + '">remove</a>' : '') + '</div>';
    }).join('');
    var by = function (role) { return ps.find(function (p) { return p.role === role; }); }, kid = by('child'), mum = by('mother'), dad = by('father');
    if (kid && mum && dad) {
      var sh = G.people.sharing(kid, mum, dad), pc = function (x) { return (100 * x / Math.max(1, sh.n)).toFixed(1) + '%'; };
      html += '<div style="margin-top:4px"><b>Parent check</b>: of ' + n(sh.n) + ' alleles in ' + esc(kid.name) + ', ' + pc(sh.mother) + ' seen only in ' + esc(mum.name) + ', ' + pc(sh.father) + ' only in ' + esc(dad.name) + ', ' + pc(sh.both) + ' in both, ' + pc(sh.neither) + ' in neither.' +
        ' <span class="dim">Similar shares from each parent and few in neither fit a child of both. "Neither" mixes new mutations with calls a parent\'s file lacks' + (kid.data.isGvcf ? '' : ' (plain VCFs list variants only, so a missing site is unknown, not reference)') + '.</span></div>';
    }
    var shared = G.people.sharedRoh(ps.filter(function (p) { return p.visible; }));
    if (shared.length) html += '<div style="margin-top:4px"><b>Runs of homozygosity in more than one person</b>: ' + shared.slice(0, 8).map(function (s) {
      return '<a href="#" data-sroh="' + s.chrom + ':' + s.start + ':' + s.end + '">' + esc(s.chrom) + ':' + G.fmtBp(s.start) + '-' + G.fmtBp(s.end) + '</a> <span class="dim">(' + s.who.map(function (w) { return esc(ps[w].name); }).join(', ') + ')</span>';
    }).join('; ') + '</div>';
    box.innerHTML = html;
    Array.prototype.forEach.call(box.querySelectorAll('[data-pvis]'), function (el) { el.onchange = function () { ps[+el.dataset.pvis].visible = el.checked; renderPeople(); }; });
    Array.prototype.forEach.call(box.querySelectorAll('[data-prole]'), function (el) { el.onchange = function () { ps[+el.dataset.prole].role = el.value; renderPeople(); }; });
    Array.prototype.forEach.call(box.querySelectorAll('[data-prm]'), function (el) {
      el.onclick = function (e) { e.preventDefault(); ps.splice(+el.dataset.prm, 1); ps.forEach(function (p, i) { p.index = i; p.color = G.people.PALETTE[i]; }); renderPeople(); };
    });
    Array.prototype.forEach.call(box.querySelectorAll('[data-sroh]'), function (el) {
      el.onclick = function (e) { e.preventDefault(); var x = el.dataset.sroh.split(':'); G.app.setMode('arcs'); view.goTo(x[0], +x[1], +x[2]); };
    });
    if (G.app.onExtrasChanged) G.app.onExtrasChanged();
  }
  G.app.renderPeople = renderPeople;

  // ----- structural variants from a separate VCF, added to the loaded genome
  async function loadSv(file) {
    var base = current;
    view.setStatus('reading ' + file.name + ' (structural variants, added to ' + base.fileName + ')');
    try {
      var sv = await G.vcf.parse(file, {});
      if (current !== base) return;
      if (sv.build !== base.build) throw new Error(file.name + ' is ' + sv.build + ', the genome is ' + base.build + ': coordinates are never compared across builds.');
      var added = sv.arcs.filter(function (a) { return a.type !== 'junction'; });
      base.arcs = (base.arcs || []).concat(added);
      Object.keys(sv.tracks).forEach(function (k) {
        var to = base.tracks[k], from = sv.tracks[k];
        if (!to || !from) return;
        var a = from.sv.levels[0], b = to.sv._fit(a.length - 1);
        for (var i = 0; i < a.length; i++) b[i] += a[i];
        to.sv.buildPyramid();
      });
      base.svFile = file.name; base.svAdded = added.length;
      view.setStatus(null);
      view.setData(base); describe(base);
      if (G.app.atrium) { // rebuild with the new arches: now if the Atrium is showing, else next time
        if (view.mode === 'atrium' && G.app.atrium.group) G.app.atrium.build(base); else G.app.atrium.built = null;
      }
    } catch (err) {
      console.error(err);
      view.setStatus(null);
      $('summary').innerHTML += '<br><span class="warn">' + esc(err.message) + '</span>';
    }
  }

  // ----- polygenic scores (prs.js): Asclepius's rules; the percentile only from a supplied reference
  var prsTimer = null;
  $('prsSearch').addEventListener('input', function () {
    var q = this.value.trim();
    clearTimeout(prsTimer);
    if (/^PGS\d{6}$/i.test(q)) { $('prsResults').innerHTML = '<a href="#" data-pgs="' + q.toUpperCase() + '">score ' + q.toUpperCase() + '</a>'; bindPgs(); return; }
    if (q.length < 3) { $('prsResults').innerHTML = ''; return; }
    prsTimer = setTimeout(function () {
      $('prsResults').innerHTML = '<span class="dim">searching...</span>';
      G.prs.searchTraits(q).then(function (ts) {
        $('prsResults').innerHTML = ts.length ? ts.slice(0, 8).map(function (t) {
          return '<div>' + esc(t.label) + ': ' + t.scores.slice(0, 10).map(function (id) { return '<a href="#" data-pgs="' + esc(id) + '">' + esc(id) + '</a>'; }).join(' ') + (t.scores.length > 10 ? ' <span class="dim">+' + (t.scores.length - 10) + '</span>' : '') + '</div>';
        }).join('') : '<span class="dim">no trait matches</span>';
        bindPgs();
      }, function (err) { $('prsResults').innerHTML = '<span class="warn">' + esc(err.message) + '</span>'; });
    }, 350);
  });
  function bindPgs() {
    Array.prototype.forEach.call($('prsResults').querySelectorAll('[data-pgs]'), function (el) { el.onclick = function (e) { e.preventDefault(); loadScore(el.dataset.pgs); }; });
  }
  async function loadScore(id) {
    var out = $('prsOut'), base = current, say = function (t) { out.innerHTML = '<div class="dim">' + t + '</div>'; };
    if (!base || base.format !== 'vcf') { say('Open a VCF or gVCF first.'); return; }
    if (base.build !== 'GRCh38') { out.innerHTML = '<div class="warn">Scores are harmonised to GRCh38; this file is ' + esc(base.build || 'of unknown build') + '. Coordinates are never compared across builds.</div>'; return; }
    try {
      say('reading ' + esc(id) + ' from the PGS Catalog...');
      var info = await G.prs.scoreInfo(id);
      if (!info.url) throw new Error(id + ' has no GRCh38 harmonised file');
      say('downloading ' + esc(id) + ' (' + n(info.n) + ' variants)...');
      var sf = await G.prs.parseScoreFile(await (await fetch(info.url)).blob());
      var refNote;
      var local = await fetch('data/pgs/' + id + '.refbase.tsv.gz').catch(function () { return null; });
      if (local && local.ok) { say('reading reference bases from data/pgs...'); sf.ref = await G.prs.refBasesFile(sf, await local.blob()); refNote = 'reference bases from tools/pgs_refbases.py (UCSC hg38)'; }
      else if (base.isGvcf && sf.n <= 20000) {
        sf.ref = await G.prs.refBasesEnsembl(sf, function (f) { say('reference bases from Ensembl for all ' + n(sf.n) + ' score sites: ' + Math.round(100 * f) + '%'); });
        refNote = 'reference bases from Ensembl (asked for every score site)';
      } else refNote = base.isGvcf ? 'no reference bases: run <code>python3 tools/pgs_refbases.py ' + esc(id) + '</code> to score the hom-ref sites' : 'a plain VCF has no reference blocks, so sites without a call are unknown';
      if (current !== base) return;
      var r = G.prs.compute(base, sf);
      r.info = info; r.refNote = refNote;
      G.app.prs = view.prs = base.prs = r;
      renderPrs();
    } catch (err) { console.error(err); out.innerHTML = '<div class="warn">' + esc(err.message) + '</div>'; }
  }
  function renderPrs() {
    var r = G.app.prs, out = $('prsOut');
    if (!r || !current || current.prs !== r) { out.innerHTML = ''; return; }
    var pc = function (x) { return (100 * x).toFixed(1) + '%'; }, why = Object.keys(r.unknownWhy).filter(function (k) { return r.unknownWhy[k]; }).map(function (k) { return n(r.unknownWhy[k]) + ' ' + k; }).join(', ');
    out.innerHTML = '<div class="fhead" style="margin-top:4px">' + esc(r.info.id) + ' ' + esc(r.info.name || '') + ' <span class="dim">' + esc(r.info.trait || '') + (r.info.publication ? '; ' + esc(r.info.publication) : '') + '</span></div>' +
      '<div>Raw score <b>' + r.score.toFixed(4) + '</b> <span class="warn">uncalibrated</span></div>' +
      '<div>' + n(r.called) + ' sites called, ' + n(r.homref) + ' confident hom-ref, ' + n(r.unknown) + ' unknown (' + esc(why) + '): ' + pc(r.coverage) + ' of ' + n(r.n) + ' sites used. Unknown sites are left out, never counted as 0.</div>' +
      '<div class="dim">' + r.refNote + '. At called sites the effect allele is the reference ' + (r.effectRefRate === null ? '-' : pc(r.effectRefRate)) + ' of the time, so "effect = non-reference" cannot be assumed.</div>' +
      '<div>Percentile against your ancestry-matched reference: mean <input id="prsMean" type="number" step="any" style="width:70px"> sd <input id="prsSd" type="number" step="any" style="width:60px"> <span id="prsPct" class="dim">needs both</span></div>' +
      '<div class="dim">A raw score has no meaning on its own. The percentile, and any clinical reading, belong to Asclepius.</div>' +
      '<div style="margin-top:4px">Largest contributions:</div>' + r.top.slice(0, 12).map(function (x, i) {
        return '<div><a href="#" data-prs="' + i + '">' + esc(x.chrom) + ':' + n(x.pos) + '</a> ' + esc(x.effect) + ' x' + x.dose + ' (' + (x.how === 'homref' ? 'hom-ref' : 'called') + ') weight ' + x.weight.toFixed(3) + ' = <b style="color:' + (x.contrib > 0 ? 'rgb(255,110,90)' : 'rgb(110,170,255)') + '">' + (x.contrib > 0 ? '+' : '') + x.contrib.toFixed(3) + '</b></div>';
      }).join('');
    var upd = function () {
      var m = parseFloat($('prsMean').value), sd = parseFloat($('prsSd').value);
      $('prsPct').textContent = isFinite(m) && sd > 0 ? Math.round(G.prs.prsPercentile(r.score, m, sd)) + 'th percentile (only as good as the reference you gave)' : 'needs both';
    };
    $('prsMean').oninput = $('prsSd').oninput = upd;
    Array.prototype.forEach.call(out.querySelectorAll('[data-prs]'), function (el) { el.onclick = function (e) { e.preventDefault(); var x = r.top[+el.dataset.prs]; G.app.setMode('arcs'); view.goTo(x.chrom, x.pos - 50, x.pos + 50); }; });
    describe(current);
  }

  // ----- the depression card: antidepressant pharmacogenomics (Asclepius pgx, CPIC),
  // depression polygenic scores (PGS Catalog) and GWAS loci. Evidence, cited; no
  // metaboliser status (Asclepius D7) and no diagnosis.
  var CPIC_ANTIDEPRESSANTS = [
    ['SSRIs and SNRIs (CYP2D6, CYP2C19, CYP2B6)', 'https://cpicpgx.org/guidelines/cpic-guideline-for-ssri-and-snri-antidepressants/'],
    ['Tricyclic antidepressants (CYP2D6, CYP2C19)', 'https://cpicpgx.org/guidelines/guideline-for-tricyclic-antidepressants-and-cyp2d6-and-cyp2c19/']];
  var depressionScores = null;
  $('depressionOpen').onclick = function () { renderDepression(true); };
  function renderDepression(open) {
    var out = $('depressionOut'), d = current;
    if (!open && !out.innerHTML) return;
    if (!d || d.format !== 'vcf') { out.innerHTML = '<div class="dim">Open a VCF or gVCF first.</div>'; return; }
    if (d.build !== 'GRCh38') { out.innerHTML = '<div class="warn">The card uses GRCh38 positions (CPIC sites, GWAS Catalog, PGS); this file is ' + esc(d.build || 'of unknown build') + '.</div>'; return; }
    var STATE = function (g) {
      return g.state === 'copies' ? '<b>' + g.copies + (g.copies === 1 ? ' copy' : ' copies') + '</b>' : g.state === 'reference' ? '0 copies <span class="dim">(reference, gVCF block)</span>'
        : g.state === 'other' ? '<span class="warn">other allele</span> <span class="dim">(' + esc(g.why) + ')</span>' : '<span class="warn">' + (g.state === 'unknown' ? 'unknown' : g.state === 'lowdp' ? 'low depth' : 'not called') + '</span> <span class="dim">(' + esc(g.why || '') + ')</span>';
    };
    var pgxRows = ['CYP2C19', 'CYP2D6'].map(function (gene) {
      var gd = G.pgx.pgxGuidance(gene), vs = G.pgx.PGX_VARIANTS.filter(function (v) { return v.gene === gene; });
      return '<div style="margin-top:3px"><b>' + gene + '</b> <span class="dim">CPIC level ' + esc(gd.cpic_level) + '</span>' + vs.map(function (v) {
        return '<div>&nbsp; ' + esc(v.star) + ' <span class="dim">' + esc(v.rsid) + ' ' + v.ref + '>' + v.alt + '</span>: ' + STATE(G.pgx.genotypeAt(d, v, 0)) + '</div>';
      }).join('') + '<div class="dim">&nbsp; ' + esc(gd.recommendation) + '</div></div>';
    }).join('');
    var html = '<div><b>Antidepressants</b> <span class="dim">(pharmacogenomics, from Asclepius pgx ' + esc(G.pgx.KB_VERSION) + ')</span></div>' + pgxRows +
      '<div class="dim" style="margin-top:3px">Copies of the variant that defines each star allele. A metaboliser status needs the full diplotype, and CYP2D6 also its copy number (whole-gene deletions such as *5, duplications, CYP2D7 hybrids), which a VCF does not hold; dedicated callers exist (Cyrius, Aldy). None is assigned here (Asclepius D7). CYP2B6 is in the SSRI guideline but not yet in the Asclepius knowledge base.</div>' +
      '<div>CPIC guidelines: ' + CPIC_ANTIDEPRESSANTS.map(function (c) { return '<a href="' + c[1] + '" target="_blank" rel="noopener">' + esc(c[0]) + '</a>'; }).join('; ') + '</div>';
    // polygenic scores
    var r = G.app.prs && current.prs === G.app.prs && /depress/i.test(G.app.prs.info.trait || '') ? G.app.prs : null;
    html += '<div style="margin-top:6px"><b>Polygenic scores</b> <span class="dim">(PGS Catalog, major depressive disorder)</span></div>';
    if (r) html += '<div>' + esc(r.info.id) + ': raw score ' + r.score.toFixed(4) + ' <span class="warn">uncalibrated</span>, ' + (100 * r.coverage).toFixed(1) + '% of ' + n(r.n) + ' sites used. <span class="dim">Without an ancestry-matched reference the number has no meaning on its own; these scores explain only a few percent of the variation in depression.</span></div>';
    html += '<div id="depScores">' + (depressionScores ? depressionScores.map(function (s) { return '<a href="#" data-dpgs="' + esc(s.id) + '">' + esc(s.id) + '</a> <span class="dim">' + esc(s.name || '') + ', ' + n(s.n) + ' variants</span>'; }).join('<br>') : '<span class="dim">loading the scores...</span>') + '</div>';
    // GWAS loci
    var gw = G.app.gwas;
    html += '<div style="margin-top:6px"><b>GWAS Catalog</b> <span class="dim">(depression traits, p &le; 5e-8)</span></div>';
    if (gw) {
      // depression itself only: interaction tests (lipids x depression), combined disorders
      // ("bipolar or major depressive") and MTAG analyses are left out
      var keep = /^(lifetime |recurrent )?(major depressive disorder|depression|depressive symptoms|depressed affect)( \((broad|narrow)\))?$/i;
      var traits = gw.searchTraits('depress', 200).filter(function (t) { return keep.test(t.trait.trim()); }), loci = gw.lociFor(traits.map(function (t) { return t.idx; }));
      var carried = 0, readable = 0, copies = 0, rows = [];
      loci.forEach(function (l) {
        var ds = gw.dosage(d, l.k, l.i, 0);
        if (ds.dosage !== null) { readable++; copies += ds.dosage; if (ds.dosage > 0) carried++; }
        if (rows.length < 10) rows.push({ l: l, ds: ds });
      });
      html += '<div>' + n(loci.length) + ' loci across ' + traits.length + ' depression traits (' + esc(traits.map(function (t) { return t.trait; }).join('; ')) + '; interaction, combined-disorder and MTAG traits left out). Readable here: ' + n(readable) + '; this genome carries the reported risk allele at ' + n(carried) + ' of them (' + n(copies) + ' copies). <span class="dim">Risk alleles at these loci are common: most people carry many, and each moves risk very little.</span></div>' +
        rows.map(function (x) {
          var g = gw.byContig[x.l.k];
          return '<div class="dim">&nbsp; <a href="#" data-dgw="' + x.l.k + ':' + g.pos[x.l.i] + '">' + esc(g.rsid[x.l.i]) + '</a> ' + esc(g.gene[x.l.i] || '') + ' risk ' + esc(g.risk[x.l.i]) + ', p ' + x.l.p.toExponential(0) + ': ' + (x.ds.dosage === null ? 'unknown (' + esc(x.ds.reason) + ')' : x.ds.dosage + ' copies') + '</div>';
        }).join('');
    } else html += '<div class="dim">GWAS Catalog not loaded (python3 tools/fetch_annotations.py).</div>';
    html += '<div class="warn" style="margin-top:6px">Not a diagnosis and not a prescription. Depression is common and mostly not predictable from DNA; drug choice and dose belong to a clinician, with the CPIC guidelines.</div>';
    out.innerHTML = html;
    Array.prototype.forEach.call(out.querySelectorAll('[data-dpgs]'), function (el) { el.onclick = function (e) { e.preventDefault(); loadScore(el.dataset.dpgs).then(function () { renderDepression(true); }); }; });
    Array.prototype.forEach.call(out.querySelectorAll('[data-dgw]'), function (el) { el.onclick = function (e) { e.preventDefault(); var x = el.dataset.dgw.split(':'); G.app.setMode('arcs'); view.goTo(current.genome.get(x[0]).name, +x[1] - 5000, +x[1] + 5000); }; });
    if (!depressionScores) G.prs.searchTraits('major depressive').then(async function (ts) {
      var ids = [];
      ts.forEach(function (t) { if (/major depressive disorder/i.test(t.label)) t.scores.forEach(function (id) { if (ids.indexOf(id) < 0) ids.push(id); }); });
      var infos = await Promise.all(ids.slice(0, 6).map(function (id) { return G.prs.scoreInfo(id).catch(function () { return null; }); }));
      depressionScores = infos.filter(Boolean).map(function (s) { return { id: s.id, name: s.name, n: s.n }; });
      renderDepression(true);
    }, function () { depressionScores = []; $('depScores').innerHTML = '<span class="warn">PGS Catalog unreachable</span>'; });
  }

  // ----- methylation overlay (methylation.js): only when a methylation file is opened
  async function loadMethylation(file) {
    var box = $('methylBox');
    if (!current || current.format !== 'vcf' && current.format !== 'bam') { box.innerHTML = '<div class="warn">Open the genome (VCF, gVCF or BAM) first, then the methylation file: it is drawn on that genome.</div>'; return; }
    box.innerHTML = '<div class="fhead">Methylation</div><div class="dim">reading ' + esc(file.name) + '...</div>';
    try {
      var m = await G.methylation.parse(file, current, { name: file.name, onProgress: function (st) { box.innerHTML = '<div class="fhead">Methylation</div><div class="dim">reading ' + esc(file.name) + ': ' + n(st.sites) + ' sites...</div>'; } });
      // a haplotype file pairs with the other haplotype already loaded for this genome
      var prev = current.methyl, other = prev && (prev.hap ? prev : null);
      if (m.hap && other && other.hap !== m.hap) m = G.methylation.combine(other, m);
      G.app.methyl = view.methyl = current.methyl = m;
      describe(current); renderMethylation();
      if (G.app.onExtrasChanged) G.app.onExtrasChanged();
    } catch (err) {
      console.error(err);
      box.innerHTML = '<div class="fhead">Methylation</div><div class="warn">' + esc(err.message) + '</div>';
    }
  }
  function renderMethylation() {
    var m = G.app.methyl, box = $('methylBox');
    if (!m || !current || current.methyl !== m) { box.innerHTML = ''; return; }
    var s = m.summary, pc = function (x) { return x === null ? '-' : Math.round(100 * x) + '%'; };
    var rows = current.build === 'GRCh38' ? m.checks(genes) : []; // the gene table is GRCh38: never read it on another build
    var hapNote = m.parts ? '<div class="dim">Two haplotypes: copies are compared within phase blocks, so which copy is 1 or 2 is arbitrary; the difference is what counts.</div>'
      : m.hap ? '<div class="dim">Haplotype ' + m.hap + ' only; open the other haplotype file to compare the two copies.</div>' : '';
    box.innerHTML = '<div class="fhead">Methylation <span class="dim">' + esc(m.fileName) + ', ' + esc(m.stats.format) + ', ' + n(m.stats.sites) + ' sites</span></div>' +
      (m.codes.length > 1 ? '<div>Mark: ' + m.codes.map(function (c) { return '<a href="#" data-mcode="' + c + '"' + (c === m.code ? ' style="font-weight:bold"' : '') + '>' + esc(G.methylation.Methylation.CODE_NAMES[c] || c) + '</a>'; }).join(' ') + '</div>' : '') +
      hapNote + '<div>Genome-wide ' + pc(s.mean) + ' methylated; of the well-covered 1 kb regions ' + pc(s.low) + ' low (under 20%), ' + pc(s.mid) + ' middle, ' + pc(s.high) + ' high (over 70%). <span class="dim">Blood and tissue: most regions high, CpG island promoters low. Cell lines (HG002 DNA comes from a lymphoblastoid line) often show many partly methylated regions.</span></div>' +
      (rows.length ? '<div style="margin-top:4px">Checks <span class="dim">(' + (m.parts ? 'per copy; promoter window or the locus scan below' : 'promoter windows, a stand-in for the exact control regions') + ')</span>:</div>' + rows.map(function (r) {
        return '<div><a href="#" data-mgene="' + esc(r.gene) + '">' + esc(r.gene) + '</a> ' + (r.hap ? 'copies ' + pc(r.hap[0]) + ' / ' + pc(r.hap[1]) : pc(r.frac)) + ' <span class="' + (r.unusual ? 'warn' : 'dim') + '">' + esc(r.flag) + '</span> <span class="dim">' + esc(r.role) + (r.where ? ', ' + esc(r.where) : '') + '</span></div>';
      }).join('') : '<div class="dim">Imprinting and FMR1 checks need a GRCh38 genome (the gene table is GRCh38)' + (genes ? '' : ' and the gene table') + '.</div>') +
      (m.parts ? '<div class="dim">With two haplotypes each imprinted locus (gene +-10 kb) is also scanned for its strongest stretch where the copies differ by 0.4 or more. Haplotype files hold only phased reads: unphased stretches, and a male X, have no calls.</div>' : '') +
      '<div class="dim">These are levels, not diagnoses. Imprinting and Fragile X are confirmed with dedicated assays.</div>';
    Array.prototype.forEach.call(box.querySelectorAll('[data-mcode]'), function (el) { el.onclick = function (e) { e.preventDefault(); m.code = el.dataset.mcode; m.summary = m.globalSummary(); renderMethylation(); if (G.app.onExtrasChanged) G.app.onExtrasChanged(); }; });
    Array.prototype.forEach.call(box.querySelectorAll('[data-mgene]'), function (el) {
      el.onclick = function (e) { e.preventDefault(); var r = rows.find(function (x) { return x.gene === el.dataset.mgene; }); G.app.setMode('arcs'); view.goTo(r.chrom, r.start - 3000, r.end + 3000); };
    });
  }

  // ----- runs of homozygosity (roh.js)
  function refreshRoh() {
    var r = G.app.roh = current && current.format === 'vcf' ? G.roh.call(current) : null, box = $('rohBox');
    if (!r) { box.innerHTML = ''; return; }
    var auto = r.segments.filter(function (s) { return G.roh.isAutosome(s.chrom); }).sort(function (a, b) { return b.length - a.length; });
    var other = r.segments.filter(function (s) { return !G.roh.isAutosome(s.chrom); });
    var mb = function (x) { return (x / 1e6).toFixed(1) + ' Mb'; };
    box.innerHTML = '<div class="fhead">Runs of homozygosity <span class="dim">' + (r.method === 'gvcf' ? 'het calls per callable kb, 100 kb windows' : 'het share of called sites, 1 Mb windows (exome: coarse)') + '</span></div>' +
      '<div>F<sub>ROH</sub> ' + (100 * r.fRoh).toFixed(2) + '% of the autosomes in runs of 1.5 Mb or more; ' + auto.length + ' autosomal runs of 1 Mb or more' + (r.longCount ? ', <b>' + r.longCount + ' over 5 Mb</b> (recent shared ancestry)' : '') + '.</div>' +
      (r.mostlyRoh.length ? '<div class="warn">Mostly homozygous: ' + esc(r.mostlyRoh.join(', ')) + ' (consider uniparental isodisomy).</div>' : '') +
      (other.length ? '<div class="dim">' + esc(other.map(function (s) { return s.chrom + ' ' + mb(s.length); }).join(', ')) + ': X and Y are single copy in a male, so long runs there are expected and left out of F<sub>ROH</sub>.</div>' : '') +
      auto.slice(0, 10).map(function (s, i) { return '<div><a href="#" data-roh="' + r.segments.indexOf(s) + '">' + esc(s.chrom) + ':' + G.fmtBp(s.start) + '-' + G.fmtBp(s.end) + '</a> ' + mb(s.length) + '</div>'; }).join('');
    Array.prototype.forEach.call(box.querySelectorAll('[data-roh]'), function (el) {
      el.onclick = function (e) { e.preventDefault(); var s = r.segments[+el.dataset.roh]; G.app.setMode('arcs'); view.goTo(s.chrom, s.start, s.end); };
    });
    if (G.app.onExtrasChanged) G.app.onExtrasChanged();
  }

  // ----- gene panels (panels.js)
  var chosenPanels = [];
  try { chosenPanels = JSON.parse(localStorage.getItem('genomeatrium.panels') || '[]'); } catch (e) { chosenPanels = []; }
  function savePanels() { try { localStorage.setItem('genomeatrium.panels', JSON.stringify(chosenPanels)); } catch (e) { /* not kept */ } }
  $('panelSearch').addEventListener('input', function () {
    var q = this.value;
    if (q.trim().length < 2) { $('panelResults').innerHTML = ''; return; }
    $('panelResults').innerHTML = '<span class="dim">loading the panel lists...</span>';
    G.panels.listPanels().then(function (all) {
      if ($('panelSearch').value !== q) return;
      var hits = G.panels.searchPanels(all, q);
      $('panelResults').innerHTML = hits.length ? hits.map(function (p, i) {
        return '<div><a href="#" data-pi="' + i + '">' + esc(p.name) + '</a> <span class="dim">' + p.nGenes + ' genes, ' + esc(p.sourceName) + ' v' + esc(p.version) + '</span></div>';
      }).join('') : '<span class="dim">no panel matches</span>';
      Array.prototype.forEach.call($('panelResults').querySelectorAll('[data-pi]'), function (el) {
        el.onclick = function (e) {
          e.preventDefault(); var p = hits[+el.dataset.pi];
          if (!chosenPanels.some(function (c) { return c.source === p.source && c.id === p.id; })) chosenPanels.push({ source: p.source, id: p.id });
          savePanels(); $('panelSearch').value = ''; $('panelResults').innerHTML = ''; refreshPanels();
        };
      });
    }, function (err) { $('panelResults').innerHTML = '<span class="warn">' + esc(err.message) + '</span>'; });
  });

  function refreshPanels() {
    var box = $('panelList');
    if (!chosenPanels.length) { box.innerHTML = '<div class="dim">Pick a panel to check its genes: coverage gaps, recessive genes in runs of homozygosity, findings.</div>'; G.app.panelRows = []; if (G.app.onExtrasChanged) G.app.onExtrasChanged(); return; }
    if (!current) { box.innerHTML = ''; return; }
    var same = current.build === 'GRCh38'; // PanelApp places genes on GRCh38; never compare across builds
    Promise.all(chosenPanels.map(function (c) { return G.panels.loadPanel(c.source, c.id).catch(function (e) { return { error: e.message, source: c.source, id: c.id }; }); })).then(function (panels) {
      var all = [];
      box.innerHTML = panels.map(function (p, pi) {
        if (p.error) return '<div class="warn">' + esc(p.error) + ' <a href="#" data-rm="' + pi + '">remove</a></div>';
        var a = G.panels.assess(same ? current : null, p, { genes: genes, roh: same ? G.app.roh : null, findings: view.findings });
        a.rows.forEach(function (r) { r.panel = p.name; all.push(r); });
        var s = a.summary, gene = function (r) { return '<a href="#" data-gene="' + esc(r.symbol) + '">' + esc(r.symbol) + '</a>'; };
        return '<div class="fhead" style="margin-top:6px">' + esc(p.name) + ' <span class="dim">' + esc(p.sourceName) + ' v' + esc(p.version) + '</span> <a href="#" data-rm="' + pi + '">remove</a></div>' +
          (!same ? '<div class="dim">This file is ' + esc(current.build || 'of unknown build') + '; panel genes are placed on GRCh38, so coverage and runs are not compared.</div>'
            : !a.canCover ? '<div class="dim">Coverage needs a gVCF (reference blocks); this file has variant sites only.</div>'
            : '<div>' + s.complete + ' of ' + s.green + ' green genes fully callable (95% of the span or more)</div>' +
              (s.gaps.length ? '<div>Gaps: ' + s.gaps.slice(0, 20).map(function (r) { return gene(r) + ' <span class="dim">' + Math.round(100 * r.callable) + '%</span>'; }).join(', ') + (s.gaps.length > 20 ? ' ...' : '') + '</div>' : '')) +
          (s.recessiveInRoh.length ? '<div>Recessive genes in a run of homozygosity: ' + s.recessiveInRoh.map(gene).join(', ') + '</div>' : '') +
          (s.withFindings.length ? '<div>Findings in panel genes: ' + s.withFindings.map(gene).join(', ') + '</div>' : '') +
          (s.unplaced ? '<div class="dim">' + s.unplaced + ' genes without a GRCh38 location</div>' : '');
      }).join('');
      G.app.panelRows = all;
      Array.prototype.forEach.call(box.querySelectorAll('[data-rm]'), function (el) { el.onclick = function (e) { e.preventDefault(); chosenPanels.splice(+el.dataset.rm, 1); savePanels(); refreshPanels(); }; });
      Array.prototype.forEach.call(box.querySelectorAll('[data-gene]'), function (el) { el.onclick = function (e) { e.preventDefault(); focusGene(el.dataset.gene, true); }; });
      if (G.app.onExtrasChanged) G.app.onExtrasChanged();
    });
  }
  G.app.refreshPanels = refreshPanels;

  function renderGwas() {
    if (!gwas) return;
    var q = $('gwasSearch').value, sel = G.app.gwasTraitLoci, html = [];
    if (sel) {
      var d = current, counts = [0, 0, 0], unknown = 0, viaEnh = 0;
      sel.loci.forEach(function (l) {
        var dz = d && d.variants ? gwas.dosage(d, l.k, l.i).dosage : null;
        if (dz === null) unknown++; else counts[dz]++;
      });
      html.push('<div style="margin-top:4px"><b>' + esc(sel.name) + '</b> <button id="gwasClear">clear</button><br><span class="dim">' + sel.loci.length + ' loci. Risk allele in this genome: 2 copies ' +
        counts[2] + ', 1 copy ' + counts[1] + ', none ' + counts[0] + ', unknown ' + unknown + '. A count, not a risk score.</span></div>');
      var seenRs = {};
      sel.loci.forEach(function (l, i) { l.shown = false; });
      sel.loci.filter(function (l) { var r = gwas.byContig[l.k].rsid[l.i]; if (seenRs[r]) return false; seenRs[r] = 1; return true; }).slice(0, 40).forEach(function (l) {
        var i = sel.loci.indexOf(l);
        var gc = gwas.byContig[l.k], dz = d && d.variants ? gwas.dosage(d, l.k, l.i).dosage : null, enh = '';
        var reg = G.app.reg;
        if (reg && reg.links) {
          var t = reg.linksInRange(l.k, gc.pos[l.i], gc.pos[l.i]).filter(function (x) { return x.start <= gc.pos[l.i] && x.end >= gc.pos[l.i] && !x.self; });
          if (t.length) { viaEnh++; enh = ' <span class="warn">enhancer of ' + esc(t.map(function (x) { return x.gene; }).filter(function (v, j, a) { return a.indexOf(v) === j; }).slice(0, 3).join(', ')) + '</span>'; }
        }
        html.push('<div class="gl" data-i="' + i + '">' + (dz === null ? '?' : dz) + ' &middot; ' + esc(gc.rsid[l.i]) + ' <span class="dim">' + esc((gc.gene[l.i] || '').slice(0, 24)) + ' p ' + fmtP(l.p) + '</span>' + enh + '</div>');
      });
    } else if (q.trim().length >= 3) {
      gwas.searchTraits(q, 25).forEach(function (t) { html.push('<div class="gl" data-t="' + t.idx + '">' + esc(t.trait) + ' <span class="dim">' + t.n + ' loci</span></div>'); });
      if (!html.length) html.push('<div class="dim">no trait matches</div>');
    }
    $('gwasList').innerHTML = html.join('');
    Array.prototype.forEach.call(document.querySelectorAll('#gwasList [data-t]'), function (el) {
      el.onclick = function () { pickTrait(+el.dataset.t); };
    });
    Array.prototype.forEach.call(document.querySelectorAll('#gwasList [data-i]'), function (el) {
      el.onclick = function () { var l = G.app.gwasTraitLoci.loci[+el.dataset.i], gc = gwas.byContig[l.k]; setMode('arcs'); view.goTo(current.genome.get(l.k).name, gc.pos[l.i] - 5000, gc.pos[l.i] + 5000); };
    });
    if ($('gwasClear')) $('gwasClear').onclick = function () { G.app.gwasTraitLoci = null; G.app.gwasGenes = new Set(); renderGwas(); };
  }
  function fmtP(p) { return p > 0 ? p.toExponential(0) : '< 1e-300'; }
  G.app.fmtP = fmtP;

  function pickTrait(idx) {
    var loci = gwas.lociFor([idx]), set = new Set(), genesSet = new Set();
    loci.forEach(function (l) { set.add(l.k + ':' + l.i); });
    // genes of loci where this genome carries a risk allele, for the pathway view
    loci.forEach(function (l) {
      var dz = current && current.variants ? gwas.dosage(current, l.k, l.i).dosage : null;
      if (dz) String(gwas.byContig[l.k].gene[l.i]).split(/[,;-]\s*|\s+-\s+/).forEach(function (g) { g = g.trim(); if (g && view.genes && view.genes.get(g)) genesSet.add(g); });
    });
    G.app.gwasTraitLoci = { name: gwas.traits[idx][0], loci: loci, set: set };
    G.app.gwasGenes = genesSet;
    renderGwas();
  }

  // ncRNA: antisense pairs (GENCODE) with GeneChords results, and variants in small RNA genes.
  var coordination = null;
  fetch('data/annotations/genechords_antisense.csv').then(function (r) { return r.ok ? r.text() : null; }).then(function (t) {
    if (!t) return;
    var lines = t.trim().split('\n'), h = lines[0].split(',');
    coordination = lines.slice(1).map(function (l) { var f = l.split(','), o = {}; h.forEach(function (k, i) { o[k] = f[i]; }); return o; });
    refreshNc();
  }).catch(function () {});
  function refreshNc() {
    if (!view.genes) return;
    var pairs = view.genes.antisensePairs(coordination), withCoord = pairs.filter(function (p) { return p.coord; });
    var html = ['<div class="fhead">ncRNA <span class="dim">' + pairs.length.toLocaleString() + ' lncRNA antisense pairs (GENCODE)</span></div>'];
    var d = current;
    if (d && d.format === 'vcf' && d.build === 'GRCh38') {
      var hits = [];
      Object.keys(view.genes.byContig).forEach(function (k) {
        var c = d.variants[k];
        if (!c) return;
        view.genes.byContig[k].list.forEach(function (gn) {
          if (gn.cls !== 'smallRNA') return;
          var lo = 0, hi = c.n;
          while (lo < hi) { var mid = (lo + hi) >> 1; if (c.pos[mid] < gn.start) lo = mid + 1; else hi = mid; }
          for (var i = lo; i < c.n && c.pos[i] <= gn.end; i++) hits.push({ gene: gn, k: k, i: i });
        });
      });
      // functional small RNAs first (miRNA, snoRNA, snRNA); misc_RNA is mostly 7SK/7SL copies
      var rank = { miRNA: 0, snoRNA: 1, snRNA: 2, scaRNA: 3 };
      hits.sort(function (a, b) { return (rank[a.gene.type] === undefined ? 9 : rank[a.gene.type]) - (rank[b.gene.type] === undefined ? 9 : rank[b.gene.type]); });
      G.app.smallRnaHits = hits;
      var byType = {};
      hits.forEach(function (h) { byType[h.gene.type] = (byType[h.gene.type] || 0) + 1; });
      html.push('<div class="dim">sample variants inside small RNA genes: ' + Object.keys(byType).sort(function (a, b) { return byType[b] - byType[a]; }).map(function (t) { return t + ' ' + byType[t]; }).join(', ') + '</div>');
      hits.slice(0, 12).forEach(function (h, j) {
        var c = d.variants[h.k], al = c.alleles(h.i), r = G.app.gnomad && G.app.gnomad.enabled ? G.app.gnomad.forVariant(h.gene.chrom, c, h.i) : undefined;
        html.push('<div class="gl" data-sr="' + j + '">' + esc(h.gene.name) + ' <span class="dim">' + h.gene.type + ' &middot; ' + (al ? esc(al.ref.slice(0, 6) + '>' + al.alts.join(',').slice(0, 8)) : '') + ' &middot; ' +
          G.vcf.ZYG_NAMES[c.zyg[h.i]] + (r !== undefined ? ' &middot; ' + G.gnomad.rarity(r).text : '') + '</span></div>');
      });
    }
    if (withCoord.length) {
      html.push('<div class="dim" style="margin-top:4px">GeneChords coordinated pairs (gold chords in Arcs)</div>');
      withCoord.sort(function (a, b) { return +b.coord.rho - +a.coord.rho; }).slice(0, 12).forEach(function (p, j) {
        html.push('<div class="gl" data-gc="' + j + '">' + esc(p.lnc.name) + ' / ' + esc(p.gene.name) + ' <span class="dim">rho ' + (+p.coord.rho).toFixed(2) + ', ' + esc(p.coord.marks) + '</span></div>');
      });
    }
    $('ncBox').innerHTML = html.join('');
    var sorted = withCoord;
    Array.prototype.forEach.call(document.querySelectorAll('[data-gc]'), function (el) {
      el.onclick = function () { var p = sorted[+el.dataset.gc]; setMode('arcs'); view.goTo(p.gene.chrom, Math.min(p.lnc.start, p.gene.start) - 2000, Math.max(p.lnc.end, p.gene.end) + 2000); };
    });
    Array.prototype.forEach.call(document.querySelectorAll('[data-sr]'), function (el) {
      el.onclick = function () { var h = G.app.smallRnaHits[+el.dataset.sr]; setMode('arcs'); view.goTo(h.gene.chrom, h.gene.start - 200, h.gene.end + 200); };
    });
  }
  G.app.refreshNc = refreshNc;

  var pathways = G.app.pathways = new G.Pathways();
  fetch('data/annotations/reactome.json').then(function (r) { return r.ok ? r.json() : null; }).then(function (doc) {
    if (!doc) return;
    var tryLoad = function () { if (view.genes) pathways.load(doc, view.genes); else setTimeout(tryLoad, 500); };
    tryLoad();
  }).catch(function () {});
  ['Findings', 'Gwas', 'Focus'].forEach(function (k) { $('ps' + k).onchange = function () { pathways.sets[k.toLowerCase()] = $('ps' + k).checked; pathways.scoredFor = null; }; });

  G.app.mergeFindings = mergeFindings;

  // ----- regulatory layer: GENCODE genes + ENCODE-rE2G links for picked tissues

  var reg = new G.regulatory.Regulatory(), genes = null;
  G.app.reg = reg; view.reg = reg;
  var SAVE_KEY = 'genomeatrium.tissues', MODEL_KEY = 'genomeatrium.models';
  var models = (function () { try { return JSON.parse(localStorage.getItem(MODEL_KEY) || '["e2g"]'); } catch (e) { return ['e2g']; } })();
  $('modelE2g').checked = models.indexOf('e2g') >= 0; $('modelAbc').checked = models.indexOf('abc') >= 0;
  var chosen = new Set((function () { try { return JSON.parse(localStorage.getItem(SAVE_KEY) || '[]'); } catch (e) { return []; } })());
  function okBlob(r) { if (!r.ok) throw new Error(r.status); return r.blob(); }
  function regNote(html) { $('regNote').innerHTML = html; }

  Promise.all([
    fetch('data/regulatory/genes_grch38.tsv.gz').then(okBlob).then(function (b) { return new G.regulatory.Genes().load(b); })
      .then(function (gn) { genes = gn; view.genes = gn; refreshNc(); refreshPanels(); }),
    fetch('data/regulatory/catalog.json').then(okBlob).then(function (b) { return reg.loadCatalog(b); })
  ]).then(function () {
    $('pickerSource').textContent = reg.catalog.sets.length.toLocaleString() + ' prediction sets from ' + reg.biosamples.length + ' biosamples (' + reg.catalog.source + '). Genes: ' + (genes.meta.source || 'GENCODE').replace(/\s*\(.*\)/, '') + '.';
    if (chosen.size) applyTissues(); else refreshRegulatory();
  }).catch(function () {
    regNote('Regulatory data not found. Run <code>python3 tools/fetch_regulatory.py</code> (genes and tissue catalog).');
  });

  async function applyTissues() {
    try {
      var done = await reg.pick(Array.from(chosen), models, function (t) { regNote(esc(t) + '...'); });
      if (!done) return; // superseded by a newer pick
      if (current) reg.intersect(current);
      refreshRegulatory();
    } catch (err) {
      console.error(err);
      regNote('<span class="warn">' + esc(err.message) + '</span>');
    }
  }

  function geneSpan(name) {
    var gn = genes && genes.get(name), links = reg.byGene ? reg.byGene.get(gn ? gn.name : name) || [] : [];
    if (!gn && !links.length) return null;
    var chrom = gn ? gn.chrom : links[0].chrom;
    var a = gn ? gn.start : Infinity, b = gn ? gn.end : -Infinity;
    links.forEach(function (l) { a = Math.min(a, l.start, l.tss); b = Math.max(b, l.end, l.tss); });
    var pad = Math.max(2000, (b - a) * 0.08);
    return { name: gn ? gn.name : name, chrom: chrom, start: Math.max(1, Math.round(a - pad)), end: Math.round(b + pad) };
  }

  function focusGene(name, jump) {
    view.focusGene = name;
    if (G.app.proteinView) G.app.proteinView.gene = null; // the protein view follows the focused gene
    if (name && jump && (view.mode === 'gene' || view.mode === 'protein')) jump = false;
    if (G.app.gnomadWants) setTimeout(G.app.gnomadWants, 0);
    if (name && jump) {
      var sp = geneSpan(name);
      if (sp) { setMode('arcs'); view.goTo(sp.chrom, sp.start, sp.end); }
    }
    refreshRegulatory();
  }
  G.app.focusGene = focusGene;
  G.app.onFocusGene = function () { refreshRegulatory(); if (G.app.gnomadWants) G.app.gnomadWants(); };

  function refreshRegulatory() {
    var d = current, list = [];
    if (!reg.catalog) return;
    if (!chosen.size) { regNote('No tissues picked. Enhancer and promoter links are per tissue: pick one or more.'); $('regList').innerHTML = ''; return; }
    var sources = {};
    (reg.used || []).forEach(function (u) { sources[u.from] = (sources[u.from] || 0) + 1; });
    var byModel = { e2g: 0, abc: 0, both: 0 };
    (reg.links || []).forEach(function (l) { if (l.agree) byModel.both++; else if (l.scores.e2g !== undefined) byModel.e2g++; else byModel.abc++; });
    var head = esc(Array.from(chosen).join(', ')) + '<br>' + (reg.links ? reg.links.length.toLocaleString() + ' links from ' + (reg.used || []).length + ' files (' +
      Object.keys(sources).map(function (k) { return sources[k] + ' ' + k; }).join(', ') + ')' : '');
    if (reg.links && (reg.models || []).length > 1) head += '<br>rE2G only ' + byModel.e2g.toLocaleString() + ', ABC only ' + byModel.abc.toLocaleString() +
      ', both <b>' + byModel.both.toLocaleString() + '</b>';
    if (d && d.build && d.build !== 'GRCh38') head += '<br><span class="warn">Not drawn: the file is ' + esc(d.build) + ', the links are GRCh38.</span>';
    else if (reg.hits) head += '<br>' + reg.hits.links.toLocaleString() + ' linked elements carry ' + reg.hits.variants.toLocaleString() + ' sample variants';
    regNote(head);
    // finding genes first: do their regulatory elements carry variants?
    var findingGenes = [];
    (view.findings || []).forEach(function (f) { String(f.gene || '').split(/[;,]/).forEach(function (x) { if (x && findingGenes.indexOf(x) < 0) findingGenes.push(x); }); });
    if (reg.byGene && findingGenes.length) {
      list.push('<div class="dim" style="margin-top:4px">finding genes: links (with a sample variant)</div>');
      findingGenes.forEach(function (gname) {
        var ls = reg.byGene.get(gname) || [], hit = ls.filter(function (l) { return l.variants; }).length;
        list.push('<div class="rgene' + (view.focusGene === gname ? ' on' : '') + '" data-g="' + esc(gname) + '"><b>' + esc(gname) + '</b> <span class="dim">' +
          ls.length + ' (' + hit + ')</span></div>');
      });
    }
    var fg = view.focusGene;
    if (fg && reg.byGene) {
      var ls = (reg.byGene.get(fg) || []).slice().sort(function (a, b) { return (b.agree - a.agree) || b.score - a.score; });
      list.push('<div class="dim" style="margin-top:6px">' + esc(fg) + ': ' + ls.length + ' links <button id="releaseGene">release</button></div>');
      ls.slice(0, 40).forEach(function (l, i) {
        list.push('<div class="rgene rlink" data-i="' + i + '"><span class="rsw" style="background:' + G.REG_COLORS[l.cls] + '"></span>' +
          (l.self ? 'own promoter' : l.cls === 'intergenic' ? 'distal enhancer' : l.cls === 'genic' ? 'intragenic enhancer' : 'promoter') +
          modelBadge(l) + ' <span class="dim">' + (l.self ? '' : G.fmtBp(Math.abs(l.tss - l.mid)) + ' &middot; ') + scoreText(l) + ' &middot; ' +
          Object.keys(l.tissues).length + ' tissue' + (Object.keys(l.tissues).length > 1 ? 's' : '') + '</span>' +
          (l.variants ? ' <span class="warn">' + l.variants.length + ' variant' + (l.variants.length > 1 ? 's' : '') + '</span>' + rareCount(l) : '') + '</div>');
      });
      $('regList').innerHTML = list.join('');
      Array.prototype.forEach.call(document.querySelectorAll('.rlink'), function (el) {
        el.onclick = function (e) { e.stopPropagation(); var l = ls[+el.dataset.i]; setMode('arcs'); view.goTo(l.chrom, Math.min(l.start, l.tss) - 500, Math.max(l.end, l.tss) + 500); };
      });
      $('releaseGene').onclick = function () { focusGene(null); };
    } else {
      list.push('<div class="dim" style="margin-top:4px">Type a gene in the box above, or click one in the gene lane, to list its enhancers.</div>');
      $('regList').innerHTML = list.join('');
    }
    Array.prototype.forEach.call(document.querySelectorAll('.rgene[data-g]'), function (el) {
      el.onclick = function () { focusGene(el.dataset.g === view.focusGene ? null : el.dataset.g, true); };
    });
  }
  G.app.refreshRegulatory = refreshRegulatory;

  function modelBadge(l) {
    return l.agree ? '<span class="mb both" title="predicted by both ENCODE-rE2G and ABC">E+A</span>' :
      l.scores.e2g !== undefined ? '<span class="mb" title="ENCODE-rE2G only">E</span>' : '<span class="mb" title="ABC only">A</span>';
  }
  function scoreText(l) {
    return [l.scores.e2g !== undefined ? 'rE2G ' + l.scores.e2g.toFixed(2) : null, l.scores.abc !== undefined ? 'ABC ' + l.scores.abc.toFixed(3) : null]
      .filter(Boolean).join(' / ');
  }

  // "(1 rare)" after a link's variant count, once its gnomAD tile is in.
  function rareCount(l) {
    if (!gnomad.enabled || !current || !current.variants) return '';
    var cols = current.variants[l.key], rare = 0, unknown = 0;
    l.variants.forEach(function (i) {
      var r = gnomad.forVariant(l.chrom, cols, i);
      if (r === undefined) unknown++; else if (r.absent || r.af < 0.01) rare++;
    });
    return unknown ? ' <span class="af">AF ...</span>' : rare ? ' <span class="af rare">' + rare + ' rare</span>' : ' <span class="af common">all common</span>';
  }

  // Tissue picker
  function renderPicker() {
    if (!reg.biosamples) return;
    var q = $('pickerSearch').value.trim().toLowerCase(), cls = $('pickerClass').value;
    var groups = {};
    reg.biosamples.forEach(function (b) {
      if (cls && b.classification !== cls) return;
      var hay = (b.biosample + ' ' + b.organs.join(' ') + ' ' + b.classification).toLowerCase();
      if (q && hay.indexOf(q) < 0) return;
      var org = b.organ;
      (groups[org] = groups[org] || []).push(b);
    });
    var html = Object.keys(groups).sort().map(function (org) {
      return '<div class="org">' + esc(org) + '</div>' + groups[org].map(function (b) {
        return '<label><input type="checkbox" data-b="' + esc(b.biosample) + '"' + (chosen.has(b.biosample) ? ' checked' : '') + '>' + esc(b.biosample) +
          ' <span class="dim">' + esc(b.classification) + ', ' + b.sets.length + ' set' + (b.sets.length > 1 ? 's' : '') + '</span></label>';
      }).join('');
    }).join('');
    $('pickerList').innerHTML = html || '<span class="dim">nothing matches</span>';
    $('pickerChosen').textContent = chosen.size ? 'selected: ' + Array.from(chosen).join(', ') : 'nothing selected';
    Array.prototype.forEach.call(document.querySelectorAll('#pickerList input'), function (el) {
      el.onchange = function () { if (el.checked) chosen.add(el.dataset.b); else chosen.delete(el.dataset.b); $('pickerChosen').textContent = chosen.size ? 'selected: ' + Array.from(chosen).join(', ') : 'nothing selected'; };
    });
  }
  $('openPicker').onclick = function () { $('picker').hidden = false; renderPicker(); $('pickerSearch').focus(); };
  $('closePicker').onclick = function () { $('picker').hidden = true; };
  $('pickerSearch').oninput = renderPicker;
  $('pickerClass').onchange = renderPicker;
  $('pickerClear').onclick = function () { chosen.clear(); renderPicker(); };
  $('pickerApply').onclick = function () {
    models = [];
    if ($('modelE2g').checked) models.push('e2g');
    if ($('modelAbc').checked) models.push('abc');
    if (!models.length) { models = ['e2g']; $('modelE2g').checked = true; }
    try { localStorage.setItem(SAVE_KEY, JSON.stringify(Array.from(chosen))); localStorage.setItem(MODEL_KEY, JSON.stringify(models)); } catch (e) { /* private window: selection is not kept */ }
    $('picker').hidden = true;
    applyTissues();
  };


  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }
  function n(x) { return Number(x).toLocaleString(); }

  function describe(d) {
    var rows = ['<span class="file">' + esc(d.fileName) + '</span>'];
    var meta = [d.format.toUpperCase() + (d.isGvcf ? ' (gVCF)' : ''), d.build || 'build unknown', 'read in ' + d.seconds.toFixed(1) + ' s'];
    rows.push('<span class="dim">' + meta.join(' &middot; ') + '</span>');
    if (d.aborted) rows.push('<span class="warn">Stopped early: showing the part read so far.</span>');
    var hidden = d.genome.hidden || [];
    if (hidden.length) rows.push('<span class="dim">' + hidden.length + ' unplaced/alt/decoy contigs with data not drawn</span>');
    var legend = [];

    if (d.format === 'vcf') {
      var s = d.stats, T = G.vcf.T, Z = G.vcf.Z;
      if (d.samples.length) rows.push('sample' + (d.samples.length > 1 ? 's' : '') + ': ' + esc(d.samples.slice(0, 3).join(', ')) + (d.samples.length > 3 ? ' +' + (d.samples.length - 3) : '') +
        (d.samples.length > 1 ? ' <span class="dim">(genotypes from the first)</span>' : ''));
      rows.push(n(s.variants) + ' variant records, ' + n(s.filtered) + ' failing FILTER' + (s.refBlocks ? ', ' + n(s.refBlocks) + ' reference blocks' : ''));
      rows.push('<span class="dim">het ' + n(s.byZyg[Z.HET]) + ' &middot; hom ' + n(s.byZyg[Z.HOM]) + ' &middot; missing ' + n(s.byZyg[Z.MISSING]) + ' &middot; 0/0 ' + n(s.byZyg[Z.REF]) + '</span>');
      legend.push(['snv', G.COLORS.snv, 'SNV/MNV row (' + n(s.byType[T.SNV] + s.byType[T.MNV]) + ')']);
      legend.push(['indel', G.COLORS.indel, 'indel row (' + n(s.byType[T.INS] + s.byType[T.DEL] + s.byType[T.COMPLEX]) + ')']);
      legend.push(['het', 'rgb(190,140,255)', 'het fraction row: dips are homozygous runs or hemizygous X']);
      legend.push(['similar', 'hsl(200,70%,62%)', similarText()]);
      legend.push(['arcs', G.COLORS.sv, 'SV arcs: SVs, breakends, indels of 50 bp or more (' + n(view.arcs.length) + ')']);
      legend.push(['ticks', 'rgb(255,255,255)', 'variant ticks when zoomed in (tall = hom, short = het, hollow = 0/0 or missing)']);
      if (d.isGvcf) legend.push(['callable', 'rgb(230,170,40)', 'callable strip: grey called, amber low depth, red not called']);
      if (s.labelsDropped) rows.push('<span class="dim">Over 1.5M variants: hover details kept for the first 1.5M only.</span>');
    } else if (d.stats.aligned) {
      var b = d.stats;
      rows.push(n(b.reads) + ' records &middot; ' + n(b.mapped) + ' mapped &middot; ' + n(b.dup) + ' duplicates');
      if (b.paired) rows.push('<span class="dim">insert median ' + d.insert.median + ' bp, discordant cut-off ' + n(d.insert.cutoff) + ' bp</span>');
      var junc = d.arcs.filter(function (a) { return a.type === 'junction'; }).length;
      legend.push(['depth', G.COLORS.depth, 'mean depth row']);
      legend.push(['similar', 'hsl(200,70%,62%)', similarText()]);
      legend.push(['arcs', G.COLORS.arc.pair_inter, 'pairs across chromosomes (' + n(b.discordantInter) + ')']);
      legend.push(['arcs', G.COLORS.arc.pair_long, 'pairs with long insert (' + n(b.discordantLong) + ')']);
      legend.push(['arcs', G.COLORS.arc.junction, 'splice junctions (' + n(junc) + ' distinct, ' + n(b.spliced) + ' spliced reads)']);
      if (b.cgCigars) rows.push('<span class="dim">' + n(b.cgCigars) + ' reads with CIGAR in CG tag</span>');
    } else {
      rows.push(n(d.stats.reads) + ' unaligned reads');
    }
    if (d.format === 'vcf') {
      if (d.svFile) rows.push('<span class="dim">+ ' + n(d.svAdded) + ' structural variants from ' + esc(d.svFile) + '</span>');
      rows.push('<span class="dim">ClinVar: ' + esc(d.clinvarStatus) + '</span>');
      var cal = G.app.clinvar && d.clinvarStatus === 'matched' ? G.app.clinvar.callability(d) : null;
      if (cal && cal.total) rows.push('<span class="dim">ClinVar sites in this gVCF: ' + pct(cal.called, cal.total) + ' called, ' +
        pct(cal.lowdp, cal.total) + ' low depth, <span class="warn">' + pct(cal.nocall, cal.total) + ' no call</span> (a missed finding there is not a negative)</span>');
      if (d.clinvarStatus === 'matched') legend.push(['clinvar', 'rgb(255,80,80)', 'ClinVar P/LP site density (strip under the line)']);
      if (d.build === 'GRCh38') {
        legend.push(['genes', 'rgba(255,255,255,0.6)', 'genes (GENCODE), under the line']);
        legend.push(['lncRNA', 'rgb(80,210,190)', 'lncRNA genes']);
        legend.push(['smallRNA', 'rgb(255,110,190)', 'small RNA genes (miRNA, snoRNA, snRNA...)']);
        legend.push(['pseudogene', 'rgba(255,255,255,0.3)', 'pseudogenes']);
        legend.push(['antisense', 'rgb(80,210,190)', 'antisense chords under the line (gold: GeneChords coordinated pairs)']);
        legend.push(['gwas', 'rgb(190,140,255)', 'GWAS Catalog SNPs: strip, and diamonds when zoomed (fill = risk allele copies here)']);
        legend.push(['regulatory', 'rgb(255,190,70)', 'enhancer arcs: amber distal, blue intragenic, pink promoter; solid rE2G, dashed ABC, thick both; white box = sample variant in element, red = a rare one (gnomAD on)']);
      }
    }
    if (d.prs) legend.push(['prs', 'linear-gradient(90deg,rgb(110,170,255),rgb(255,110,90))', 'polygenic score band: net contribution per region (' + esc(d.prs.info.id) + '; red raises, blue lowers)']);
    if (d.methyl) legend.push(['methyl', 'linear-gradient(90deg,rgb(60,160,255),rgb(255,60,50))', 'methylation band: blue unmethylated to red methylated (' + esc(G.methylation.Methylation.CODE_NAMES[d.methyl.code] || d.methyl.code) + ')']);
    $('summary').innerHTML = rows.join('<br>');

    $('legend').innerHTML = legend.map(function (l) {
      return '<label><input type="checkbox" data-layer="' + l[0] + '"' + (view.layers[l[0]] ? ' checked' : '') +
        '><span class="sw" style="background:' + l[1] + '"></span>' + l[2] + '</label>';
    }).join('');
    Array.prototype.forEach.call(document.querySelectorAll('#legend input'), function (el) {
      el.addEventListener('change', function () { setLayer(el.dataset.layer, el.checked); });
    });
  }

  // ----- controls

  // The side panel: minimised to its heading or open, the same in every view, remembered.
  function setPanelMin(min) {
    $('info').classList.toggle('min', min);
    $('toggleMark').textContent = min ? 'show panel' : 'hide panel';
    try { localStorage.setItem('genomeatrium.panel', min ? 'min' : 'open'); } catch (e) { /* not kept */ }
  }
  (function () { var m = false; try { m = localStorage.getItem('genomeatrium.panel') === 'min'; } catch (e) { /* open */ } setPanelMin(m); })();
  $('toggleInfo').onclick = function () { setPanelMin(!$('info').classList.contains('min')); };
  document.addEventListener('keydown', function (e) {
    var t = e.target && e.target.tagName;
    if (t === 'INPUT' || t === 'TEXTAREA' || t === 'SELECT' || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === 'i') setPanelMin(!$('info').classList.contains('min'));
    if (e.key === 'Escape' && view.mode === 'atrium' && atrium) atrium.escape();
  });
  $('placeBack').onclick = function () { if (atrium) atrium.goPlace('atrium'); };
  $('open').onclick = function () { $('file').click(); };
  $('file').onchange = function (e) { loadMany(Array.prototype.slice.call(e.target.files)); e.target.value = ''; };
  $('addPerson').onclick = function () { $('fileAdd').click(); };
  $('fileAdd').onchange = function (e) {
    var fs = Array.prototype.slice.call(e.target.files); e.target.value = '';
    fs.reduce(function (p, f) { return p.then(function () { return addPerson(f); }); }, Promise.resolve());
  };
  $('stop').onclick = function () { if (controller) controller.abort(); };

  document.addEventListener('dragover', function (e) { e.preventDefault(); document.body.classList.add('dragging'); });
  document.addEventListener('dragleave', function (e) { if (!e.relatedTarget) document.body.classList.remove('dragging'); });
  document.addEventListener('drop', function (e) {
    e.preventDefault(); document.body.classList.remove('dragging');
    loadMany(Array.prototype.filter.call(e.dataTransfer.files, function (f) { return !/\.(tbi|csi|crai)$/i.test(f.name); }));
  });

  $('goto').addEventListener('keydown', function (e) {
    if (e.key !== 'Enter' || !current || !current.genome) return;
    var m = /^\s*([^:\s]+)(?::([\d,]+)(?:-([\d,]+))?)?\s*$/.exec(e.target.value);
    var ok = false;
    if (m) {
      var a = m[2] ? +m[2].replace(/,/g, '') : null, b = m[3] ? +m[3].replace(/,/g, '') : null;
      if (a && !b) { b = a + 500; a = Math.max(1, a - 500); }
      setMode('arcs');
      ok = view.goTo(m[1], a, b);
    }
    if (!ok && genes && !/:/.test(e.target.value)) { // a gene name
      var sp = geneSpan(e.target.value.trim());
      if (sp) { focusGene(sp.name, true); ok = true; }
    }
    e.target.style.borderColor = ok ? '' : 'rgb(255,120,100)';
  });
  document.addEventListener('keydown', function (e) {
    if (e.target.tagName === 'INPUT') return;
    if (e.key === 'r' && view.mode === 'arcs') view.reset();
  });

  // Hilbert map, Circos ring and Gene view (Moebio canvas).
  G.app.hilbert = new G.Hilbert(); G.app.hilbert.install(view.g.canvas);
  G.app.circos = new G.Circos();
  G.app.geneView = new G.GeneView();
  G.app.proteinView = new G.ProteinView();
  function renderHilbertLayers() {
    var h = G.app.hilbert;
    $('hilbertLayers').innerHTML = (h.layers || []).map(function (k) {
      return '<button data-l="' + k + '"' + (k === h.layer ? ' class="on"' : '') + '>' + G.HILBERT_LAYERS[k].label + '</button>';
    }).join('');
    Array.prototype.forEach.call(document.querySelectorAll('#hilbertLayers button'), function (b) {
      b.onclick = function () { h.setLayer(b.dataset.l); renderHilbertLayers(); };
    });
  }


  // Hi-C view: ENCODE contact maps read by range requests.
  var hic = G.app.hicView = new G.HicView(), lastArcsSpan = null;
  function hicRegion() {
    var sel = view.selection;
    if (sel) return sel;
    if (view.mode === 'atrium' && view.visibleSpan && view.zoom > 1.01) { // the Atrium window's region
      var vs = view.visibleSpan();
      if (vs && vs.end - vs.start < 12e6) return vs;
    }
    if (view.focusGene && G.app.reg && G.app.reg.byGene) {
      var gs = geneSpan(view.focusGene);
      if (gs) return { chrom: gs.chrom, start: gs.start, end: gs.end };
    }
    if (lastArcsSpan && lastArcsSpan.end - lastArcsSpan.start < 12e6) return lastArcsSpan;
    var f = (view.findings || [])[0];
    if (f) return { chrom: f.chrom, start: f.pos - 2.5e6, end: f.pos + 2.5e6 };
    return { chrom: current.genome.contigs[0].name, start: 1, end: 5e6 };
  }
  async function hicLoad() {
    if (!current || current.build !== 'GRCh38') return;
    hic.error = null;
    try {
      if (!hic.files) {
        var files = await hic.listFiles();
        $('hicFile').innerHTML = files.map(function (f) {
          return '<option value="' + f.acc + '">' + esc(f.biosample) + ' (' + f.acc + ', ' + Math.round(f.size / 1e9) + ' GB)</option>';
        }).join('');
        var gm = files.find(function (f) { return f.biosample === 'GM12878'; });
        if (gm) $('hicFile').value = gm.acc;
      }
      if (hic.file !== $('hicFile').value) {
        $('hicNote').textContent = 'opening ' + $('hicFile').value + '...';
        await hic.open($('hicFile').value);
        $('hicNorm').innerHTML = hic.norms.map(function (n) { return '<option' + (n === hic.norm ? ' selected' : '') + '>' + n + '</option>'; }).join('');
      }
      hic.norm = $('hicNorm').value || hic.norm;
      var r = hicRegion();
      $('hicNote').textContent = 'reading contacts for ' + r.chrom + ':' + G.fmtBp(r.start) + '-' + G.fmtBp(r.end) + '...';
      await hic.load(r);
      $('hicNote').textContent = 'only byte ranges of the file are read; contacts are from the chosen cell type, not the sample';
    } catch (err) {
      console.error(err);
      hic.error = err.message; $('hicNote').textContent = '';
    }
  }
  $('hicReload').onclick = hicLoad;
  $('hicFile').onchange = hicLoad;
  $('hicNorm').onchange = hicLoad;



  // One setter for every layer switch: legend, Atrium panel and VR card stay in step, and it is saved.
  function setLayer(key, on) {
    view.layers[key] = on;
    try { localStorage.setItem('genomeatrium.layers', JSON.stringify(view.layers)); } catch (e) { /* not kept */ }
    Array.prototype.forEach.call(document.querySelectorAll('input[data-layer="' + key + '"]'), function (el) { el.checked = on; });
  }
  G.app.setLayer = setLayer;
  function renderAtriumFilters() {
    $('atriumFilters').innerHTML = '<div class="fhead" id="filtersHead" style="cursor:pointer" title="Fold or open the filters">Filters <span class="dim" id="filtersMark"></span></div>' + G.Atrium.FILTERS.map(function (f) {
      return '<label><input type="checkbox" data-layer="' + f[0] + '"' + (view.layers[f[0]] !== false ? ' checked' : '') + '> ' + f[1] + '</label>';
    }).join('') + '<div class="dim" style="margin-top:4px">In VR: click the left stick for these on your wrist.</div>';
    Array.prototype.forEach.call(document.querySelectorAll('#atriumFilters input'), function (el) {
      el.onchange = function () { setLayer(el.dataset.layer, el.checked); };
    });
    // folds to its heading like the side panel, remembered
    var fold = function (min) {
      $('atriumFilters').classList.toggle('min', min); $('filtersMark').textContent = min ? 'show' : 'hide';
      try { localStorage.setItem('genomeatrium.filters', min ? 'min' : 'open'); } catch (e) { /* not kept */ }
    };
    var startMin = false; try { startMin = localStorage.getItem('genomeatrium.filters') === 'min'; } catch (e) { /* open */ }
    fold(startMin);
    $('filtersHead').onclick = function () { fold(!$('atriumFilters').classList.contains('min')); };
  }

  // Atrium (WebXR), created on first use.
  var atrium = null;
  async function enterAtrium() {
    if (!current || !current.genome || !current.genome.contigs.length) return;
    $('xr').hidden = false; document.body.classList.add('atriumMode');
    renderAtriumFilters(); $('atriumFilters').hidden = false;
    view.g.stop && view.g.stop();
    try {
      atrium = G.app.atrium = atrium || new G.Atrium($('xr'));
      await atrium.open(current);
      $('help').innerHTML = (atrium.xrSupported ? 'Enter VR with the button below; ' : '') +
        'drag to orbit, wheel to zoom, hover to read an object, click it to open a window with view tabs (Landscape is a portal into the Landscape room)';
      atriumIntro();
      $('placeBack').hidden = atrium.where === 'atrium'; // back in a room: offer the way out again
    } catch (err) { // no WebGL, or three.js could not load: Arcs instead, and say why
      console.error(err);
      setMode('arcs');
      $('help').innerHTML = '<span class="warn">The Atrium could not open (' + esc(err.message) + '); it needs WebGL and three.js from cdn.jsdelivr.net. Showing Arcs.</span>';
    }
  }
  // A first-visit card on the desktop: what the Atrium shows and how to use it.
  function atriumIntro() {
    var seen = false;
    try { seen = localStorage.getItem('genomeatrium.atriumIntro') === 'seen'; } catch (e) { /* show it */ }
    if (seen || (atrium && atrium.renderer && atrium.renderer.xr.isPresenting)) return;
    var el = $('atriumIntro');
    el.hidden = false;
    var done = function () { el.hidden = true; try { localStorage.setItem('genomeatrium.atriumIntro', 'seen'); } catch (e) { /* not kept */ } };
    $('introOk').onclick = done;
    $('introArcs').onclick = function () { done(); setMode('arcs'); };
  }
  function leaveAtrium() {
    $('atriumIntro').hidden = true;
    $('xr').hidden = true; document.body.classList.remove('atriumMode');
    $('placeBack').hidden = true;
    $('atriumFilters').hidden = true;
    if (atrium) { if (atrium.panel) atrium.closePanel(); atrium.close(); }
    view.g.start && view.g.start();
  }

  // Tracks view (GenomeSpy), created on first use.
  var tracks = null;
  function arcsRange() {
    if (!view.visibleSpan || view.zoom <= 1.01) return null;
    return view.visibleSpan();
  }
  G.app.onTracksInfo = function (t) {
    var gnOn = G.app.gnomad && G.app.gnomad.enabled;
    $('tracksInfo').innerHTML = (t.sampled ? '<span class="warn">variants thinned: every ' + t.sampled + 'th shown, zoom in for all.</span> ' : '') +
      'Variants: ' + (gnOn ? '<span style="color:#ff3c3c">not in gnomAD</span>, <span style="color:#ffa03c">rare</span>, <span style="color:#c8c878">low frequency</span>, <span style="color:#5a6e8c">common</span>'
        : '<span style="color:#78b4ff">SNV</span>, <span style="color:#ffaa50">indel</span>, <span style="color:#eb5ac8">SV</span>') +
      '; faint = failed FILTER. Allele fraction: <span style="color:#be8cff">het</span>, <span style="color:#78b4ff">hom</span>. ' +
      'Enhancers: <span style="color:#ffbe46">distal</span>, <span style="color:#5ac8ff">intragenic</span>, <span style="color:#ff64aa">promoter</span>, thick = both models, red box = rare variant in element. ' +
      (t.span > 20e6 ? 'Genes and enhancers appear below 20 Mb. ' : '') + 'Drawn with GenomeSpy (genomespy.app).';
  };
  async function enterTracks(prev) {
    if (!current || !current.genome || !current.genome.contigs.length) return;
    $('tracks').hidden = false; $('tracksInfo').hidden = false;
    document.body.classList.add('tracksMode');
    view.g.stop && view.g.stop();
    try {
      tracks = G.app.tracks = tracks || new G.Tracks($('tracks'));
      await tracks.open(current, prev === 'arcs' ? arcsRange() : null);
    } catch (err) {
      console.error(err);
      $('tracksInfo').innerHTML = '<span class="warn">' + esc(err.message) + '</span>';
    }
  }
  function leaveTracks(next) {
    $('tracks').hidden = true; $('tracksInfo').hidden = true;
    document.body.classList.remove('tracksMode');
    view.g.start && view.g.start();
    if (next === 'arcs' && tracks && tracks.visible()) { // carry the tracks' region back to Arcs
      var dom = tracks.visible();
      if (dom[0].chrom === dom[1].chrom && tracks.span < current.genome.get(dom[0].chrom).length * 0.9) view.goTo(dom[0].chrom, dom[0].pos, dom[1].pos, { instant: true });
    }
  }

  function setMode(mode) {
    if (mode === 'arena') mode = 'atrium'; // the Atrium's earlier name
    var prev = view.mode;
    document.body.classList.toggle('viewMode', mode !== 'arcs');

    if (prev === 'arcs' && view.visibleSpan && view.zoom > 1.01) lastArcsSpan = view.visibleSpan();
    if (prev === 'tracks' && mode !== 'tracks') leaveTracks(mode);
    if (prev === 'atrium' && mode !== 'atrium') leaveAtrium();
    $('modeAtrium').classList.toggle('on', mode === 'atrium');
    $('modeTracks').classList.toggle('on', mode === 'tracks');
    ['Circos', 'Hilbert', 'Gene', 'Protein'].forEach(function (m) { $('mode' + m).classList.toggle('on', mode === m.toLowerCase()); });
    $('hilbertLayers').hidden = mode !== 'hilbert';
    $('hicBar').hidden = mode !== 'hic';
    $('pathBar').hidden = mode !== 'pathways';
    $('modePathways').classList.toggle('on', mode === 'pathways');
    $('modeMito').classList.toggle('on', mode === 'mito');
    $('modeHic').classList.toggle('on', mode === 'hic');
    prepareView(mode, prev);
    $('modeArcs').classList.toggle('on', mode === 'arcs');
    $('mode3d').classList.toggle('on', mode === '3d');
    $('modeMatrix').classList.toggle('on', mode === 'matrix');
    view.mode = mode;
    $('help').innerHTML = {
      arcs: 'wheel or pinch zoom &middot; drag or swipe pan &middot; shift+drag zoom to &middot; alt+drag select &middot; arrows, +/- &middot; Ctrl+K search &middot; double click resets',
      matrix: 'wheel zoom &middot; drag pan &middot; hover a cell for both windows &middot; click to open in Arcs &middot; double click resets',
      '3d': 'drag to rotate &middot; wheel zoom &middot; hover a window &middot; click it to open that region in Arcs',
      tracks: 'wheel zoom &middot; drag pan &middot; hover for details &middot; switch to Arcs to keep the region',
      hilbert: 'wheel zoom &middot; drag pan &middot; hover a cell &middot; click to open it in Arcs &middot; buttons switch the layer',
      circos: 'drag to rotate &middot; hover chords and chromosomes &middot; click a chromosome to open it in Arcs',
      gene: 'hover an element or a heatmap row &middot; click to open it in Arcs &middot; Ctrl+K picks another gene',
      protein: 'hover lollipops and domains &middot; click a finding in the list, or search a gene, to switch protein',
      atrium: 'loading three.js...',
      mito: 'hover a gene or a variant &middot; click a protein-coding gene for its protein &middot; stems: blue homoplasmic, orange heteroplasmic',
      pathways: 'click a box to open it, the title to go back &middot; the list ranks pathways by enrichment &middot; switches above pick the gene set',
      hic: 'hover contacts and rings &middot; select a region in Arcs (alt+drag) or focus a gene, then load region'
    }[mode];
    if (mode === 'tracks' && prev !== 'tracks') enterTracks(prev);
    if (mode === 'atrium' && prev !== 'atrium') enterAtrium();
  }
  // Setup a view needs before it draws; also used by the Atrium's window tabs.
  function prepareView(mode, prev) {
    if (mode === 'hic' && prev !== 'hic') hicLoad();
    if (mode === 'circos' && current && G.app.circos.data !== current) G.app.circos.setData(current);
  }
  G.app.prepareView = prepareView;
  $('modeArcs').onclick = function () { setMode('arcs'); };
  $('mode3d').onclick = function () { setMode('3d'); };
  $('modeMatrix').onclick = function () { setMode('matrix'); };
  $('modeTracks').onclick = function () { setMode('tracks'); };
  $('modeCircos').onclick = function () { setMode('circos'); };
  $('modeHilbert').onclick = function () { setMode('hilbert'); };
  $('modeGene').onclick = function () { setMode('gene'); };
  $('modeProtein').onclick = function () { setMode('protein'); };
  $('modeHic').onclick = function () { setMode('hic'); };
  $('modeAtrium').onclick = function () { setMode('atrium'); };
  $('modePathways').onclick = function () { setMode('pathways'); };
  $('modeMito').onclick = function () { setMode('mito'); };
  G.app.mitoView = new G.MitoView();
  // Views call this to switch view (a click in Circos opens Arcs, say). Inside
  // the Atrium with its window open, that switches the window's tab instead.
  G.app.setMode = function (mode) {
    if (view.mode === 'atrium' && atrium && atrium.panel && G.Atrium.PANEL_VIEWS.indexOf(mode) >= 0) { atrium.setPanelView(mode); return; }
    setMode(mode);
  };
  view.mode = 'arcs';

  // Dev convenience when served over http: ?url=fixtures/synthetic.vcf.gz
  // Several files load in order: ?url=local/a.g.vcf.gz,local/a.findings.json
  var params = new URLSearchParams(location.search), q = params.get('url');
  // The public HG002 demo (tools/fetch_hg002.py): small variants, SVs, methylation per haplotype.
  // reads for the HG002 demo: GIAB's 300x Illumina BAM (public; read by byte range, CORS open)
  var GIAB_HG002_BAM = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam';
  var DEMO_READS = { hg002: GIAB_HG002_BAM, hg002trio: GIAB_HG002_BAM };
  var DEMO_ROLES = { hg002trio: { hg002: 'child', hg003: 'father', hg004: 'mother' } };
  var DEMO = { hg002trio: ['hg002.wf_snp.vcf.gz', 'hg003.wf_snp.vcf.gz', 'hg004.wf_snp.vcf.gz'].map(function (f) { return 'data/demo/hg002/' + f; }).join(','),
    hg002: ['hg002.wf_snp.vcf.gz', 'hg002.wf_sv.vcf.gz', 'hg002.hap1.methyl.1kb.cov.gz', 'hg002.hap2.methyl.1kb.cov.gz'].map(function (f) { return 'data/demo/hg002/' + f; }).join(',') };
  if (!q && DEMO[params.get('demo')]) q = DEMO[params.get('demo')];
  if (!q) fetch('data/demo/hg002/hg002.wf_snp.vcf.gz', { method: 'HEAD' }).then(function (r) {
    if (r.ok && !current) $('summary').innerHTML += '<br>Or open the <a href="?demo=hg002">public HG002 demo</a> (Genome in a Bottle, nanopore 60x: variants, SVs, methylation per haplotype).';
  }, function () { /* no demo data here */ });
  if (q) {
    var files = q.split(',');
    files.reduce(function (prev, u) {
      return prev.then(function () { return fetch(u); }).then(function (r) { return r.blob(); })
        .then(function (b) { var f = new File([b], u.split('/').pop()); return isPersonFile(f) && current ? addPerson(f) : load(f); });
    }, Promise.resolve()).then(function () {
      var readsUrl = params.get('reads') || DEMO_READS[params.get('demo')];
      if (readsUrl && current) attachReads(G.reads.urlSource(readsUrl), G.reads.urlSource(readsUrl + '.bai'), readsUrl.split('/').pop());
      var roles = DEMO_ROLES[params.get('demo')];
      if (roles) { G.app.people.forEach(function (p) { if (roles[p.name]) p.role = roles[p.name]; }); renderPeople(); }
    });
  }
})(globalThis.G = globalThis.G || {});
