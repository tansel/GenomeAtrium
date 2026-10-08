/*
 * The protein in the Atrium: the AlphaFold model of the protein shown in the
 * Protein tab, as a cartoon (helices and strands as ribbons, loops as tubes),
 * as atoms (ball and stick), or both. Coloured by AlphaFold confidence
 * (pLDDT), by secondary structure, or by Pfam domain.
 *
 * Variants are marked on the structure: the sample's findings in red,
 * ClinVar P/LP residues in white (missense) or orange (truncating), sized by
 * how many ClinVar entries fall there. Picking one (or Previous/Next on the
 * card) zooms to the site and shows normal against variant:
 *  - missense: the residue's own side chain from the model (normal) and the
 *    new residue's side chain modelled on the same backbone (variant,
 *    magenta), the residues around it, atoms the new side chain would clash
 *    with, and the change in class, size, hydropathy, burial and confidence;
 *  - stop-gain and frameshift: the part of the chain after the site greyed.
 * These are physical facts about the model. There is no pathogenicity call.
 *
 * Side by side (the default at a site): two linked copies, the normal model
 * on the left and the variant on the right (the modelled side chain, or the
 * chain cut at a stop or frameshift). They turn and zoom together. The right
 * copy is the same model with that change only: AlphaFold predicts one
 * structure, of the normal protein, so for a missense change the two differ
 * at that residue alone. "Overlay" puts both on one copy instead.
 *
 * AlphaMissense (Cheng et al. 2023, served by AlphaFold DB with the model):
 * a predicted pathogenicity for every possible missense change. As a colour,
 * each residue's mean over its 19 substitutions; at a site, the score and
 * class of that change. A published prediction, shown as such.
 */
(function (G) {
  var CPK = { C: 0x9a9a9a, N: 0x3050f8, O: 0xff2a2a, S: 0xffd02a };
  var SS_COL = { H: '#e0507a', E: '#f2c14e', C: '#9aa4b0' };
  var DOMAIN_COLORS = ['#6f5bd3', '#4e79a7', '#59a14f', '#f28e2b', '#e15759', '#76b7b2', '#edc948'];
  var SAMPLES = 5, PROFILE = 8, WHOLE_RADIUS = 0.3, FOCUS_SPAN = 0.4; // metres: the whole protein (held near the hands); 14 A around a site

  function Molecule(atrium) {
    var T = G.THREE;
    this.atrium = atrium;
    this.root = new T.Group();   // placed in the room; grabbed by the grip
    // left copy (normal): anchor -> spin (turns about the focus) -> inner (model coordinates
    // in Angstrom, scaled and shifted so the focus sits at the anchor); right copy (variant)
    // the same, linked to the left each frame
    this.anchorL = new T.Group(); this.spin = new T.Group(); this.inner = new T.Group();
    this.anchorR = new T.Group(); this.spinR = new T.Group(); this.innerR = new T.Group();
    this.root.add(this.anchorL); this.anchorL.add(this.spin); this.spin.add(this.inner);
    this.root.add(this.anchorR); this.anchorR.add(this.spinR); this.spinR.add(this.innerR);
    this.anchorR.visible = false;
    this.style = 'cartoon'; this.colorBy = 'plddt'; this.show = 'both'; this.layout = 'side';
    this.scale = 0.01; this.targetScale = 0.01;
    this.focus = new T.Vector3(); this.targetFocus = new T.Vector3();
    this.labels = [];
    this.makeCard();
  }

  // ----- loading

  Molecule.prototype.load = async function (gene, prot, extra) {
    this.gene = gene; this.prot = prot; this.acc = prot.acc; this.extra = extra || {};
    this.status = 'Loading the AlphaFold model of ' + gene + ' (' + prot.acc + ')...';
    this.clearModel(); this.drawCard();
    try {
      var m = await G.structure.fetchModel(prot.acc);
      if (this.acc !== prot.acc) return; // another protein was picked meanwhile
      this.entry = m; this.model = m.model;
      this.variants = this.collectVariants();
      this.build();
      this.status = null; this.cur = -1;
      this.whole();
    } catch (err) {
      if (this.acc !== prot.acc) return;
      this.status = err.message;
    }
    this.drawCard();
  };

  // Kinds the model can show: a changed residue (site) or a chain that ends early (truncation).
  var SHOWABLE = { missense: 1, nonsense: 1, frameshift: 1, inframe: 1, 'stop lost': 1, 'start lost': 1 };
  var WHY = { synonymous: 'synonymous (same amino acid)', intron: 'intronic', UTR: 'in the untranslated ends', 'splice site': 'at a splice site (changes splicing, not a residue the model can show)',
    'coding block change': 'a multi-base change the page does not translate' };

  // The sample's variants in this gene first (findings, then every protein-changing variant
  // placed on the canonical transcript), then ClinVar P/LP residues by count. sampleNote says
  // what of the sample's variants could not be shown, and why.
  Molecule.prototype.collectVariants = function () {
    var gene = this.gene, out = [], seen = {}, pc = G.proteinChange, ONE = G.structure.ONE, bySeq = this.model.bySeq, x = this.extra || {};
    var add = function (v) { var k = v.pos + v.short; if (seen[k]) return; seen[k] = 1; out.push(v); };
    var altOf = function (short) { var m = /^[A-Z](\d+)([A-Z])$/.exec(short || ''); return m ? ONE[m[2]] : null; };
    var notShown = [], nSample = 0;
    (G.app.view.findings || []).forEach(function (f) {
      if (String(f.gene).split(/[;,]/).indexOf(gene) < 0) return;
      nSample++;
      var c = pc(f.variant_name);
      if (!c) { notShown.push(f.variant_name + ': no protein change named'); return; }
      if (!bySeq[c.pos]) { notShown.push(c.short + ': residue ' + c.pos + ' is outside the AlphaFold model'); return; }
      add({ pos: c.pos, kind: c.kind, short: c.short, alt: altOf(c.short), source: 'sample', n: 1, label: c.short + ' (this sample, ' + (f.classification || 'finding') + ')' });
    });
    var pv = G.app.proteinView, sv = pv ? pv.sampleVariants(gene) : [], why = {};
    sv.forEach(function (v) {
      var cq = v.cq, where = v.chrom + ':' + v.pos.toLocaleString() + ' ' + v.ref + '>' + v.alt;
      if (cq.check === 'mismatch') { notShown.push(where + ': the transcript base differs from the VCF REF'); nSample++; return; }
      if (!SHOWABLE[cq.kind]) { why[cq.kind] = (why[cq.kind] || 0) + 1; nSample++; return; }
      nSample++;
      if (!bySeq[cq.residue]) { notShown.push((cq.short || cq.kind) + ': residue ' + cq.residue + ' is outside the AlphaFold model'); return; }
      add({ pos: cq.residue, kind: cq.kind, short: cq.short || cq.kind, alt: cq.kind === 'missense' ? altOf(cq.short) : null, source: 'sample', n: 1,
        label: (cq.hgvs || cq.short || cq.kind) + ' (this sample, ' + v.zyg + ', ' + cq.kind + ')' });
    });
    var nShown = out.length;
    var cv = pv ? pv.clinvarFor(gene) : [];
    cv.slice().sort(function (a, b) { return b.n - a.n; }).forEach(function (s) {
      s.names.forEach(function (nm) {
        add({ pos: s.pos, kind: s.kind, short: nm, alt: s.kind === 'missense' ? altOf(nm) : null, source: 'clinvar', n: s.n, label: nm + ' (ClinVar P/LP, ' + s.n + ' at this residue)' });
      });
    });
    // what to say when no variant copy is shown
    var other = Object.keys(why).map(function (k) { return why[k] + ' ' + (WHY[k] || k); });
    var noTx = !x.cds && G.app.view.data && /^GRCh3[78]$/.test(G.app.view.data.build) ? 'The ' + gene + ' transcript could not be read (' + (x.cdsError || x.txError || 'Ensembl gave no answer') + '), so only findings were placed. ' : '';
    if (nShown) this.sampleNote = notShown.length ? 'Not shown: ' + notShown.join('; ') + '.' : null;
    else if (!nSample) this.sampleNote = noTx + 'This sample has no variant in ' + gene + ': the normal protein only.';
    else this.sampleNote = noTx + 'No variant copy: this sample\'s ' + nSample + ' variant' + (nSample > 1 ? 's' : '') + ' in ' + gene + ' do' + (nSample > 1 ? '' : 'es') +
      ' not change the protein as the model can show it (' + other.concat(notShown).join('; ') + '). Normal protein only.';
    return out.filter(function (v) { return bySeq[v.pos]; });
  };

  // The transcript arrived after the model: place the sample's variants again and rebuild.
  Molecule.prototype.setExtra = function (extra) {
    this.extra = extra || {};
    this.variants = this.collectVariants();
    this.build();
  };

  // The default view after a load: side by side at the sample's first variant, or the whole
  // normal protein with sampleNote saying why there is no variant copy.
  Molecule.prototype.showSample = function () {
    if (!this.model || !this.variants) return;
    var i = this.variants.findIndex(function (v) { return v.source === 'sample'; });
    if (i >= 0) { this.layout = 'side'; this.show = 'both'; this.applyShow(); this.focusVariant(i); }
    else this.whole();
  };

  Molecule.prototype.clearModel = function () {
    var self = this;
    this.atrium.pickables = this.atrium.pickables.filter(function (o) { return !o.userData.molecule; });
    while (this.inner.children.length) this.inner.remove(this.inner.children[0]);
    while (this.innerR.children.length) this.innerR.remove(this.innerR.children[0]);
    this.labels = []; this.am = null; this.amLoading = null;
    ['cartoon', 'atoms', 'bondMesh', 'markers', 'site'].forEach(function (k) { self[k] = null; });
  };

  // ----- building

  Molecule.prototype.build = function () {
    var T = G.THREE, m = this.model, A = m.atoms, R = m.residues, self = this;
    this.clearModel();
    // centre on the CA atoms; size by their 95th percentile distance
    var c = new T.Vector3();
    R.forEach(function (r) { c.x += A.x[r.ca]; c.y += A.y[r.ca]; c.z += A.z[r.ca]; });
    c.multiplyScalar(1 / R.length); this.center = c;
    var d = R.map(function (r) { return Math.hypot(A.x[r.ca] - c.x, A.y[r.ca] - c.y, A.z[r.ca] - c.z); }).sort(function (a, b) { return a - b; });
    this.r95 = d[Math.floor(0.95 * (d.length - 1))] || 20;
    this.rMax = d[d.length - 1] || this.r95; // loose loops can reach well past r95
    this.resIndex = {}; R.forEach(function (r, i) { self.resIndex[r.seq] = i; });
    this.items = R.map(function (r) {
      var v = (self.variantsAt || {})[r.seq];
      return { info: self.gene + ' ' + r.name[0] + r.name.slice(1).toLowerCase() + r.seq + '  pLDDT ' + Math.round(r.plddt) + ', ' + { H: 'helix', E: 'strand', C: 'loop' }[r.ss], residue: r.seq };
    });
    this.buildCartoon(); this.buildAtoms(); this.buildMarkers(); this.buildTwin();
    this.applyStyle(); this.recolor();
  };

  Molecule.prototype.pos = function (i) { var A = this.model.atoms, c = this.center; return new G.THREE.Vector3(A.x[i] - c.x, A.y[i] - c.y, A.z[i] - c.z); };

  // Cartoon: a spline through the CA atoms, with a cross-section that is a flat
  // ribbon in helices and strands (oriented by the peptide C=O) and a round tube elsewhere.
  Molecule.prototype.buildCartoon = function () {
    var T = G.THREE, m = this.model, R = m.residues, n = R.length, self = this;
    if (n < 2) return;
    var pts = R.map(function (r) { return self.pos(r.ca); });
    var guide = R.map(function (r, i) {
      var g = r.o >= 0 && r.c >= 0 ? self.pos(r.o).sub(self.pos(r.c)).normalize() : new T.Vector3(0, 1, 0);
      return g;
    });
    for (var i = 1; i < n; i++) if (guide[i].dot(guide[i - 1]) < 0) guide[i].negate(); // no half-turn flips
    var curve = new T.CatmullRomCurve3(pts, false, 'centripetal');
    var shape = function (i, f) { // width, thickness at residue i, fraction f along it
      var r = R[i], ss = r.ss;
      if (ss === 'E' && (i === n - 1 || R[i + 1].ss !== 'E')) return [2.2 * (1 - f) + 0.25, 0.35]; // arrow head
      return ss === 'H' ? [1.4, 0.25] : ss === 'E' ? [1.6, 0.35] : [0.3, 0.3];
    };
    var N = (n - 1) * SAMPLES + 1, pos = new Float32Array(N * PROFILE * 3), nor = new Float32Array(N * PROFILE * 3), col = new Float32Array(N * PROFILE * 3);
    var resOf = new Int32Array(N), prevW = null;
    for (var s = 0; s < N; s++) {
      var u = s / (N - 1), ri = Math.min(n - 1, Math.floor(s / SAMPLES)), f = (s % SAMPLES) / SAMPLES;
      var p = curve.getPoint(u), t = curve.getTangent(u);
      var gv = guide[ri].clone().lerp(guide[Math.min(n - 1, ri + 1)], f);
      var side = gv.sub(t.clone().multiplyScalar(gv.dot(t))).normalize(), up = new T.Vector3().crossVectors(t, side);
      var wh = shape(ri, f);
      if (prevW) { wh = [prevW[0] + (wh[0] - prevW[0]) * 0.5, prevW[1] + (wh[1] - prevW[1]) * 0.5]; } // ease shape changes
      prevW = wh;
      resOf[s] = Math.min(n - 1, Math.round(s / SAMPLES));
      for (var k = 0; k < PROFILE; k++) {
        var a = k / PROFILE * Math.PI * 2, cx = Math.cos(a) * wh[0], cy = Math.sin(a) * wh[1];
        var nx = Math.cos(a) / wh[0], ny = Math.sin(a) / wh[1], nl = Math.hypot(nx, ny);
        var o = (s * PROFILE + k) * 3;
        pos[o] = p.x + side.x * cx + up.x * cy; pos[o + 1] = p.y + side.y * cx + up.y * cy; pos[o + 2] = p.z + side.z * cx + up.z * cy;
        nor[o] = (side.x * nx + up.x * ny) / nl; nor[o + 1] = (side.y * nx + up.y * ny) / nl; nor[o + 2] = (side.z * nx + up.z * ny) / nl;
      }
    }
    var idx = [];
    for (s = 0; s < N - 1; s++) for (k = 0; k < PROFILE; k++) {
      var a0 = s * PROFILE + k, a1 = s * PROFILE + (k + 1) % PROFILE, b0 = a0 + PROFILE, b1 = a1 + PROFILE;
      idx.push(a0, b0, a1, a1, b0, b1);
    }
    var g = new T.BufferGeometry();
    g.setAttribute('position', new T.BufferAttribute(pos, 3)); g.setAttribute('normal', new T.BufferAttribute(nor, 3)); g.setAttribute('color', new T.BufferAttribute(col, 3));
    g.setIndex(idx);
    var mesh = new T.Mesh(g, new T.MeshStandardMaterial({ vertexColors: true, roughness: 0.45, metalness: 0.05, side: T.DoubleSide }));
    mesh.userData.molecule = true; mesh.userData.items = this.items;
    mesh.userData.itemOf = function (hit) { return resOf[Math.floor(hit.faceIndex / (PROFILE * 2))]; };
    this.cartoon = mesh; this.cartoonRes = resOf; this.inner.add(mesh);
  };

  // Atoms: one instanced sphere per atom, one instanced cylinder per bond.
  Molecule.prototype.buildAtoms = function () {
    var T = G.THREE, A = this.model.atoms, self = this, mtx = new T.Matrix4();
    var atoms = new T.InstancedMesh(new T.IcosahedronGeometry(0.42, 1), new T.MeshStandardMaterial({ roughness: 0.4 }), A.n);
    for (var i = 0; i < A.n; i++) { var p = this.pos(i); mtx.makeTranslation(p.x, p.y, p.z); atoms.setMatrixAt(i, mtx); }
    var atomRes = new Int32Array(A.n);
    for (i = 0; i < A.n; i++) atomRes[i] = this.resIndex[A.seq[i]] || 0;
    atoms.userData.molecule = true; atoms.userData.items = this.items;
    atoms.userData.itemOf = function (hit) { return atomRes[hit.instanceId]; };
    var bl = G.structure.bonds(A), nb = bl.length / 2;
    var bonds = new T.InstancedMesh(new T.CylinderGeometry(0.16, 0.16, 1, 6, 1, true), new T.MeshStandardMaterial({ color: 0xb8b8b8, roughness: 0.5 }), nb);
    var up = new T.Vector3(0, 1, 0), q = new T.Quaternion(), sc = new T.Vector3();
    for (var b = 0; b < nb; b++) {
      var p0 = this.pos(bl[2 * b]), p1 = this.pos(bl[2 * b + 1]), dir = p1.clone().sub(p0), len = dir.length();
      q.setFromUnitVectors(up, dir.normalize()); sc.set(1, len, 1);
      mtx.compose(p0.clone().add(p1).multiplyScalar(0.5), q, sc); bonds.setMatrixAt(b, mtx);
    }
    this.atoms = atoms; this.atomRes = atomRes; this.bondMesh = bonds;
    this.inner.add(atoms); this.inner.add(bonds);
  };

  // The right copy: the same geometry, colours and materials (no extra memory for the
  // model), with its own draw range so a truncated chain can stop at the site.
  Molecule.prototype.buildTwin = function () {
    var T = G.THREE, self = this;
    if (!this.cartoon) return;
    var g0 = this.cartoon.geometry, g = new T.BufferGeometry();
    Object.keys(g0.attributes).forEach(function (k) { g.setAttribute(k, g0.attributes[k]); });
    g.setIndex(g0.index);
    this.cartoonR = new T.Mesh(g, this.cartoon.material);
    this.atomsR = new T.InstancedMesh(this.atoms.geometry, this.atoms.material, this.atoms.count);
    this.atomsR.instanceMatrix = this.atoms.instanceMatrix; this.atomsR.instanceColor = this.atoms.instanceColor;
    this.bondsR = new T.InstancedMesh(this.bondMesh.geometry, this.bondMesh.material, this.bondMesh.count);
    this.bondsR.instanceMatrix = this.bondMesh.instanceMatrix;
    [[this.cartoonR, this.cartoon], [this.atomsR, this.atoms]].forEach(function (pair) { pair[0].userData = pair[1].userData; });
    this.bondPairs = G.structure.bonds(this.model.atoms); // bond i joins atoms 2i, 2i+1
    [this.cartoonR, this.atomsR, this.bondsR].forEach(function (o) { o.frustumCulled = false; self.innerR.add(o); });
  };

  // The right copy ends at residue seq (a stop or frameshift); null shows it whole.
  Molecule.prototype.cutTwin = function (seq) {
    if (!this.cartoonR) return;
    var ri = seq === null ? Infinity : this.resIndex[seq], resOf = this.cartoonRes, A = this.model.atoms;
    var sMax = resOf.length; for (var s0 = 0; s0 < resOf.length; s0++) if (resOf[s0] >= ri) { sMax = s0; break; }
    this.cartoonR.geometry.setDrawRange(0, seq === null ? Infinity : Math.max(0, sMax - 1) * PROFILE * 6);
    var aMax = A.n; for (var a = 0; a < A.n; a++) if (this.atomRes[a] >= ri) { aMax = a; break; }
    this.atomsR.count = aMax;
    var bp = this.bondPairs, bMax = bp.length / 2; for (var b = 0; b < bp.length / 2; b++) if (bp[2 * b] >= aMax) { bMax = b; break; }
    this.bondsR.count = bMax;
  };

  Molecule.prototype.buildMarkers = function () {
    var T = G.THREE, self = this, V = this.variants, A = this.model.atoms;
    if (!V.length) return;
    var per = {}; // one marker per residue: the sample's wins over ClinVar
    V.forEach(function (v, i) { var k = v.pos; if (!per[k] || (v.source === 'sample' && per[k].v.source !== 'sample')) per[k] = { v: v, i: i }; });
    var list = Object.keys(per).map(function (k) { return per[k]; });
    var mk = new T.InstancedMesh(new T.IcosahedronGeometry(1, 2), new T.MeshStandardMaterial({ roughness: 0.3, emissive: 0x222222 }), list.length);
    var mtx = new T.Matrix4(), items = [];
    list.forEach(function (e, j) {
      var v = e.v, r = self.model.bySeq[v.pos], p = self.pos(r.ca), rad = v.source === 'sample' ? 2.2 : 1.0 + 0.25 * Math.sqrt(v.n);
      mtx.makeScale(rad, rad, rad).setPosition(p); mk.setMatrixAt(j, mtx);
      mk.setColorAt(j, new T.Color(v.source === 'sample' ? 0xff3b3b : v.kind === 'missense' ? 0xffffff : 0xffa03c));
      items.push({ info: v.label, variant: e.i });
      if (v.source === 'sample') { var lb = G.Atrium.label(T, v.short, '#ffd0d0', 1); lb.position.copy(p).add(new T.Vector3(0, 3.5, 0)); lb.userData.base = 0.07; self.labels.push(lb); self.inner.add(lb); }
    });
    mk.userData.molecule = true; mk.userData.items = items; mk.userData.itemOf = function (hit) { return hit.instanceId; };
    this.markers = mk; this.inner.add(mk);
    this.atrium.pickables.push(mk);
  };

  // ----- style and colour

  // AlphaMissense: blue (likely benign) through grey to red (likely pathogenic).
  function amColor(x) {
    var T = G.THREE, lo = new T.Color(0x3b6fd8), mid = new T.Color(0x9a9aa2), hi = new T.Color(0xe8423a);
    return x < 0.5 ? lo.clone().lerp(mid, x / 0.5) : mid.clone().lerp(hi, (x - 0.5) / 0.5);
  }
  var AM_CLASS = { LBen: 'likely benign', Amb: 'ambiguous', LPath: 'likely pathogenic' };

  Molecule.prototype.loadAm = function () {
    var self = this;
    if (this.am || this.amLoading || !this.entry || !this.entry.amUrl) return this.amLoading || Promise.resolve(this.am);
    this.amLoading = G.structure.fetchAlphaMissense(this.entry.amUrl).then(function (am) {
      self.am = am; self.amLoading = null; self.recolor(); self.drawCard(); return am;
    }, function (err) { self.amLoading = null; self.amError = err.message; self.drawCard(); });
    return this.amLoading;
  };

  Molecule.prototype.applyStyle = function () {
    var a = this.atrium, self = this;
    if (!this.cartoon) return;
    this.cartoon.visible = this.style !== 'atoms';
    this.atoms.visible = this.bondMesh.visible = this.style !== 'cartoon';
    if (this.cartoonR) { this.cartoonR.visible = this.cartoon.visible; this.atomsR.visible = this.bondsR.visible = this.atoms.visible; }
    a.pickables = a.pickables.filter(function (o) { return o !== self.cartoon && o !== self.atoms && o !== self.cartoonR && o !== self.atomsR; });
    if (this.cartoon.visible) a.pickables.push(this.cartoon, this.cartoonR);
    if (this.atoms.visible) a.pickables.push(this.atoms, this.atomsR);
    if (this.cur >= 0) this.applyFocusLook(); // keeps the side-by-side framing
  };

  Molecule.prototype.residueColor = function (r) {
    var T = G.THREE;
    if (this.cut && this.layout === 'overlay' && r.seq > this.cut) return new T.Color(0x3a3d44); // overlay: lost part greyed (side by side cuts the right copy)
    if (this.colorBy === 'am') {
      var m = this.am && this.am.mean[r.seq];
      return m === undefined || m === null ? new T.Color(0x5a5f6a) : amColor(m);
    }
    if (this.colorBy === 'ss') return new T.Color(SS_COL[r.ss]);
    if (this.colorBy === 'domain') {
      var ds = this.prot.domains || [];
      for (var i = 0; i < ds.length; i++) if (r.seq >= ds[i].start && r.seq <= ds[i].end) return new T.Color(DOMAIN_COLORS[i % DOMAIN_COLORS.length]);
      return new T.Color(0x6b7280);
    }
    return new T.Color(G.structure.plddtColor(r.plddt));
  };

  Molecule.prototype.recolor = function () {
    if (!this.cartoon) return;
    var R = this.model.residues, A = this.model.atoms, self = this, rc = R.map(function (r) { return self.residueColor(r); });
    var col = this.cartoon.geometry.attributes.color, resOf = this.cartoonRes;
    for (var s = 0; s < resOf.length; s++) { var c = rc[resOf[s]]; for (var k = 0; k < PROFILE; k++) col.setXYZ(s * PROFILE + k, c.r, c.g, c.b); }
    col.needsUpdate = true;
    var T = G.THREE, tmp = new T.Color();
    for (var i = 0; i < A.n; i++) { // carbon takes the residue colour, other elements their own (as in PyMOL)
      var r = R[this.atomRes[i]];
      if (this.cut && this.layout === 'overlay' && r.seq > this.cut) tmp.set(0x3a3d44); else if (A.el[i] === 'C') tmp.copy(rc[this.atomRes[i]]); else tmp.set(CPK[A.el[i]] || 0xdddddd);
      this.atoms.setColorAt(i, tmp);
    }
    this.atoms.instanceColor.needsUpdate = true;
  };

  // ----- focus: the whole protein, or one variant site

  Molecule.prototype.whole = function () {
    this.cur = -1; this.cut = null; this.cutTwin(null); this.twinMode = null;
    if (this.cartoon) this.cartoon.geometry.setDrawRange(0, Infinity);
    this.targetScale = WHOLE_RADIUS / this.r95; this.targetFocus.set(0, 0, 0);
    this.clearSite(); this.recolor(); this.applyFocusLook(); this.drawCard();
  };

  Molecule.prototype.step = function (d) {
    if (!this.variants || !this.variants.length) return;
    var n = this.variants.length;
    this.focusVariant(((this.cur < 0 ? (d > 0 ? -1 : 0) : this.cur) + d + n) % n);
  };

  Molecule.prototype.focusVariant = function (i) {
    var v = this.variants[i], r = this.model.bySeq[v.pos];
    if (!r) return;
    this.cur = i; this.spinning = false;
    var local = v.kind === 'missense' || v.kind === 'inframe' || v.kind === 'stop lost';
    this.cut = local ? null : v.pos;
    // side by side: a missense change is compared up close (the local backbone and the
    // residue's neighbourhood in each copy); a truncation as whole proteins, full and cut
    this.twinMode = local ? 'site' : 'trunc';
    if (this.layout === 'side' && this.twinMode === 'trunc') { this.targetScale = WHOLE_RADIUS / this.r95; this.targetFocus.set(0, 0, 0); }
    else { this.targetScale = FOCUS_SPAN / 14; this.targetFocus.copy(this.pos(r.ca)); }
    this.siteSeq = v.pos;
    this.change = v.kind === 'missense' && v.alt ? G.structure.mutate(this.model, v.pos, v.alt) : null;
    this.cutTwin(this.cut);
    if (v.kind === 'missense') this.loadAm(); // the change's AlphaMissense score, for the card
    this.recolor(); this.buildSite(v, r); this.applyFocusLook(); this.drawCard();
  };

  Molecule.prototype.clearSite = function () {
    if (this.site) { this.inner.remove(this.site); this.site = null; }
    if (this.siteR) { this.innerR.remove(this.siteR); this.siteR = null; }
    var self = this;
    ['capL', 'capR'].forEach(function (k) { if (self[k]) { self[k].parent.remove(self[k]); self[k] = null; } });
  };

  // The site: the residue's own side chain (normal), the modelled one (variant),
  // neighbours within 6 A as thin sticks, clashing atoms marked.
  Molecule.prototype.buildSite = function (v, r) {
    var T = G.THREE, A = this.model.atoms, self = this;
    this.clearSite();
    var site = this.site = new T.Group(); this.inner.add(site);
    site.renderOrder = 5;
    var ball = function (p, rad, color, opacity) {
      var m = new T.Mesh(new T.IcosahedronGeometry(rad, 2), new T.MeshStandardMaterial({ color: color, transparent: opacity < 1, opacity: opacity, roughness: 0.35 }));
      m.position.copy(p); return m;
    };
    var stick = function (p0, p1, rad, color, opacity) {
      var dir = p1.clone().sub(p0), m = new T.Mesh(new T.CylinderGeometry(rad, rad, dir.length(), 8, 1, true), new T.MeshStandardMaterial({ color: color, transparent: opacity < 1, opacity: opacity }));
      m.position.copy(p0).add(p1).multiplyScalar(0.5); m.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), dir.normalize()); return m;
    };
    var bonded = function (list) { // pairs within bonding distance in a small atom list
      var out = [];
      for (var i = 0; i < list.length; i++) for (var j = i + 1; j < list.length; j++) if (list[i].p.distanceTo(list[j].p) < 1.95) out.push([list[i], list[j]]);
      return out;
    };
    // neighbours: residues with an atom within 6 A of this residue's CA, as thin grey sticks
    var ca = this.pos(r.ca), nearRes = {};
    for (var i = 0; i < A.n; i++) if (A.seq[i] !== r.seq && this.pos(i).distanceTo(ca) < 7) nearRes[A.seq[i]] = 1;
    var env = new T.Group();
    Object.keys(nearRes).forEach(function (s) {
      var list = self.model.bySeq[s].atoms.map(function (a) { return { p: self.pos(a), el: A.el[a] }; });
      bonded(list).forEach(function (b) { env.add(stick(b[0].p, b[1].p, 0.12, 0x8899aa, 0.55)); });
    });
    site.add(env);
    // normal: the model's own residue, CPK sticks and balls
    var normal = new T.Group(), own = r.atoms.map(function (a) { return { p: self.pos(a), el: A.el[a], name: A.name[a] }; });
    own.forEach(function (o) { normal.add(ball(o.p, 0.45, CPK[o.el] || 0xdddddd, 1)); });
    bonded(own).forEach(function (b) { normal.add(stick(b[0].p, b[1].p, 0.2, 0xdddddd, 1)); });
    site.add(normal); this.normalGroup = normal;
    // variant: the modelled side chain on the same backbone, magenta; clashes red
    var variant = new T.Group(), ch = this.change;
    if (ch) {
      var c = this.center, sc = ch.side.map(function (s) { return { p: new T.Vector3(s.p[0] - c.x, s.p[1] - c.y, s.p[2] - c.z), el: s.el, name: s.name }; });
      var bb = own.filter(function (o) { return o.name === 'CA'; });
      sc.forEach(function (o) { variant.add(ball(o.p, 0.5, 0xe040fb, 0.85)); });
      bonded(bb.concat(sc)).forEach(function (b) { variant.add(stick(b[0].p, b[1].p, 0.22, 0xe040fb, 0.85)); });
      var hit = {};
      ch.clashes.forEach(function (k) { if (hit[k.with]) return; hit[k.with] = 1; variant.add(ball(self.pos(k.with), 0.7, 0xff2020, 0.6)); });
    }
    site.add(variant); this.variantGroup = variant;
    var lb = G.Atrium.label(T, v.short + (ch ? '  ' + ch.ref + ' > ' + ch.alt : v.kind === 'missense' ? '' : '  ' + v.kind + (this.layout === 'side' ? ' at ' + v.pos : ': chain greyed after ' + v.pos)), '#fff', 1);
    lb.position.copy(ca).add(new T.Vector3(0, 5, 0)); lb.userData.base = 0.06; this.labels.push(lb); site.add(lb);
    // the right copy: the same neighbourhood, with the variant (the modelled side chain;
    // for a truncation the cut chain itself is the difference)
    var siteR = this.siteR = new T.Group(); this.innerR.add(siteR);
    siteR.add(env.clone());
    siteR.add(ch ? variant.clone() : normal.clone());
    // captions over each copy, in the room's scale (not the model's)
    var capText = v.kind === 'missense' ? 'Variant ' + v.short + (ch ? ' (modelled)' : '') : 'Variant ' + v.short + ': ends at ' + v.pos;
    var capY = this.twinMode === 'trunc' ? WHOLE_RADIUS + 0.06 : 0.26;
    this.capL = G.Atrium.label(T, 'Normal', '#cfe', 0.03); this.capL.position.set(0, capY, 0); this.anchorL.add(this.capL);
    this.capR = G.Atrium.label(T, capText, '#f6c', 0.03); this.capR.position.set(0, capY, 0); this.anchorR.add(this.capR);
    this.applyShow();
  };

  // At a site the rest of the protein fades and the other markers hide, so the site's
  // atoms stand out; the whole view brings them back.
  // Cartoon samples covering residues seq a to b (sequence numbers), as an index draw range.
  Molecule.prototype.segmentRange = function (a, b) {
    var resOf = this.cartoonRes, lo = this.resIndex[a], hi = this.resIndex[b], s0 = 0, s1 = resOf.length;
    if (lo === undefined) lo = 0; if (hi === undefined) hi = this.model.residues.length - 1;
    for (var s = 0; s < resOf.length; s++) if (resOf[s] >= lo) { s0 = s; break; }
    for (s = s0; s < resOf.length; s++) if (resOf[s] > hi) { s1 = s; break; }
    return [s0 * PROFILE * 6, Math.max(0, s1 - 1 - s0) * PROFILE * 6];
  };

  Molecule.prototype.applyFocusLook = function () {
    var site = this.cur >= 0, twin = this.twinOn(), local = twin && this.twinMode === 'site';
    var fade = site && !twin; // overlay: the rest fades; side by side: nothing fades
    [this.cartoon, this.atoms, this.bondMesh].forEach(function (o) {
      if (!o) return;
      o.material.transparent = fade; o.material.opacity = fade ? 0.22 : 1; o.material.depthWrite = !fade; o.material.needsUpdate = true;
    });
    // up close side by side: only the local backbone (15 residues either side) in each copy,
    // plus the site groups; whole proteins would overlap each other at this scale
    if (this.cartoon) {
      if (local) {
        var R = this.segmentRange(this.siteSeq - 15, this.siteSeq + 15);
        this.cartoon.geometry.setDrawRange(R[0], R[1]); this.cartoonR.geometry.setDrawRange(R[0], R[1]);
      } else { this.cartoon.geometry.setDrawRange(0, Infinity); this.cutTwin(twin ? this.cut : null); }
      var showAtoms = this.style !== 'cartoon' && !local;
      this.atoms.visible = this.bondMesh.visible = showAtoms;
      if (this.atomsR) this.atomsR.visible = this.bondsR.visible = showAtoms;
    }
    if (this.markers) this.markers.visible = !site;
    var self = this;
    this.labels.forEach(function (l) { if (!self.site || l.parent !== self.site) l.visible = !site; });
  };

  Molecule.prototype.applyShow = function () {
    var side = this.layout === 'side';
    if (this.normalGroup) this.normalGroup.visible = side || this.show !== 'variant';
    if (this.variantGroup) this.variantGroup.visible = !side && this.show !== 'normal'; // side by side: the variant is on the right copy
    if (this.capL) this.capL.visible = side;
    this.recolor();
    if (this.cur >= 0) { // the layout changes how the site is framed
      var v = this.variants[this.cur];
      if (side && this.twinMode === 'trunc') { this.targetScale = WHOLE_RADIUS / this.r95; this.targetFocus.set(0, 0, 0); }
      else if (v) { this.targetScale = FOCUS_SPAN / 14; this.targetFocus.copy(this.pos(this.model.bySeq[v.pos].ca)); }
      this.applyFocusLook();
    }
  };
  Molecule.prototype.twinOn = function () { return this.layout === 'side' && this.cur >= 0 && !!this.cartoonR; };
  // Half the gap between the two copies' centres.
  Molecule.prototype.twinSep = function () {
    if (!this.twinOn()) return 0;
    return this.twinMode === 'site' ? 0.3 : this.rMax * WHOLE_RADIUS / this.r95 + 0.06; // whole proteins: never overlapping
  };

  // ----- the card: what is shown, and buttons

  var CW = 1280, CH = 580, BTN = [
    { row: 0, k: 'style:cartoon', t: 'Cartoon' }, { row: 0, k: 'style:atoms', t: 'Atoms' }, { row: 0, k: 'style:both', t: 'Both' },
    { row: 0, k: 'color:plddt', t: 'Confidence' }, { row: 0, k: 'color:ss', t: 'Structure' }, { row: 0, k: 'color:domain', t: 'Domains' }, { row: 0, k: 'color:am', t: 'AlphaMissense' },
    { row: 1, k: 'whole', t: 'Whole protein' }, { row: 1, k: 'prev', t: '< Variant' }, { row: 1, k: 'next', t: 'Variant >' },
    { row: 1, k: 'close', t: 'Close' },
    { row: 2, k: 'layout:side', t: 'Side by side' }, { row: 2, k: 'layout:overlay', t: 'Overlay' },
    { row: 2, k: 'show:normal', t: 'Normal' }, { row: 2, k: 'show:variant', t: 'Variant' }, { row: 2, k: 'show:both', t: 'Both' }, { row: 2, k: 'room', t: 'Back to Protein view' }
  ];
  (function layout() {
    var x = [24, 24, 24];
    BTN.forEach(function (b) { b.w = b.t.length * 17 + 34; b.x = x[b.row]; x[b.row] += b.w + 10; b.y = CH - (3 - b.row) * 74 + 6; b.h = 60; });
  })();

  Molecule.prototype.makeCard = function () {
    var T = G.THREE, c = document.createElement('canvas'); c.width = CW; c.height = CH;
    var tex = new T.CanvasTexture(c); tex.colorSpace = T.SRGBColorSpace;
    var m = new T.Mesh(new T.PlaneGeometry(0.62, 0.62 * CH / CW), new T.MeshBasicMaterial({ map: tex, transparent: true, side: T.DoubleSide, depthTest: false }));
    m.position.set(0, -WHOLE_RADIUS - 0.16, 0); m.userData.molecule = true; m.renderOrder = 18; // over the zoomed protein
    this.cardMesh = m; this.cardCanvas = c; this.cardTex = tex;
    this.root.add(m);
  };

  // In the Protein room the card stands beside the viewer at normal size; cardHome puts it
  // back under the protein.
  Molecule.prototype.cardAway = function (parent, pos) {
    parent.add(this.cardMesh); this.cardMesh.position.copy(pos); this.cardMesh.rotation.set(0, 0.35, 0); this.cardMesh.scale.setScalar(1.4);
  };
  Molecule.prototype.cardHome = function () {
    this.root.add(this.cardMesh); this.cardMesh.position.set(0, -WHOLE_RADIUS - 0.16, 0); this.cardMesh.rotation.set(0, 0, 0); this.cardMesh.scale.setScalar(1);
  };

  Molecule.prototype.drawCard = function () {
    var ctx = this.cardCanvas.getContext('2d'), self = this;
    ctx.fillStyle = 'rgba(14,14,22,0.93)'; ctx.fillRect(0, 0, CW, CH);
    ctx.strokeStyle = '#5ad2be'; ctx.lineWidth = 4; ctx.strokeRect(2, 2, CW - 4, CH - 4);
    var lines = [], noteLines = [];
    var wrap = function (t, n) { var out = [], cur = ''; t.split(' ').forEach(function (w) { if ((cur + ' ' + w).length > n && cur) { out.push(cur); cur = w; } else cur = cur ? cur + ' ' + w : w; }); if (cur) out.push(cur); return out.slice(0, 3); };
    var head = this.gene ? this.gene + (this.prot && this.prot.name ? ': ' + this.prot.name : '') : 'Protein';
    if (this.status) lines.push(this.status);
    else if (this.model) {
      lines.push('AlphaFold ' + this.entry.entryId + ' v' + this.entry.version + ', ' + this.model.residues.length + ' residues. ' +
        (this.variants.length ? this.variants.length + ' variants on it: red this sample, white ClinVar missense, orange truncating.' : 'No variants to mark.'));
      var v = this.cur >= 0 ? this.variants[this.cur] : null, ch = this.change;
      if (v) {
        lines.push((this.cur + 1) + '/' + this.variants.length + '  ' + v.label);
        if (ch) {
          lines.push(ch.ref + ' (' + ch.refClass + ') to ' + ch.alt + ' (' + ch.altClass + '); size ' + (ch.dVolume >= 0 ? '+' : '') + Math.round(ch.dVolume) + ' A3, hydropathy ' + (ch.dHydropathy >= 0 ? '+' : '') + ch.dHydropathy.toFixed(1));
          lines.push('Site ' + ch.buried + ' (' + ch.near + ' atoms within 10 A), ' + { H: 'helix', E: 'strand', C: 'loop' }[ch.ss] + ', pLDDT ' + Math.round(ch.plddt) +
            '. Modelled ' + ch.alt + ': ' + (ch.clashes.length ? ch.clashes.length + ' contacts under 3 A (red)' : 'no contacts under 3 A'));
          var amv = this.am && this.am.byChange[v.short];
          lines.push('AlphaMissense ' + v.short + ': ' + (amv ? amv.score.toFixed(2) + ' (' + AM_CLASS[amv.cls] + '), a published prediction' : this.amLoading ? 'loading...' : this.amError ? 'unavailable' : 'not listed') +
            '. Variant copy: this model with the modelled side chain only; AlphaFold predicts one structure, the normal one.');
        } else if (v.kind === 'inframe' || v.kind === 'stop lost') lines.push(v.kind + ' at ' + v.pos + ': site shown; the model cannot show the changed chain (residues added or removed).');
        else if (v.kind !== 'missense') lines.push(v.kind + ' at ' + v.pos + ': ' + Math.round(100 * (1 - v.pos / this.model.residues.length)) + '% of the chain lost (' + (this.layout === 'side' ? 'the right copy ends there' : 'greyed') + '). A frameshift may add residues first; not shown.');
        else lines.push('Missense with no single new residue named; site shown.');
      } else if (this.variants.length) lines.push('Pick a marker, or use the Variant buttons, to zoom to a site and compare normal and variant.');
      if (this.sampleNote) noteLines = wrap(this.sampleNote, 92);
    }
    ctx.fillStyle = '#5ad2be'; ctx.font = 'bold 34px Helvetica, Arial, sans-serif'; ctx.textBaseline = 'top';
    ctx.fillText(head.length > 64 ? head.slice(0, 62) + '..' : head, 24, 16);
    ctx.fillStyle = '#fff'; ctx.font = '25px Helvetica, Arial, sans-serif';
    lines.forEach(function (t, i) { ctx.fillText(t.length > 96 ? t.slice(0, 94) + '..' : t, 24, 62 + i * 33); });
    ctx.fillStyle = '#ffd27a'; // the sample note: why there is no variant copy, or what was left out
    noteLines.slice(0, Math.max(0, 9 - lines.length)).forEach(function (t, i) { ctx.fillText(t, 24, 62 + (lines.length + i) * 33); });
    BTN.forEach(function (b) {
      var kv = b.k.split(':'), on = kv[1] && self[{ style: 'style', color: 'colorBy', show: 'show', layout: 'layout' }[kv[0]]] === kv[1];
      ctx.fillStyle = on ? '#5ad2be' : b.k === 'close' ? 'rgba(120,40,40,0.9)' : 'rgba(40,44,58,0.95)'; ctx.fillRect(b.x, b.y, b.w, b.h);
      ctx.fillStyle = on ? '#101014' : '#fff'; ctx.font = (on ? 'bold ' : '') + '25px Helvetica, Arial, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(b.k === 'room' ? (self.atrium.where === 'protein' ? 'Back to ' + self.atrium.backLabel() : 'Protein room') : b.t, b.x + b.w / 2, b.y + b.h / 2); ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    });
    this.cardTex.needsUpdate = true;
  };

  // A trigger (or click) on the card: true when it hit a button.
  Molecule.prototype.click = function () {
    var ray = this.atrium.raycaster, hit = ray.intersectObject(this.cardMesh, false)[0];
    if (!hit || !hit.uv) return false;
    var x = hit.uv.x * CW, y = (1 - hit.uv.y) * CH, b = BTN.find(function (b) { return x >= b.x && x <= b.x + b.w && y >= b.y && y <= b.y + b.h; });
    if (!b) return true;
    var kv = b.k.split(':');
    if (kv[0] === 'style') { this.style = kv[1]; this.applyStyle(); }
    else if (kv[0] === 'color') { this.colorBy = kv[1]; if (kv[1] === 'am') this.loadAm(); this.recolor(); }
    else if (kv[0] === 'show') { this.show = kv[1]; this.applyShow(); }
    else if (kv[0] === 'layout') { this.layout = kv[1]; this.applyShow(); }
    else if (b.k === 'room') { if (this.atrium.where === 'protein') this.atrium.back(); else this.atrium.goPlace('protein'); return true; }
    else if (b.k === 'whole') this.whole();
    else if (b.k === 'prev') this.step(-1);
    else if (b.k === 'next') this.step(1);
    else if (b.k === 'close') { this.atrium.closeMolecule(); return true; }
    this.drawCard();
    return true;
  };

  // A pick on the structure: a marker focuses its variant; a residue just reads.
  Molecule.prototype.onPick = function (it) {
    if (it.variant !== undefined) { this.focusVariant(it.variant); return true; }
    return false;
  };

  // Turn the protein by ax about the viewer's up axis and ay about their right axis
  // (radians), the way a hand or the mouse turns a globe.
  Molecule.prototype.rotateBy = function (ax, ay) {
    var T = G.THREE, cam = this.atrium.camera, cq = cam.getWorldQuaternion(new T.Quaternion());
    var up = new T.Vector3(0, 1, 0).applyQuaternion(cq), right = new T.Vector3(1, 0, 0).applyQuaternion(cq);
    var q = new T.Quaternion().setFromAxisAngle(up, ax).multiply(new T.Quaternion().setFromAxisAngle(right, ay));
    var pq = this.spin.parent.getWorldQuaternion(new T.Quaternion()), pqi = pq.clone().invert();
    this.spin.quaternion.premultiply(pqi.multiply(q).multiply(pq));
    this.spinning = false;
  };

  // The molecule's own objects a ray can hit (structure, markers, card).
  Molecule.prototype.hitTest = function () {
    var objs = [this.cardMesh];
    [this.cartoon, this.atoms, this.markers].forEach(function (o) { if (o && o.visible) objs.push(o); });
    return this.atrium.raycaster.intersectObjects(objs, false)[0] || null;
  };

  // Is the ray (already aimed) pointing at the molecule? Tests its bounding sphere only.
  Molecule.prototype.aimed = function () {
    var T = G.THREE, c = new T.Vector3(); this.root.getWorldPosition(c);
    var r = WHOLE_RADIUS + this.twinSep(); // two copies side by side are wider
    return this.atrium.raycaster.ray.intersectsSphere(new T.Sphere(c, r * this.root.getWorldScale(new T.Vector3()).x));
  };

  Molecule.prototype.update = function (dt) {
    var k = Math.min(1, dt * 4), self = this;
    this.scale += (this.targetScale - this.scale) * k;
    this.focus.lerp(this.targetFocus, k);
    this.inner.scale.setScalar(this.scale);
    this.inner.position.copy(this.focus).multiplyScalar(-this.scale);
    var twin = this.twinOn(), sep = this.twinSep();
    this.anchorR.visible = twin;
    this.anchorL.position.x += (-sep - this.anchorL.position.x) * k;
    this.anchorR.position.x = -this.anchorL.position.x;
    this.spinR.quaternion.copy(this.spin.quaternion);
    this.innerR.scale.copy(this.inner.scale); this.innerR.position.copy(this.inner.position);
    if (this.spinning !== false && this.cur < 0 && !this.held) this.spin.rotation.y += dt * 0.12;
    var rs = this.root.scale.x || 1; // the Protein room scales the whole protein; labels keep their size
    this.labels.forEach(function (l) { var b = l.userData.base / self.scale / rs; l.scale.set(b * l.userData.aspect, b, 1); });
  };

  G.Molecule = Molecule;
})(globalThis.G = globalThis.G || {});
