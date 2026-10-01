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
    this.spin = new T.Group();   // turns about the focus point
    this.inner = new T.Group();  // model coordinates (Angstrom), scaled and shifted so the focus sits at the root
    this.root.add(this.spin); this.spin.add(this.inner);
    this.style = 'cartoon'; this.colorBy = 'plddt'; this.show = 'both';
    this.scale = 0.01; this.targetScale = 0.01;
    this.focus = new T.Vector3(); this.targetFocus = new T.Vector3();
    this.labels = [];
    this.makeCard();
  }

  // ----- loading

  Molecule.prototype.load = async function (gene, prot) {
    this.gene = gene; this.prot = prot; this.acc = prot.acc;
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

  // The sample's findings in this gene, then ClinVar P/LP residues by count.
  Molecule.prototype.collectVariants = function () {
    var gene = this.gene, out = [], seen = {}, pc = G.proteinChange, ONE = G.structure.ONE;
    var add = function (v) { var k = v.pos + v.short; if (seen[k]) return; seen[k] = 1; out.push(v); };
    var altOf = function (short) { var m = /^[A-Z](\d+)([A-Z])$/.exec(short || ''); return m ? ONE[m[2]] : null; };
    (G.app.view.findings || []).forEach(function (f) {
      if (String(f.gene).split(/[;,]/).indexOf(gene) < 0) return;
      var c = pc(f.variant_name);
      if (c) add({ pos: c.pos, kind: c.kind, short: c.short, alt: altOf(c.short), source: 'sample', n: 1, label: c.short + ' (this sample, ' + (f.classification || 'finding') + ')' });
    });
    var cv = G.app.proteinView ? G.app.proteinView.clinvarFor(gene) : [];
    cv.slice().sort(function (a, b) { return b.n - a.n; }).forEach(function (s) {
      s.names.forEach(function (nm) {
        add({ pos: s.pos, kind: s.kind, short: nm, alt: s.kind === 'missense' ? altOf(nm) : null, source: 'clinvar', n: s.n, label: nm + ' (ClinVar P/LP, ' + s.n + ' at this residue)' });
      });
    });
    var bySeq = this.model.bySeq;
    return out.filter(function (v) { return bySeq[v.pos]; });
  };

  Molecule.prototype.clearModel = function () {
    var self = this;
    this.atrium.pickables = this.atrium.pickables.filter(function (o) { return !o.userData.molecule; });
    while (this.inner.children.length) this.inner.remove(this.inner.children[0]);
    this.labels = [];
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
    this.resIndex = {}; R.forEach(function (r, i) { self.resIndex[r.seq] = i; });
    this.items = R.map(function (r) {
      var v = (self.variantsAt || {})[r.seq];
      return { info: self.gene + ' ' + r.name[0] + r.name.slice(1).toLowerCase() + r.seq + '  pLDDT ' + Math.round(r.plddt) + ', ' + { H: 'helix', E: 'strand', C: 'loop' }[r.ss], residue: r.seq };
    });
    this.buildCartoon(); this.buildAtoms(); this.buildMarkers();
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

  Molecule.prototype.applyStyle = function () {
    var a = this.atrium, self = this;
    if (!this.cartoon) return;
    this.cartoon.visible = this.style !== 'atoms';
    this.atoms.visible = this.bondMesh.visible = this.style !== 'cartoon';
    a.pickables = a.pickables.filter(function (o) { return o !== self.cartoon && o !== self.atoms; });
    if (this.cartoon.visible) a.pickables.push(this.cartoon);
    if (this.atoms.visible) a.pickables.push(this.atoms);
  };

  Molecule.prototype.residueColor = function (r) {
    var T = G.THREE;
    if (this.cut && r.seq > this.cut) return new T.Color(0x3a3d44); // lost after a stop or frameshift
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
      if (this.cut && r.seq > this.cut) tmp.set(0x3a3d44); else if (A.el[i] === 'C') tmp.copy(rc[this.atomRes[i]]); else tmp.set(CPK[A.el[i]] || 0xdddddd);
      this.atoms.setColorAt(i, tmp);
    }
    this.atoms.instanceColor.needsUpdate = true;
  };

  // ----- focus: the whole protein, or one variant site

  Molecule.prototype.whole = function () {
    this.cur = -1; this.cut = null;
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
    this.targetScale = FOCUS_SPAN / 14; this.targetFocus.copy(this.pos(r.ca));
    this.cut = v.kind === 'missense' ? null : v.pos;
    this.change = v.kind === 'missense' && v.alt ? G.structure.mutate(this.model, v.pos, v.alt) : null;
    this.recolor(); this.buildSite(v, r); this.applyFocusLook(); this.drawCard();
  };

  Molecule.prototype.clearSite = function () {
    if (this.site) { this.inner.remove(this.site); this.site = null; }
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
    var lb = G.Atrium.label(T, v.short + (ch ? '  ' + ch.ref + ' > ' + ch.alt : v.kind === 'missense' ? '' : '  ' + v.kind + ': chain greyed after ' + v.pos), '#fff', 1);
    lb.position.copy(ca).add(new T.Vector3(0, 5, 0)); lb.userData.base = 0.06; this.labels.push(lb); site.add(lb);
    this.applyShow();
  };

  // At a site the rest of the protein fades and the other markers hide, so the site's
  // atoms stand out; the whole view brings them back.
  Molecule.prototype.applyFocusLook = function () {
    var site = this.cur >= 0;
    [this.cartoon, this.atoms, this.bondMesh].forEach(function (o) {
      if (!o) return;
      o.material.transparent = site; o.material.opacity = site ? 0.22 : 1; o.material.depthWrite = !site; o.material.needsUpdate = true;
    });
    if (this.markers) this.markers.visible = !site;
    var self = this;
    this.labels.forEach(function (l) { if (!self.site || l.parent !== self.site) l.visible = !site; });
  };

  Molecule.prototype.applyShow = function () {
    if (this.normalGroup) this.normalGroup.visible = this.show !== 'variant';
    if (this.variantGroup) this.variantGroup.visible = this.show !== 'normal';
  };

  // ----- the card: what is shown, and buttons

  var CW = 1280, CH = 470, BTN = [
    { row: 0, k: 'style:cartoon', t: 'Cartoon' }, { row: 0, k: 'style:atoms', t: 'Atoms' }, { row: 0, k: 'style:both', t: 'Both' },
    { row: 0, k: 'color:plddt', t: 'Confidence' }, { row: 0, k: 'color:ss', t: 'Structure' }, { row: 0, k: 'color:domain', t: 'Domains' },
    { row: 1, k: 'whole', t: 'Whole protein' }, { row: 1, k: 'prev', t: '< Variant' }, { row: 1, k: 'next', t: 'Variant >' },
    { row: 1, k: 'show:normal', t: 'Normal' }, { row: 1, k: 'show:variant', t: 'Variant' }, { row: 1, k: 'show:both', t: 'Both' }, { row: 1, k: 'close', t: 'Close' }
  ];
  (function layout() {
    var x = [24, 24];
    BTN.forEach(function (b) { b.w = b.t.length * 17 + 34; b.x = x[b.row]; x[b.row] += b.w + 10; b.y = CH - (2 - b.row) * 74 + 6; b.h = 60; });
  })();

  Molecule.prototype.makeCard = function () {
    var T = G.THREE, c = document.createElement('canvas'); c.width = CW; c.height = CH;
    var tex = new T.CanvasTexture(c); tex.colorSpace = T.SRGBColorSpace;
    var m = new T.Mesh(new T.PlaneGeometry(0.62, 0.62 * CH / CW), new T.MeshBasicMaterial({ map: tex, transparent: true, side: T.DoubleSide, depthTest: false }));
    m.position.set(0, -WHOLE_RADIUS - 0.16, 0); m.userData.molecule = true; m.renderOrder = 18; // over the zoomed protein
    this.cardMesh = m; this.cardCanvas = c; this.cardTex = tex;
    this.root.add(m);
  };

  Molecule.prototype.drawCard = function () {
    var ctx = this.cardCanvas.getContext('2d'), self = this;
    ctx.fillStyle = 'rgba(14,14,22,0.93)'; ctx.fillRect(0, 0, CW, CH);
    ctx.strokeStyle = '#5ad2be'; ctx.lineWidth = 4; ctx.strokeRect(2, 2, CW - 4, CH - 4);
    var lines = [];
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
          lines.push('Normal: the model\'s own residue. Variant (magenta): ideal geometry, common rotamer, not a prediction.');
        } else if (v.kind !== 'missense') lines.push(v.kind + ' at ' + v.pos + ': residues after it greyed, ' + Math.round(100 * (1 - v.pos / this.model.residues.length)) + '% of the chain.');
        else lines.push('Missense with no single new residue named; site shown.');
      } else lines.push('Pick a marker, or use the Variant buttons, to zoom to a site and compare normal and variant.');
    }
    ctx.fillStyle = '#5ad2be'; ctx.font = 'bold 34px Helvetica, Arial, sans-serif'; ctx.textBaseline = 'top';
    ctx.fillText(head.length > 64 ? head.slice(0, 62) + '..' : head, 24, 16);
    ctx.fillStyle = '#fff'; ctx.font = '25px Helvetica, Arial, sans-serif';
    lines.forEach(function (t, i) { ctx.fillText(t.length > 96 ? t.slice(0, 94) + '..' : t, 24, 62 + i * 33); });
    BTN.forEach(function (b) {
      var kv = b.k.split(':'), on = kv[1] && self[{ style: 'style', color: 'colorBy', show: 'show' }[kv[0]]] === kv[1];
      ctx.fillStyle = on ? '#5ad2be' : b.k === 'close' ? 'rgba(120,40,40,0.9)' : 'rgba(40,44,58,0.95)'; ctx.fillRect(b.x, b.y, b.w, b.h);
      ctx.fillStyle = on ? '#101014' : '#fff'; ctx.font = (on ? 'bold ' : '') + '25px Helvetica, Arial, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(b.t, b.x + b.w / 2, b.y + b.h / 2); ctx.textAlign = 'left'; ctx.textBaseline = 'top';
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
    else if (kv[0] === 'color') { this.colorBy = kv[1]; this.recolor(); }
    else if (kv[0] === 'show') { this.show = kv[1]; this.applyShow(); }
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
    return this.atrium.raycaster.ray.intersectsSphere(new T.Sphere(c, WHOLE_RADIUS * this.root.getWorldScale(new T.Vector3()).x));
  };

  Molecule.prototype.update = function (dt) {
    var k = Math.min(1, dt * 4), self = this;
    this.scale += (this.targetScale - this.scale) * k;
    this.focus.lerp(this.targetFocus, k);
    this.inner.scale.setScalar(this.scale);
    this.inner.position.copy(this.focus).multiplyScalar(-this.scale);
    if (this.spinning !== false && this.cur < 0 && !this.held) this.spin.rotation.y += dt * 0.12;
    this.labels.forEach(function (l) { var b = l.userData.base / self.scale; l.scale.set(b * l.userData.aspect, b, 1); });
  };

  G.Molecule = Molecule;
})(globalThis.G = globalThis.G || {});
