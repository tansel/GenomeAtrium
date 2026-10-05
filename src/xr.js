/*
 * Atrium: the genome as a room, in WebXR (three.js). Works on a desktop with
 * the mouse, and in a headset (Quest browser, Vision Pro Safari) through
 * the "Enter VR" button when the browser reports immersive-vr support.
 *
 *  - The genome is a ring around the viewer at waist height, chromosomes as
 *    coloured arcs, with a variant density (or depth) fringe below it.
 *  - Arches overhead: the strongest similarity links and SV/read-pair links,
 *    the Arcs view turned into a dome.
 *  - Findings float as red spheres with their gene names.
 *  - In the middle, the Landscape: genome windows placed by PCA, threads
 *    per chromosome, at table height.
 * Point at (desktop: hover) a sphere, arch or window to read what it is.
 *
 * The Atrium is the hub. Pulling the trigger on an object (desktop: click)
 * opens a window in front of the viewer with tabs for the page's 2D views
 * (Arcs, Circos, Hilbert, Matrix, Gene, Protein, Hi-C, Pathways): the page's
 * own canvas drawn onto a plane, so each tab is the real view. The
 * Landscape tab turns the window into a portal: a live view into the
 * Landscape room, a room-sized version of the Landscape cloud placed far
 * from the Atrium. The trigger on the portal steps through; in the room the
 * same tab is a portal back.
 * three.js loads from jsDelivr on first use.
 */
(function (G) {
  var THREE_URL = 'https://cdn.jsdelivr.net/npm/three@0.170.0/+esm';
  var ORBIT_URL = 'https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/controls/OrbitControls.js/+esm';
  var VRBTN_URL = 'https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/webxr/VRButton.js/+esm';
  var RING_R = 3.2, RING_Y = 1.0, TAU = Math.PI * 2;
  var MOVE_SPEED = 0.9, TURN_SPEED = 0.7, SLIDE_SPEED = 0.6, SCALE_SPEED = 0.5; // m/s, rad/s (40 deg/s), log-scale/s
  var ROOM_Z = -200, ROOM_START = 3.2;
  var PROT_X = 200, PROT_SCALE = 9; // the Protein room, out of sight too; the protein at 9x (about 5 m across) // the Landscape room, out of sight of the Atrium; viewer starts at its edge

  function hsl(i, n) { return new (G.THREE.Color)().setHSL((i / Math.max(1, n)) * 300 / 360, 0.7, 0.62); }

  function label(THREE, text, color, scale) {
    var c = document.createElement('canvas'), ctx = c.getContext('2d'), fs = 48;
    ctx.font = fs + 'px Helvetica, Arial, sans-serif';
    c.width = Math.ceil(ctx.measureText(text).width) + 24; c.height = fs + 20;
    ctx.font = fs + 'px Helvetica, Arial, sans-serif';
    ctx.fillStyle = 'rgba(0,0,0,0.55)'; ctx.fillRect(0, 0, c.width, c.height);
    ctx.fillStyle = color || 'white'; ctx.textBaseline = 'middle'; ctx.fillText(text, 12, c.height / 2);
    var sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), depthTest: false, transparent: true }));
    var s = scale || 0.12;
    sp.scale.set(s * c.width / c.height, s, 1);
    sp.renderOrder = 10; sp.userData.aspect = c.width / c.height;
    return sp;
  }

  // A flat card with text lines (a Mesh, so controller rays can pick it; sprites cannot).
  function card(THREE, lines, width, color) {
    var c = document.createElement('canvas'), ctx = c.getContext('2d'), fs = 40;
    ctx.font = fs + 'px Helvetica, Arial, sans-serif';
    c.width = Math.ceil(Math.max.apply(null, lines.map(function (t) { return ctx.measureText(t).width; }))) + 48; c.height = lines.length * (fs + 14) + 30;
    ctx.fillStyle = 'rgba(14,14,22,0.9)'; ctx.fillRect(0, 0, c.width, c.height);
    ctx.strokeStyle = color || '#5ad2be'; ctx.lineWidth = 6; ctx.strokeRect(3, 3, c.width - 6, c.height - 6);
    ctx.font = fs + 'px Helvetica, Arial, sans-serif'; ctx.fillStyle = '#fff'; ctx.textBaseline = 'top';
    lines.forEach(function (t, i) { ctx.fillText(t, 24, 18 + i * (fs + 14)); });
    var tex = new THREE.CanvasTexture(c); tex.colorSpace = THREE.SRGBColorSpace;
    return new THREE.Mesh(new THREE.PlaneGeometry(width, width * c.height / c.width), new THREE.MeshBasicMaterial({ map: tex, transparent: true, side: THREE.DoubleSide }));
  }

  function Atrium(el) { this.el = el; this.built = null; }

  Atrium.prototype.open = async function (data) {
    var T = G.THREE = G.THREE || await import(THREE_URL);
    if (!this.renderer) {
      var Orbit = (await import(ORBIT_URL)).OrbitControls, VRButton = (await import(VRBTN_URL)).VRButton;
      var r = this.renderer = new T.WebGLRenderer({ antialias: true });
      r.setPixelRatio(window.devicePixelRatio); r.xr.enabled = true;
      r.xr.setFramebufferScaleFactor(1.25); // sharper text; fixed foveation below pays for it at the edges
      r.xr.setFoveation(1);
      this.el.appendChild(r.domElement);
      this.scene = new T.Scene(); this.scene.background = new T.Color(0x101014);
      this.camera = new T.PerspectiveCamera(70, 1, 0.05, 100);
      this.camera.position.set(0, 1.6, 0.01);
      // The dolly carries the camera and controllers: moving it moves the viewer in VR.
      this.dolly = new T.Group(); this.dolly.add(this.camera); this.scene.add(this.dolly);
      this.clock = new T.Clock();
      this.controls = new Orbit(this.camera, r.domElement);
      this.controls.enableDamping = true;
      this.controls.maxPolarAngle = Math.PI * 0.49; // never below the floor
      this.homeView();
      this.scene.add(new T.AmbientLight(0xffffff, 0.8));
      var dl = new T.DirectionalLight(0xffffff, 0.8); dl.position.set(2, 5, 3); this.scene.add(dl);
      var floor = new T.Mesh(new T.CircleGeometry(RING_R + 1, 64), new T.MeshBasicMaterial({ color: 0x1a1a22 }));
      floor.rotation.x = -Math.PI / 2; this.scene.add(floor);
      this.xrSupported = navigator.xr ? await navigator.xr.isSessionSupported('immersive-vr').catch(function () { return false; }) : false;
      if (this.xrSupported) { var b = VRButton.createButton(r); b.style.position = 'absolute'; this.el.appendChild(b); }
      this.raycaster = new T.Raycaster(); this.raycaster.params.Line.threshold = 0.03; this.raycaster.params.Points = { threshold: 0.03 };
      this.raycaster.far = 40; // never reach from the Atrium into the Landscape room, or back
      this.where = 'atrium';
      this.pointer = new T.Vector2(-9, -9);
      var self = this;
      // Desktop: a click (no drag) does what the trigger does; on the window,
      // dragging pans its view and the wheel zooms it instead of the camera.
      var setPointer = function (e) {
        var rc = r.domElement.getBoundingClientRect();
        self.pointer.set((e.clientX - rc.left) / rc.width * 2 - 1, -(e.clientY - rc.top) / rc.height * 2 + 1);
      };
      this.el.addEventListener('pointerdown', function (e) {
        setPointer(e);
        var hp = self.screenHit(null), onMol = false;
        if (!hp && self.mol) { self.aim(null); var mh = self.mol.hitTest(); onMol = !!mh && mh.object !== self.mol.cardMesh; }
        self.press = { x: e.clientX, y: e.clientY, panel: hp, mol: onMol, lx: e.clientX, ly: e.clientY };
        if (hp || onMol) { self.controls.enabled = false; e.stopPropagation(); }
      }, true);
      this.el.addEventListener('pointermove', function (e) {
        setPointer(e);
        var pr = self.press;
        if (pr && pr.panel) {
          var hp = self.screenHit(null);
          if (hp) { self.panelDrag(hp.x - pr.panel.x, hp.y - pr.panel.y); pr.panel = hp; }
        } else if (pr && pr.mol && self.mol) { // drag turns the protein; shift-drag moves it
          var dx = e.clientX - pr.lx, dy = e.clientY - pr.ly; pr.lx = e.clientX; pr.ly = e.clientY;
          if (e.shiftKey) {
            var T = G.THREE, cq = self.camera.getWorldQuaternion(new T.Quaternion()), k = self.camera.position.distanceTo(self.mol.root.position) * 0.0012;
            self.mol.root.position.addScaledVector(new T.Vector3(1, 0, 0).applyQuaternion(cq), dx * k).addScaledVector(new T.Vector3(0, 1, 0).applyQuaternion(cq), -dy * k);
          } else self.mol.rotateBy(dx * 0.008, dy * 0.008);
        }
      }, true);
      this.el.addEventListener('pointerup', function (e) {
        var pr = self.press; self.press = null; self.controls.enabled = true;
        if (pr && Math.abs(e.clientX - pr.x) + Math.abs(e.clientY - pr.y) < 5 && e.target === r.domElement) self.select(null);
      }, true);
      this.el.addEventListener('wheel', function (e) {
        var hp = self.screenHit(null);
        if (!hp && self.mol) { // wheel over the protein resizes it
          self.aim(null);
          if (self.mol.hitTest()) { e.stopPropagation(); e.preventDefault(); self.mol.root.scale.multiplyScalar(Math.exp(-e.deltaY * 0.0015)); }
          return;
        }
        if (!hp) return;
        e.stopPropagation(); e.preventDefault();
        self.panelZoom(hp, Math.exp(-(e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY) * 0.0018));
      }, { capture: true, passive: false });
      this.setupControllers();
      window.addEventListener('resize', function () { self.resize(); });
      r.setAnimationLoop(function () { self.frame(); });
      r.xr.addEventListener('sessionstart', function () { self.resetPose(); self.showHelp(true); });
      r.xr.addEventListener('sessionend', function () { self.showHelp(false); self.resetPose(); if (self.filterCard) self.toggleFilterCard(); });
    }
    this.resize();
    var me = this;
    G.app.onExtrasChanged = function () { if (me.group) me.buildExtras(); };
    if (this.built !== data) this.build(data);
    this.running = true;
  };

  Atrium.prototype.close = function () { this.running = false; };

  // Desktop camera: above and in front of the ring, the whole of it in view.
  Atrium.prototype.homeView = function () {
    this.camera.position.set(0, 4.4, 6.0);
    this.controls.target.set(0, 0.7, 0);
    this.controls.update();
  };

  Atrium.prototype.resize = function () {
    if (!this.renderer) return;
    var w = this.el.clientWidth, h = this.el.clientHeight;
    this.renderer.setSize(w, h); this.camera.aspect = w / Math.max(1, h); this.camera.updateProjectionMatrix();
  };

  // Controllers: a ray from each; trigger shows what it points at.
  Atrium.prototype.setupControllers = function () {
    var T = G.THREE, self = this;
    this.controllers = [0, 1].map(function (i) {
      var c = self.renderer.xr.getController(i);
      var line = new T.Line(new T.BufferGeometry().setFromPoints([new T.Vector3(0, 0, 0), new T.Vector3(0, 0, -5)]), new T.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.5 }));
      c.add(line);
      c.addEventListener('selectstart', function () { self.select(c); });
      c.addEventListener('selectend', function () { if (self.molDrag && self.molDrag.c === c) self.molDrag = null; });
      c.addEventListener('squeezestart', function () {
        if (self.grabPanel(c)) return; // grip on the window holds the window, not the world
        if (self.mol) { self.aim(c); if (self.mol.aimed()) { c.attach(self.mol.root); self.mol.held = c; return; } }
        c.userData.grip = true; c.userData.last = c.position.clone();
      });
      c.addEventListener('squeezeend', function () {
        c.userData.grip = false; self.releasePanel(c);
        if (self.mol && self.mol.held === c) { self.dolly.attach(self.mol.root); self.mol.held = null; }
      });
      c.addEventListener('connected', function (e) { c.userData.hand = e.data.handedness; });
      self.dolly.add(c);
      return c;
    });
  };

  Atrium.prototype.angle = function (key, pos) {
    var l = this.lay[key];
    return l ? l.a0 + (pos / l.c.length) * (l.a1 - l.a0) : null;
  };
  Atrium.prototype.ringPoint = function (a, y, r) {
    return new G.THREE.Vector3(Math.cos(a) * (r || RING_R), y === undefined ? RING_Y : y, Math.sin(a) * (r || RING_R));
  };

  Atrium.prototype.build = function (d) {
    var T = G.THREE, self = this;
    if (this.group) this.scene.remove(this.group);
    var grp = this.group = new T.Group(); this.scene.add(grp);
    // one group per switchable layer; applyLayers() shows or hides them from the page's switches
    var L = this.layerGroups = {};
    ['ring', 'labels', 'densityHet', 'densityPlain', 'similar', 'arcs', 'findings', 'landscape', 'roh', 'panelRing', 'methyl'].forEach(function (k) { L[k] = new T.Group(); grp.add(L[k]); });
    this.pickables = [];
    this.built = d;
    var cs = d.genome.contigs, total = d.genome.totalLength(), gap = 0.004 * TAU, avail = TAU - gap * cs.length, a = 0, lay = {};
    cs.forEach(function (c, i) { var w = c.length / total * avail; lay[c.key] = { a0: a, a1: a + w, c: c, i: i }; a += w + gap; });
    this.lay = lay;

    // chromosome ring and labels
    cs.forEach(function (c, i) {
      var l = lay[c.key], pts = [];
      for (var k = 0; k <= 40; k++) pts.push(self.ringPoint(l.a0 + (l.a1 - l.a0) * k / 40));
      var tube = new T.Mesh(new T.TubeGeometry(new T.CatmullRomCurve3(pts), 40, 0.025, 6, false), new T.MeshStandardMaterial({ color: hsl(i, cs.length) }));
      tube.userData.info = c.name + '  ' + G.fmtBp(c.length);
      tube.userData.region = { chrom: c.name, start: 1, end: c.length };
      L.ring.add(tube); self.pickables.push(tube);
      if (l.a1 - l.a0 > 0.04) { var lb = label(T, c.name.replace(/^chr/i, ''), '#fff', 0.1); lb.position.copy(self.ringPoint((l.a0 + l.a1) / 2, RING_Y + 0.12, RING_R + 0.05)); L.labels.add(lb); }
    });

    // density fringe below the ring
    var BINS = 720, pos = [], col = [], plain = [];
    var circ = G.app.circos && G.app.circos.rings && G.app.circos.data === d ? G.app.circos.rings : null;
    if (!circ && G.app.circos) { G.app.circos.setData(d); circ = G.app.circos.rings; }
    if (circ) for (var b = 0; b < BINS; b++) {
      if (circ.owner[b] < 0) continue;
      var ang = (b + 0.5) / BINS * TAU, h = Math.sqrt(Math.min(1, circ.dens[b] / circ.hi)) * 0.6;
      var p0 = this.ringPoint(ang, RING_Y - 0.05, RING_R), p1 = this.ringPoint(ang, RING_Y - 0.05 - h, RING_R);
      pos.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
      var cc = circ.het[b] >= 0 ? new T.Color().setRGB((40 + 150 * circ.het[b]) / 255, (30 + 110 * circ.het[b]) / 255, (80 + 175 * circ.het[b]) / 255) : new T.Color(0x78b4ff);
      col.push(cc.r, cc.g, cc.b, cc.r, cc.g, cc.b);
      var pc = new T.Color(d.format === 'vcf' ? 0x78b4ff : 0x5ad2be);
      plain.push(pc.r, pc.g, pc.b, pc.r, pc.g, pc.b);
    }
    if (pos.length) {
      // two copies of the fringe: coloured by het fraction, or plain; the het switch picks one
      [[col, L.densityHet], [plain, L.densityPlain]].forEach(function (cfg) {
        var fg = new T.BufferGeometry(); fg.setAttribute('position', new T.Float32BufferAttribute(pos, 3)); fg.setAttribute('color', new T.Float32BufferAttribute(cfg[0], 3));
        cfg[1].add(new T.LineSegments(fg, new T.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.8 })));
      });
    }

    // arches overhead
    // Arches, merged: one line object for the similarity arches and one mesh for the SV
    // arches, so the headset draws two objects instead of hundreds (each object costs a
    // draw call; too many drop the frame rate, and dropped frames show as dark edges).
    var view = G.app.view, nn = G.genome.normName, SEG = 24, TUBE_SEG = 24, TUBE_SIDES = 5;
    var archCurve = function (a0, a1, height) {
      var p0 = self.ringPoint(a0), p1 = self.ringPoint(a1), mid = p0.clone().add(p1).multiplyScalar(0.5);
      mid.y = RING_Y + height;
      return new T.QuadraticBezierCurve3(p0, mid, p1);
    };
    var simItems = [], sp = [], sc = [];
    (view.simArcs || []).slice(0, 400).forEach(function (s) {
      var a0 = self.angle(nn(s.c0), s.p0), a1 = self.angle(nn(s.c1), s.p1);
      if (a0 === null || a1 === null) return;
      var far = lay[nn(s.c1)], col = hsl(far ? far.i : 0, cs.length), alpha = 0.15 + 0.6 * Math.max(0, (s.support - 0.85) / 0.15);
      var pts = archCurve(a0, a1, 0.4 + 2.2 * Math.abs(Math.sin((a1 - a0) / 2))).getPoints(SEG);
      for (var k = 0; k < SEG; k++) {
        sp.push(pts[k].x, pts[k].y, pts[k].z, pts[k + 1].x, pts[k + 1].y, pts[k + 1].z);
        sc.push(col.r, col.g, col.b, alpha, col.r, col.g, col.b, alpha);
      }
      simItems.push({ info: s.label, region: { chrom: s.c0, start: Math.max(1, s.p0 - 1e6), end: s.p0 + 1e6 } });
    });
    if (simItems.length) {
      var sg = new T.BufferGeometry(); sg.setAttribute('position', new T.Float32BufferAttribute(sp, 3)); sg.setAttribute('color', new T.Float32BufferAttribute(sc, 4));
      var simObj = new T.LineSegments(sg, new T.LineBasicMaterial({ vertexColors: true, transparent: true }));
      simObj.userData.items = simItems; simObj.userData.itemOf = function (hit) { return Math.floor(hit.index / (2 * SEG)); };
      L.similar.add(simObj); self.pickables.push(simObj);
    }
    var svItems = [], tp = [], tn = [], tc = [], ti = [];
    // the 300 largest: links between chromosomes first, then by span (a file's first 300
    // are usually small deletions on chr1)
    var span = function (s) { return nn(s.c0) !== nn(s.c1) ? Infinity : Math.abs(s.p1 - s.p0); };
    (d.arcs || []).filter(function (s) { return s.type !== 'junction'; }).sort(function (x, y) { return span(y) - span(x); }).slice(0, 300).forEach(function (s) {
      var a0 = self.angle(nn(s.c0), s.p0), a1 = self.angle(nn(s.c1), s.p1);
      if (a0 === null || a1 === null) return;
      var same = nn(s.c0) === nn(s.c1), lo = Math.min(s.p0, s.p1), hi = Math.max(s.p0, s.p1), pad = Math.max(200, (hi - lo) * 0.15);
      var tg = new T.TubeGeometry(archCurve(a0, a1, 0.3 + 2.6 * Math.abs(Math.sin((a1 - a0) / 2))), TUBE_SEG, 0.009, TUBE_SIDES, false);
      var base = tp.length / 3, col = new T.Color(G.COLORS.arc[s.type] || '#ccc');
      var P = tg.attributes.position.array, N = tg.attributes.normal.array, I = tg.index.array;
      for (var q = 0; q < P.length; q++) { tp.push(P[q]); tn.push(N[q]); }
      for (q = 0; q < P.length / 3; q++) tc.push(col.r, col.g, col.b);
      for (q = 0; q < I.length; q++) ti.push(I[q] + base);
      tg.dispose();
      svItems.push({ info: s.label, region: same ? { chrom: s.c0, start: Math.max(1, lo - pad), end: hi + pad } : { chrom: s.c0, start: Math.max(1, s.p0 - 50000), end: s.p0 + 50000 } });
    });
    if (svItems.length) {
      var mg = new T.BufferGeometry();
      mg.setAttribute('position', new T.Float32BufferAttribute(tp, 3)); mg.setAttribute('normal', new T.Float32BufferAttribute(tn, 3)); mg.setAttribute('color', new T.Float32BufferAttribute(tc, 3));
      mg.setIndex(ti);
      var svObj = new T.Mesh(mg, new T.MeshStandardMaterial({ vertexColors: true, transparent: true, opacity: 0.9 }));
      var TRIS = TUBE_SEG * TUBE_SIDES * 2;
      svObj.userData.items = svItems; svObj.userData.itemOf = function (hit) { return Math.floor(hit.faceIndex / TRIS); };
      L.arcs.add(svObj); self.pickables.push(svObj);
    }

    // findings (as a skyline), runs of homozygosity and gene panels: buildExtras, also rerun when they change

    // landscape in the middle
    var lm = G.landscape && G.landscape.model;
    if (!lm && G.landscape) { G.landscape.setData(d); lm = G.landscape.model; }
    if (lm) {
      var lp = [], lc = [], S = 0.9 / 260;
      lm.windows.forEach(function (w) { var c = hsl(w.ci, cs.length); lp.push(w.p[0] * S, 1.1 + w.p[1] * S * -1, w.p[2] * S); lc.push(c.r, c.g, c.b); });
      var pg = new T.BufferGeometry(); pg.setAttribute('position', new T.Float32BufferAttribute(lp, 3)); pg.setAttribute('color', new T.Float32BufferAttribute(lc, 3));
      var inst = new T.InstancedMesh(new T.SphereGeometry(0.011, 10, 8), new T.MeshStandardMaterial({ roughness: 0.5 }), lm.windows.length);
      var mtx = new T.Matrix4();
      lm.windows.forEach(function (w, wi) {
        mtx.makeTranslation(lp[3 * wi], lp[3 * wi + 1], lp[3 * wi + 2]); inst.setMatrixAt(wi, mtx);
        inst.setColorAt(wi, hsl(w.ci, cs.length));
      });
      inst.userData.windows = lm.windows; L.landscape.add(inst); self.pickables.push(inst);
      var tp = [];
      for (var i = 1; i < lm.windows.length; i++) {
        var w0 = lm.windows[i - 1], w1 = lm.windows[i];
        if (w0.ci !== w1.ci || w1.j !== w0.j + 1) continue;
        tp.push(lp[3 * (i - 1)], lp[3 * (i - 1) + 1], lp[3 * (i - 1) + 2], lp[3 * i], lp[3 * i + 1], lp[3 * i + 2]);
      }
      var tg = new T.BufferGeometry(); tg.setAttribute('position', new T.Float32BufferAttribute(tp, 3));
      L.landscape.add(new T.LineSegments(tg, new T.LineBasicMaterial({ color: 0x8899aa, transparent: true, opacity: 0.25 })));
    }
    if (this.tip) this.scene.remove(this.tip);
    this.tip = label(T, ' ', '#fff', 0.1); this.tip.visible = false; this.scene.add(this.tip);
    this.buildRoom(d);
    this.buildExtras();
    this.layerSig = null;
    this.applyLayers();
  };

  // Many tubes as one mesh (one draw call); a hit's faceIndex tells which item.
  function mergedTubes(T, specs, tubular, sides) {
    var pos = [], nor = [], col = [], idx = [], items = [];
    specs.forEach(function (sp) {
      var g = new T.TubeGeometry(sp.curve, tubular, sp.radius, sides, false), base = pos.length / 3, c = new T.Color(sp.color);
      var P = g.attributes.position.array, N = g.attributes.normal.array, I = g.index.array;
      for (var q = 0; q < P.length; q++) { pos.push(P[q]); nor.push(N[q]); }
      for (q = 0; q < P.length / 3; q++) col.push(c.r, c.g, c.b);
      for (q = 0; q < I.length; q++) idx.push(I[q] + base);
      g.dispose(); items.push(sp.item);
    });
    var mg = new T.BufferGeometry();
    mg.setAttribute('position', new T.Float32BufferAttribute(pos, 3)); mg.setAttribute('normal', new T.Float32BufferAttribute(nor, 3)); mg.setAttribute('color', new T.Float32BufferAttribute(col, 3));
    mg.setIndex(idx);
    var mesh = new T.Mesh(mg, new T.MeshStandardMaterial({ vertexColors: true, roughness: 0.5 })), TRIS = tubular * sides * 2;
    mesh.userData.items = items; mesh.userData.itemOf = function (hit) { return Math.floor(hit.faceIndex / TRIS); };
    return mesh;
  }

  // Parts of the Atrium that change after the file loads: the findings skyline (heights
  // use the chosen panels and runs of homozygosity), the ROH band and the panel ring.
  Atrium.prototype.buildExtras = function () {
    var T = G.THREE, self = this, L = this.layerGroups, view = G.app.view, nn = G.genome.normName;
    if (!L || !this.lay) return;
    var olds = [];
    ['findings', 'roh', 'panelRing', 'methyl'].forEach(function (k) { L[k].traverse(function (o) { olds.push(o); }); while (L[k].children.length) L[k].remove(L[k].children[0]); });
    this.pickables = this.pickables.filter(function (o) { return olds.indexOf(o) < 0; });

    // Findings skyline. Height adds up stated evidence, it is not a new classification:
    // Asclepius's class (pathogenic 2, likely pathogenic 1.5, other 0.5), +1 when the gene
    // is green in a chosen panel (your phenotype), +0.5 when the panel gene is recessive,
    // the call homozygous, and it lies in a run of homozygosity.
    var green = {};
    (G.app.panelRows || []).forEach(function (r) { if (r.level === 'green') green[r.symbol] = r; });
    (view.findings || []).forEach(function (f) {
      var ang = self.angle(nn(f.chrom), f.pos);
      if (ang === null) return;
      var cls = f.classification || '', why = [];
      var score = /^pathogenic/i.test(cls) ? 2 : /likely pathogenic/i.test(cls) ? 1.5 : 0.5;
      why.push(cls || 'finding');
      var pr = String(f.gene).split(/[;,]/).map(function (g) { return green[g]; }).find(Boolean);
      if (pr) { score += 1; why.push('green in ' + pr.panel); }
      if (pr && pr.biallelic && /hom/i.test(f.zygosity || '') && pr.roh && pr.roh.length) { score += 0.5; why.push('recessive gene, homozygous, in a run of homozygosity'); }
      var h = 0.12 + 0.3 * score, top = self.ringPoint(ang, RING_Y + h, RING_R - 0.1), color = /pathogenic/i.test(cls) && !/likely/i.test(cls) ? 0xff4646 : /pathogenic/i.test(cls) ? 0xff7a3c : 0xffb04a;
      var info = f.gene + ': ' + cls + ', ' + (f.zygosity || '') + (f.variant_name ? '\n' + f.variant_name : '') + '\nheight: ' + why.join(' + ');
      var region = { chrom: f.chrom, start: Math.max(1, f.pos - 60), end: f.pos + 60, gene: f.gene };
      var tower = new T.Mesh(new T.CylinderGeometry(0.014, 0.014, h, 8), new T.MeshStandardMaterial({ color: color, emissive: 0x200808 }));
      tower.position.copy(self.ringPoint(ang, RING_Y + h / 2, RING_R - 0.1));
      var m = new T.Mesh(new T.SphereGeometry(0.045, 16, 12), new T.MeshStandardMaterial({ color: color, emissive: 0x401010 }));
      m.position.copy(top);
      [tower, m].forEach(function (o) { o.userData.info = info; o.userData.region = region; L.findings.add(o); self.pickables.push(o); });
      var lb = label(T, f.gene, '#ffd0d0', 0.08); lb.position.copy(top).add(new T.Vector3(0, 0.1, 0)); L.findings.add(lb);
    });

    // Runs of homozygosity: a gold band just outside the ring (X and Y paler: single copy in a male).
    var roh = G.app.roh, specs = [];
    (roh ? roh.segments : []).forEach(function (sg) {
      var a0 = self.angle(sg.key, sg.start), a1 = self.angle(sg.key, sg.end);
      if (a0 === null || a1 === null) return;
      var pts = [];
      for (var k = 0; k <= 16; k++) pts.push(self.ringPoint(a0 + (a1 - a0) * k / 16, RING_Y - 0.03, RING_R + 0.09));
      specs.push({ curve: new T.CatmullRomCurve3(pts), radius: 0.02, color: G.roh.isAutosome(sg.chrom) ? (sg.long ? 0xffb000 : 0xffd060) : 0x8a7a50,
        item: { info: 'Run of homozygosity ' + sg.chrom + ':' + G.fmtBp(sg.start) + '-' + G.fmtBp(sg.end) + ' (' + (sg.length / 1e6).toFixed(1) + ' Mb)', region: { chrom: sg.chrom, start: sg.start, end: sg.end } } });
    });
    // other people (people.js): their runs in their colour, a band further out each
    var people = G.app.people || [];
    people.forEach(function (p, pi) {
      if (!pi || !p.visible) return;
      var r = p.roh();
      (r ? r.segments : []).forEach(function (sg) {
        var a0 = self.angle(sg.key, sg.start), a1 = self.angle(sg.key, sg.end);
        if (a0 === null || a1 === null) return;
        var pts = [];
        for (var k = 0; k <= 16; k++) pts.push(self.ringPoint(a0 + (a1 - a0) * k / 16, RING_Y - 0.03, RING_R + 0.09 + 0.05 * pi));
        specs.push({ curve: new T.CatmullRomCurve3(pts), radius: 0.016, color: p.color,
          item: { info: p.name + ': run of homozygosity ' + sg.chrom + ':' + G.fmtBp(sg.start) + '-' + G.fmtBp(sg.end) + ' (' + (sg.length / 1e6).toFixed(1) + ' Mb)', region: { chrom: sg.chrom, start: sg.start, end: sg.end } } });
      });
    });
    if (specs.length) { var rm = mergedTubes(T, specs, 16, 6); L.roh.add(rm); this.pickables.push(rm); }

    // other people's findings: thinner towers in their colour, on a smaller ring each
    people.forEach(function (p, pi) {
      if (!pi || !p.visible) return;
      p.findings().forEach(function (f) {
        var ang = self.angle(nn(f.chrom), f.pos);
        if (ang === null) return;
        var cls = f.classification || '', h = 0.12 + 0.3 * (/^pathogenic/i.test(cls) ? 2 : /likely pathogenic/i.test(cls) ? 1.5 : 0.5), rr = RING_R - 0.1 - 0.08 * pi;
        var tower = new T.Mesh(new T.CylinderGeometry(0.009, 0.009, h, 6), new T.MeshStandardMaterial({ color: p.color }));
        tower.position.copy(self.ringPoint(ang, RING_Y + h / 2, rr));
        var m = new T.Mesh(new T.OctahedronGeometry(0.035), new T.MeshStandardMaterial({ color: p.color, emissive: 0x111111 }));
        m.position.copy(self.ringPoint(ang, RING_Y + h, rr));
        var info = p.name + ': ' + f.gene + ', ' + cls + ', ' + (f.zygosity || '') + (f.variant_name ? '\n' + f.variant_name : '');
        [tower, m].forEach(function (o) { o.userData.info = info; o.userData.region = { chrom: f.chrom, start: Math.max(1, f.pos - 60), end: f.pos + 60, gene: f.gene }; L.findings.add(o); self.pickables.push(o); });
      });
    });

    // Gene panel ring inside the main ring: one post per panel gene, coloured by how much
    // of it is callable (green complete, amber gaps, red poor, grey not judged).
    var pc = { complete: 0x4cd07a, gaps: 0xffb43c, poor: 0xff4a4a, unknown: 0x7a8090 }, seen = {};
    specs = [];
    (G.app.panelRows || []).forEach(function (r) {
      if (!r.key || !r.span || seen[r.symbol]) return;
      seen[r.symbol] = 1;
      var ang = self.angle(r.key, (r.span.start + r.span.end) / 2);
      if (ang === null) return;
      var tall = r.level === 'green' ? 0.14 : 0.08;
      specs.push({ curve: new T.LineCurve3(self.ringPoint(ang, RING_Y - 0.06, RING_R - 0.32), self.ringPoint(ang, RING_Y - 0.06 + tall, RING_R - 0.32)), radius: 0.011, color: pc[r.status],
        item: { info: r.symbol + ' (' + r.level + ', ' + (r.moi || 'inheritance not given') + ', ' + r.panel + '): ' + (r.callable === null ? 'coverage not judged' : Math.round(100 * r.callable) + '% callable') +
          (r.roh && r.roh.length ? ', in a run of homozygosity' : '') + (r.findings.length ? ', has a finding' : ''), region: { chrom: r.span.chrom, start: r.span.start, end: r.span.end, gene: r.symbol } } });
    });
    if (specs.length) { var pm = mergedTubes(T, specs, 1, 6); L.panelRing.add(pm); this.pickables.push(pm); }

    // Methylation (when a file is loaded): a ring of short posts outside the fringe, one per
    // 1/720 of the genome, coloured blue (unmethylated) to red (methylated).
    var me = G.app.methyl, d = this.built;
    if (me && d && d.methyl === me) {
      var BINS = 720, total = d.genome.totalLength(), pos = [], col = [], c3 = new T.Color();
      for (var b = 0; b < BINS; b++) {
        var ang = (b + 0.5) / BINS * TAU, hit = null;
        d.genome.contigs.forEach(function (c) { var l = self.lay[c.key]; if (l && ang >= l.a0 && ang < l.a1) hit = l; });
        if (!hit) continue;
        var f0 = (ang - hit.a0) / (hit.a1 - hit.a0), f1 = (ang + TAU / BINS - hit.a0) / (hit.a1 - hit.a0);
        var lv = me.level(hit.c.key, Math.max(1, Math.round(f0 * hit.c.length)), Math.round(Math.min(1, f1) * hit.c.length));
        if (lv.frac === null) continue;
        var p0 = this.ringPoint(ang, RING_Y - 0.72, RING_R + 0.02), p1 = this.ringPoint(ang, RING_Y - 0.72 - 0.12, RING_R + 0.02);
        pos.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
        c3.setStyle(G.methylation.Methylation.color(lv.frac));
        col.push(c3.r, c3.g, c3.b, c3.r, c3.g, c3.b);
      }
      if (pos.length) {
        var mgm = new T.BufferGeometry(); mgm.setAttribute('position', new T.Float32BufferAttribute(pos, 3)); mgm.setAttribute('color', new T.Float32BufferAttribute(col, 3));
        L.methyl.add(new T.LineSegments(mgm, new T.LineBasicMaterial({ vertexColors: true })));
      }
    }
    this.layerSig = null;
  };

  // ----- filters: the page's layer switches (G.app.view.layers) drive the Atrium

  var FILTERS = [
    ['atriumSimilar', 'similarity arches (busy; off by default here)'], ['arcs', 'SV and read-pair arches'], ['findings', 'findings'],
    ['snv', 'SNV density fringe'], ['indel', 'indel density fringe'], ['het', 'colour fringe by het fraction'],
    ['landscape', 'landscape cloud'], ['labels', 'chromosome labels'],
    ['roh', 'runs of homozygosity (gold band)'], ['panelRing', 'gene panel ring (coverage)'], ['methyl', 'methylation ring (if a file is loaded)'],
    ['snapTurn', 'snap turning (comfort)'], ['vignette', 'dim edges when moving (comfort)']
  ];
  Atrium.FILTERS = FILTERS;

  Atrium.prototype.applyLayers = function () {
    var Ly = G.app.view.layers, L = this.layerGroups, on = function (k) { return Ly[k] !== false; };
    var sig = FILTERS.map(function (f) { return on(f[0]) ? 1 : 0; }).join('') + (on('depth') ? 1 : 0);
    if (sig === this.layerSig || !L) return;
    this.layerSig = sig;
    var fringe = this.built && this.built.format === 'vcf' ? (on('snv') || on('indel')) : on('depth');
    L.labels.visible = on('labels');
    L.densityHet.visible = fringe && on('het') && this.built.format === 'vcf';
    L.densityPlain.visible = fringe && !L.densityHet.visible;
    L.similar.visible = Ly.atriumSimilar === true; L.arcs.visible = on('arcs'); // similarity arches: their own switch here, off by default
    L.findings.visible = on('findings'); L.landscape.visible = on('landscape');
    L.roh.visible = on('roh'); L.panelRing.visible = on('panelRing'); L.methyl.visible = on('methyl');
    if (this.filterCard) this.drawFilterCard();
  };

  // Pickable objects that are currently shown (three.js rays also hit hidden ones).
  Atrium.prototype.visiblePickables = function () {
    return this.pickables.filter(function (o) { for (var p = o; p; p = p.parent) if (!p.visible) return false; return true; });
  };

  // VR filter card on the left wrist: left thumbstick click shows or hides it,
  // the right trigger on a line toggles that layer.
  var CARD_W = 512, CARD_ROW = 52;
  Atrium.prototype.toggleFilterCard = function () {
    var T = G.THREE;
    if (this.filterCard) { this.filterCard.mesh.parent.remove(this.filterCard.mesh); this.filterCard = null; return; }
    var left = this.controllers.find(function (c) { return c.userData.hand === 'left'; }) || this.controllers[0];
    var cv = document.createElement('canvas'); cv.width = CARD_W; cv.height = CARD_ROW * (FILTERS.length + 1);
    var tex = new T.CanvasTexture(cv); tex.colorSpace = T.SRGBColorSpace;
    var h = 0.36 * cv.height / cv.width;
    var mesh = new T.Mesh(new T.PlaneGeometry(0.36, h), new T.MeshBasicMaterial({ map: tex, transparent: true, side: T.DoubleSide, depthTest: false }));
    mesh.renderOrder = 30;
    mesh.position.set(0, 0.06 + h / 2, -0.06); mesh.rotation.x = -0.5; // above the wrist, tilted toward the face
    left.add(mesh);
    this.filterCard = { mesh: mesh, tex: tex, cv: cv };
    this.drawFilterCard();
  };

  Atrium.prototype.drawFilterCard = function () {
    var fc = this.filterCard, ctx = fc.cv.getContext('2d'), Ly = G.app.view.layers;
    ctx.fillStyle = 'rgba(14,14,22,0.92)'; ctx.fillRect(0, 0, fc.cv.width, fc.cv.height);
    ctx.fillStyle = '#5ad2be'; ctx.font = 'bold 30px Helvetica, Arial, sans-serif'; ctx.fillText('Filters', 18, 36);
    FILTERS.forEach(function (f, i) {
      var y = CARD_ROW * (i + 1), on = Ly[f[0]] !== false;
      ctx.strokeStyle = '#ddd'; ctx.lineWidth = 3; ctx.strokeRect(18, y + 12, 28, 28);
      if (on) { ctx.fillStyle = '#5ad2be'; ctx.fillRect(23, y + 17, 18, 18); }
      ctx.fillStyle = on ? '#fff' : '#888'; ctx.font = '26px Helvetica, Arial, sans-serif'; ctx.fillText(f[1], 62, y + 36);
    });
    fc.tex.needsUpdate = true;
  };

  // Right trigger on the card: which line, or -1.
  Atrium.prototype.filterCardHit = function (controller) {
    if (!this.filterCard) return -1;
    var T = G.THREE, m = new T.Matrix4().identity().extractRotation(controller.matrixWorld);
    this.raycaster.ray.origin.setFromMatrixPosition(controller.matrixWorld);
    this.raycaster.ray.direction.set(0, 0, -1).applyMatrix4(m);
    var hit = this.raycaster.intersectObject(this.filterCard.mesh, false)[0];
    if (!hit || !hit.uv) return -1;
    var row = Math.floor((1 - hit.uv.y) * this.filterCard.cv.height / CARD_ROW) - 1;
    return row >= 0 && row < FILTERS.length ? row : -1;
  };

  // What a ray hit is: { info, region, action }. Merged objects (arches, the
  // landscape spheres) hold many items; the hit's index says which one.
  Atrium.prototype.pick = function (hit) {
    var u = hit.object.userData;
    if (u.windows) {
      var w = u.windows[hit.instanceId !== undefined ? hit.instanceId : hit.index];
      return { info: w.contig.name + ':' + G.fmtBp(w.start) + '-' + G.fmtBp(w.end), region: { chrom: w.contig.name, start: w.start, end: w.end } };
    }
    if (u.items) return u.items[u.itemOf(hit)] || {};
    return { info: u.info, region: u.region, action: u.action };
  };
  Atrium.prototype.describe = function (hit) { return hit ? this.pick(hit).info || null : null; };

  Atrium.prototype.showTip = function (text, point) {
    var T = G.THREE;
    if (!text) { this.tip.visible = false; return; }
    if (this.tipText !== text) {
      this.scene.remove(this.tip);
      this.tip = label(T, text.split('\n')[0].slice(0, 90), '#fff', 0.09); this.scene.add(this.tip); this.tipText = text;
    }
    this.tip.position.copy(point).add(new T.Vector3(0, 0.12, 0));
    this.tip.visible = true;
  };

  Atrium.prototype.pickFrom = function (controller) {
    var T = G.THREE, m = new T.Matrix4().identity().extractRotation(controller.matrixWorld);
    this.raycaster.ray.origin.setFromMatrixPosition(controller.matrixWorld);
    this.raycaster.ray.direction.set(0, 0, -1).applyMatrix4(m);
    var hit = this.raycaster.intersectObjects(this.visiblePickables(), false)[0];
    this.showTip(this.describe(hit), hit ? hit.point : null);
  };

  // ----- VR navigation (thumbsticks on the xr-standard mapping: axes 2 and 3)
  //  right stick: up/down moves along where you look, left/right turns 30 degrees
  //  left stick: up/down zooms the genome around you, left/right slides sideways
  //  A or X: back to the start; trigger: read what the ray points at
  Atrium.prototype.resetPose = function () {
    this.dolly.rotation.set(0, 0, 0);
    if (this.renderer.xr.isPresenting) this.dolly.position.copy(this.origin(this.where)); else this.dolly.position.set(0, 0, 0);
    if (this.group) { this.group.scale.setScalar(1); this.group.position.set(0, 0, 0); }
    if (this.cloud) { this.cloud.scale.setScalar(1); this.cloud.position.set(0, 1.5, 0); }
    this.turnReady = true; this.vel = null;
  };

  Atrium.prototype.navigateXR = function (dt) {
    var T = G.THREE, session = this.renderer.xr.getSession(), self = this;
    if (!session) return;
    var head = new T.Vector3(), fwd = new T.Vector3();
    this.camera.getWorldPosition(head); this.camera.getWorldDirection(fwd); fwd.y = 0; fwd.normalize();
    var right = new T.Vector3(-fwd.z, 0, fwd.x);
    var onPanel = this.panelHit(); // pointing at the window: sticks drive its view, not the viewer
    var want = { move: 0, turn: 0, slide: 0, scale: 0 }, snap = G.app.view.layers.snapTurn === true;
    var onMol = !onPanel && this.mol && this.controllers.some(function (c) { self.aim(c); return self.mol.aimed(); });
    session.inputSources.forEach(function (src) {
      var gp = src.gamepad;
      if (!gp || gp.axes.length < 4) return;
      var x = gp.axes[2], y = gp.axes[3], dead = 0.15;
      if (onPanel) { // right stick: zoom and pan sideways; left stick: pan in 2D
        if (src.handedness === 'right') {
          if (Math.abs(y) > dead) self.panelZoom(onPanel, Math.exp(-y * 2.2 * dt));
          if (Math.abs(x) > dead) self.panelDrag(-x * 900 * dt, 0);
        } else if (Math.abs(x) > dead || Math.abs(y) > dead) self.panelDrag(Math.abs(x) > dead ? -x * 900 * dt : 0, Math.abs(y) > dead ? -y * 900 * dt : 0);
      } else if (onMol && src.handedness === 'right') { // pointing at the protein: turn it, resize it
        if (Math.abs(x) > dead) { self.mol.spin.rotation.y -= x * 1.2 * dt; self.mol.spinning = false; }
        if (Math.abs(y) > dead) self.mol.root.scale.multiplyScalar(Math.exp(-y * 0.8 * dt));
      } else if (src.handedness === 'right') {
        want.move = Math.abs(y) > dead ? -y * MOVE_SPEED : 0;
        if (snap) { // snap turn: 30 degrees per flick, no turning motion to see
          if (Math.abs(x) > 0.6 && self.turnReady) { self.turnBy(head, x > 0 ? -Math.PI / 6 : Math.PI / 6); self.turnReady = false; }
          if (Math.abs(x) < 0.3) self.turnReady = true;
        } else want.turn = Math.abs(x) > dead ? -x * TURN_SPEED : 0;
      } else if (src.handedness === 'left') {
        want.scale = Math.abs(y) > dead ? -y * SCALE_SPEED : 0;
        want.slide = Math.abs(x) > dead ? x * SLIDE_SPEED : 0;
      }
      if (src.handedness === 'left' && gp.buttons[3]) { // left thumbstick click: filter card
        if (gp.buttons[3].pressed && !self.stickWasDown) self.toggleFilterCard();
        self.stickWasDown = gp.buttons[3].pressed;
      }
      if (gp.buttons[4] && gp.buttons[4].pressed) self.resetPose();          // A or X
      if (gp.buttons[5]) { // B or Y: close the window; with none open, leave the Landscape room
        if (gp.buttons[5].pressed && !self.backWasDown) { if (self.panel) self.closePanel(); else if (self.where !== 'atrium') self.goPlace('atrium'); }
        self.backWasDown = gp.buttons[5].pressed;
      }
    });
    // Ease toward the wanted speeds: sudden starts and stops are what make people sick.
    var v = this.vel = this.vel || { move: 0, turn: 0, slide: 0, scale: 0 }, e = Math.min(1, dt * 5);
    Object.keys(v).forEach(function (k) { v[k] += (want[k] - v[k]) * e; if (Math.abs(v[k]) < 1e-3) v[k] = 0; });
    if (v.move) this.dolly.position.addScaledVector(fwd, v.move * dt);
    if (v.slide) this.dolly.position.addScaledVector(right, v.slide * dt);
    if (v.turn) this.turnBy(head, v.turn * dt);
    if (v.scale) this.scaleAbout(head, Math.exp(v.scale * dt));
    this.motion = Math.abs(v.move) / MOVE_SPEED + Math.abs(v.slide) / SLIDE_SPEED + Math.abs(v.turn) / TURN_SPEED + Math.abs(v.scale) / SCALE_SPEED;
    // grips: one hand drags the world, both hands pull apart or together to scale it
    var gripping = this.controllers.filter(function (c) { return c.userData.grip; });
    if (gripping.length === 2) {
      var a = gripping[0].position, b = gripping[1].position, d = a.distanceTo(b);
      if (this.gripDist) {
        var mid = new T.Vector3().addVectors(a, b).multiplyScalar(0.5);
        this.dolly.localToWorld(mid);
        this.scaleAbout(mid, d / this.gripDist);
      }
      this.gripDist = d;
    } else {
      this.gripDist = null;
      if (gripping.length === 1) {
        var c = gripping[0], delta = c.position.clone().sub(c.userData.last);
        this.dolly.position.sub(delta.applyQuaternion(this.dolly.quaternion));
      }
    }
    this.controllers.forEach(function (c) { if (c.userData.grip) c.userData.last = c.position.clone(); });
  };

  // Turn the viewer by an angle about their head.
  Atrium.prototype.turnBy = function (head, ang) {
    this.dolly.position.sub(head).applyAxisAngle(new G.THREE.Vector3(0, 1, 0), ang).add(head);
    this.dolly.rotation.y += ang;
  };

  // Comfort vignette: darkens the edge of view while moving, which cuts the sense of
  // motion that causes sickness. Off with the 'comfort vignette' switch.
  Atrium.prototype.updateVignette = function () {
    var T = G.THREE, on = G.app.view.layers.vignette !== false && this.renderer.xr.isPresenting, k = on ? Math.min(1, (this.motion || 0) * 1.5) : 0;
    if (!this.vignette) {
      var c = document.createElement('canvas'); c.width = c.height = 256;
      var ctx = c.getContext('2d'), gr = ctx.createRadialGradient(128, 128, 40, 128, 128, 128);
      gr.addColorStop(0, 'rgba(0,0,0,0)'); gr.addColorStop(0.55, 'rgba(0,0,0,0)'); gr.addColorStop(1, 'rgba(0,0,0,1)');
      ctx.fillStyle = gr; ctx.fillRect(0, 0, 256, 256);
      var m = new T.Mesh(new T.PlaneGeometry(1.1, 1.1), new T.MeshBasicMaterial({ map: new T.CanvasTexture(c), transparent: true, depthTest: false, depthWrite: false, opacity: 0 }));
      m.position.z = -0.3; m.renderOrder = 100; m.visible = false;
      this.vignette = m; this.camera.add(m);
    }
    this.vignette.material.opacity = k * 0.85; this.vignette.visible = k > 0.02;
  };

  // Scale the genome (not the viewer) by factor k about a world point.
  Atrium.prototype.scaleAbout = function (p, k) {
    var w = this.where === 'room' ? this.cloud : this.where === 'protein' && this.mol ? this.mol.root : this.group;
    if (!w) return;
    p = w.parent.worldToLocal(p.clone());
    var s0 = w.scale.x, s1 = Math.max(0.1, Math.min(20, s0 * k));
    k = s1 / s0;
    w.position.sub(p).multiplyScalar(k).add(p);
    w.scale.setScalar(s1);
  };

  // ----- the window: the page's 2D views on a plane, with tabs; the Landscape tab is a portal

  // The canvas is laid out at 1000 x 500 CSS px and drawn at 2x: text is larger relative to
  // the window than at 1400 px, and sharp. The window is 2.4 x 1.2 m at 1.6 m.
  var PANEL_W = 1000, PANEL_H = 500, PANEL_DPR = 2, PW = 2.4, PH = 1.2, TAB_H = 80, TAB_CW = 1400;
  var TABS = [['arcs', 'Arcs'], ['circos', 'Circos'], ['hilbert', 'Hilbert'], ['matrix', 'Matrix'], ['gene', 'Gene'],
    ['protein', 'Protein'], ['hic', 'Hi-C'], ['pathways', 'Pathways'], ['mito', 'Mito'], ['3d', 'Landscape'], ['close', 'Close']];
  Atrium.label = label; Atrium.card = card;
  var PANEL_VIEWS = TABS.map(function (t) { return t[0]; }).filter(function (m) { return m !== 'close'; });
  Atrium.PANEL_VIEWS = PANEL_VIEWS;

  // Where the viewer stands in each place.
  Atrium.prototype.origin = function (where) {
    var T = G.THREE;
    return where === 'room' ? new T.Vector3(0, 0, ROOM_Z + ROOM_START) : where === 'protein' ? new T.Vector3(PROT_X, 0, 3.6) : new T.Vector3(0, 0, 0);
  };

  Atrium.prototype.openPanel = function (region) {
    var T = G.THREE, view = G.app.view;
    if (!this.panel) {
      var md = document.getElementById('maindiv');
      md.style.width = PANEL_W + 'px'; md.style.height = PANEL_H + 'px';
      view.forceDpr = PANEL_DPR;
      view.g._adjustCanvas();
      view.g.start && view.g.start();
      view.panelMode = 'arcs'; // opens on Arcs
      var root = new T.Group();
      var tex = new T.CanvasTexture(view.g.canvas);
      tex.colorSpace = T.SRGBColorSpace; tex.anisotropy = this.renderer.capabilities.getMaxAnisotropy(); // crisp at a slant
      var screen = new T.Mesh(new T.PlaneGeometry(PW, PH), new T.MeshBasicMaterial({ map: tex, side: T.DoubleSide }));
      var frame = new T.Mesh(new T.PlaneGeometry(PW + 0.06, PH + 0.06), new T.MeshBasicMaterial({ color: 0x5ad2be, side: T.DoubleSide }));
      frame.position.z = -0.005; screen.add(frame);
      root.add(screen);
      // the window draws over the room (labels ignore depth, and would show through it)
      // (transparent: drawn in the last pass, after the arches, in renderOrder)
      [frame, screen].forEach(function (o) { o.material.depthTest = false; o.material.transparent = true; });
      frame.renderOrder = 15; screen.renderOrder = 16;
      // tab bar above the screen
      var tcv = document.createElement('canvas'); tcv.width = TAB_CW; tcv.height = TAB_H;
      var ttex = new T.CanvasTexture(tcv); ttex.colorSpace = T.SRGBColorSpace;
      var tabs = new T.Mesh(new T.PlaneGeometry(PW, PW * TAB_H / TAB_CW), new T.MeshBasicMaterial({ map: ttex, transparent: true, side: T.DoubleSide }));
      tabs.position.y = PH / 2 + 0.03 + PW * TAB_H / TAB_CW / 2; root.add(tabs);
      tabs.renderOrder = 16; tabs.material.depthTest = false;
      // the portal: a disc showing the other place live, with a glowing rim
      var rt = new T.WebGLRenderTarget(512, 512);
      rt.texture.colorSpace = T.SRGBColorSpace;
      var portal = new T.Mesh(new T.CircleGeometry(0.55, 64), new T.MeshBasicMaterial({ map: rt.texture, side: T.DoubleSide }));
      var rim = new T.Mesh(new T.TorusGeometry(0.56, 0.018, 12, 96), new T.MeshBasicMaterial({ color: 0x5ad2be }));
      portal.add(rim); portal.visible = false; root.add(portal);
      [portal, rim].forEach(function (o) { o.material.depthTest = false; o.material.transparent = true; });
      portal.renderOrder = 16; rim.renderOrder = 17;
      var pcam = new T.PerspectiveCamera(60, 1, 0.05, 60);
      this.panel = { root: root, screen: screen, tex: tex, tabs: tabs, tabCanvas: tcv, tabTex: ttex, portal: portal, rt: rt, pcam: pcam, last: 0, hoverTab: -1 };
      this.dolly.add(root);
    }
    if (this.panel.portalOn) this.setPanelView('arcs');
    // in front of the viewer, facing them
    var head = this.camera.position.clone(), dir = new T.Vector3(0, 0, -1).applyQuaternion(this.camera.quaternion);
    dir.y = 0; dir.normalize();
    this.panel.root.position.copy(head).addScaledVector(dir, 1.6);
    this.panel.root.position.y = Math.max(0.9, Math.min(2.2, head.y));
    this.panel.root.lookAt(this.dolly.localToWorld(head.clone()));
    if (region) {
      view.goTo(region.chrom, region.start, region.end, { instant: true });
      if (region.gene) G.app.focusGene(String(region.gene).split(/[;,]/)[0], false); // Gene, Protein and Hi-C follow it
    }
    this.drawTabs();
    this.showTip(null);
    this.panelActive();
  };

  // Switch the window's tab. '3d' (Landscape) shows the portal instead of the screen.
  Atrium.prototype.setPanelView = function (mode) {
    var p = this.panel, view = G.app.view;
    if (!p) return;
    if (mode === 'close') { this.closePanel(); return; }
    p.portalOn = mode === '3d';
    p.screen.visible = !p.portalOn; p.portal.visible = p.portalOn;
    if (!p.portalOn && view.panelMode !== mode) {
      var prev = view.panelMode;
      view.panelMode = mode;
      G.app.prepareView(mode, prev);
    }
    this.drawTabs();
    this.panelActive();
  };

  Atrium.prototype.drawTabs = function () {
    var p = this.panel, ctx = p.tabCanvas.getContext('2d'), w = TAB_CW / TABS.length, cur = p.portalOn ? '3d' : G.app.view.panelMode, self = this;
    ctx.clearRect(0, 0, TAB_CW, TAB_H);
    TABS.forEach(function (t, i) {
      var on = t[0] === cur, hot = i === p.hoverTab, name = t[0] === '3d' && self.where !== 'atrium' ? 'Atrium' : t[1];
      ctx.fillStyle = on ? '#5ad2be' : hot ? 'rgba(90,210,190,0.35)' : 'rgba(20,20,30,0.9)';
      ctx.fillRect(i * w + 3, 6, w - 6, TAB_H - 10);
      ctx.fillStyle = on ? '#101014' : t[0] === 'close' ? '#ff9a8a' : '#fff';
      ctx.font = (on ? 'bold ' : '') + '25px Helvetica, Arial, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.fillText(name, i * w + w / 2, TAB_H / 2 + 1);
    });
    p.tabTex.needsUpdate = true;
  };

  Atrium.prototype.closePanel = function () {
    if (!this.panel) return;
    this.panel.root.parent.remove(this.panel.root);
    this.panel.rt.dispose();
    this.panel = null;
    var md = document.getElementById('maindiv');
    md.style.width = ''; md.style.height = '';
    var view = G.app.view;
    view.panelMode = null;
    view.forceDpr = null; view.g._adjustCanvas();
    view.g.stop && view.g.stop();
  };

  // Aim the raycaster: a controller's ray, or (src null) the desktop mouse.
  Atrium.prototype.aim = function (src) {
    if (!src) { this.raycaster.setFromCamera(this.pointer, this.camera); return; }
    var m = new G.THREE.Matrix4().identity().extractRotation(src.matrixWorld);
    this.raycaster.ray.origin.setFromMatrixPosition(src.matrixWorld);
    this.raycaster.ray.direction.set(0, 0, -1).applyMatrix4(m);
  };
  Atrium.prototype.sources = function () {
    return this.renderer.xr.isPresenting ? this.controllers : [null];
  };

  // Where one source's ray meets the screen, in canvas pixels, or null.
  Atrium.prototype.screenHit = function (src) {
    if (!this.panel || !this.panel.screen.visible) return null;
    this.aim(src);
    var hit = this.raycaster.intersectObject(this.panel.screen, false)[0];
    return hit && hit.uv ? { x: hit.uv.x * PANEL_W, y: (1 - hit.uv.y) * PANEL_H, distance: hit.distance, src: src } : null;
  };
  // The nearest screen hit over all sources (controllers in VR, the mouse on a desktop).
  Atrium.prototype.panelHit = function () {
    var best = null, self = this;
    if (!this.panel || !this.controllers) return null;
    this.sources().forEach(function (src) { var h = self.screenHit(src); if (h && (!best || h.distance < best.distance)) best = h; });
    return best;
  };
  Atrium.prototype.tabHit = function (src) {
    if (!this.panel) return -1;
    this.aim(src);
    var hit = this.raycaster.intersectObject(this.panel.tabs, false)[0];
    return hit && hit.uv ? Math.min(TABS.length - 1, Math.floor(hit.uv.x * TABS.length)) : -1;
  };
  Atrium.prototype.portalHit = function (src) {
    if (!this.panel || !this.panel.portalOn) return false;
    this.aim(src);
    return !!this.raycaster.intersectObject(this.panel.portal, false)[0];
  };

  // Zoom (factor k about a screen point) and pan (canvas pixels) in the window's current view.
  Atrium.prototype.panelZoom = function (hp, k) {
    var v = G.app.view, m = v.activeMode();
    this.panelActive();
    if (m === 'arcs') { v.flight = null; v.anchorX = hp.x; v.zoomTarget = Math.max(1, v.zoomTarget * k); v.lastMove = Date.now(); v.historyPending = true; }
    else if (m === 'hilbert' && G.app.hilbert) G.app.hilbert.zoomAt(hp.x, hp.y, k);
    else if (m === 'matrix' && G.matrix) G.matrix.zoomAt(hp.x, hp.y, k);
  };
  Atrium.prototype.panelDrag = function (dx, dy) {
    var v = G.app.view, m = v.activeMode();
    this.panelActive();
    if (m === 'arcs') { v.flight = null; v.x0 += dx; v.lastMove = Date.now(); v.historyPending = true; }
    else if (m === 'hilbert' && G.app.hilbert) { G.app.hilbert.ox += dx; G.app.hilbert.oy += dy; }
    else if (m === 'matrix' && G.matrix) { G.matrix.ox += dx; G.matrix.oy += dy; }
    else if (m === 'circos' && G.app.circos) G.app.circos.rot += dx * 0.005;
  };

  // Trigger (or a desktop click): wrist card, tabs, portal, screen, then the world.
  Atrium.prototype.select = function (src) {
    if (this.choiceHit(src)) return;
    if (this.help && !this.helpPending) { this.showHelp(false); this.helpShown = 0; }
    var row = src ? this.filterCardHit(src) : -1;
    if (row >= 0) {
      var key = FILTERS[row][0];
      G.app.setLayer(key, G.app.view.layers[key] === false);
      return;
    }
    if (this.mol) { this.aim(src); if (this.mol.click()) return; }
    var tab = this.tabHit(src);
    if (tab >= 0) { this.setPanelView(TABS[tab][0]); return; }
    if (this.portalHit(src)) { this.goPlace(this.where === 'atrium' ? 'room' : 'atrium'); return; }
    var hp = this.screenHit(src);
    if (hp) {
      var g = G.app.view.g;
      g.mX = hp.x; g.mY = hp.y;
      g.NF_DOWN = g.nF - 1; g.NF_UP = g.nF; // a click on the next Moebio frame
      this.panelActive(3000); // a click may start a load (Protein, Hi-C): keep refreshing a while
      return;
    }
    this.aim(src);
    var hit = this.raycaster.intersectObjects(this.visiblePickables(), false)[0];
    if (!hit) return;
    var it = this.pick(hit);
    if (it.action === 'atrium') { this.goPlace('atrium'); return; }
    if (hit.object.userData.molecule) {
      this.showTip(it.info, hit.point);
      if (this.mol && !this.mol.onPick(it) && src) { var wp = new G.THREE.Vector3(); src.getWorldPosition(wp); this.molDrag = { c: src, last: wp }; } // hold and move to turn it
      return;
    }
    this.showTip(it.info, hit.point);
    if (it.region && it.region.gene) { // a finding or a panel gene: offer the protein too
      var self = this, gene = String(it.region.gene).split(/[;,]/)[0], region = it.region;
      this.askChoice(hit.point, gene, [
        { label: '3D protein, normal and variant side by side', fn: function () { self.showProtein(gene, { focus: true }); } },
        { label: 'Protein room: step inside ' + gene, fn: function () { self.showProtein(gene, { focus: true }).then(function () { self.goPlace('protein'); }); } },
        { label: 'Genome window (Arcs and the other views)', fn: function () { self.openPanel(region); } }]);
      return;
    }
    if (it.region) this.openPanel(it.region);
  };

  // Window upkeep each frame: hover from a controller ray or the mouse, texture refresh.
  Atrium.prototype.updatePanel = function () {
    if (!this.panel) return;
    var hp = this.panelHit(), g = G.app.view.g, self = this, tab = -1;
    if (hp) { g.mX = hp.x; g.mY = hp.y; } else { g.mX = -100; g.mY = -100; }
    this.sources().forEach(function (src) { var t = self.tabHit(src); if (t >= 0) tab = t; });
    if (tab !== this.panel.hoverTab) { this.panel.hoverTab = tab; this.drawTabs(); }
    if (hp) this.panelActive();
    if (G.app.view.panelMode === 'protein') this.followProtein();
    // Uploading the 2000 x 1000 canvas costs frame time: 12 times a second while in use,
    // once a second otherwise (late results such as a protein load still show up).
    var now = performance.now(), every = now < (this.panel.activeUntil || 0) ? 80 : 1000;
    if (this.panel.screen.visible && now - this.panel.last > every) { this.panel.tex.needsUpdate = true; this.panel.last = now; }
  };
  // The Protein tab's protein, in 3D beside the window (AlphaFold model, see molecule.js).
  Atrium.prototype.followProtein = function () {
    var pv = G.app.proteinView, view = G.app.view;
    if (!pv) return;
    var fg = (view.findings || []).map(function (f) { return String(f.gene).split(/[;,]/)[0]; })[0];
    var gene = pv.gene || view.focusGene || fg;
    if (!gene) return;
    pv.load(gene); // UniProt lookup (cached); the Protein view starts it too, but only when it draws
    var p = pv.loaded && pv.loaded[gene];
    if (!p || (this.mol && this.mol.acc === p.acc)) return;
    this.showProtein(gene);
  };

  // The protein of a gene in 3D, held near the viewer; with focus, zoomed to the sample's
  // variant in it (side by side with the normal protein). Resolves once it is built.
  Atrium.prototype.showProtein = async function (gene, opts) {
    opts = opts || {};
    var T = G.THREE, pv = G.app.proteinView;
    gene = String(gene).split(/[;,]/)[0];
    if (!this.mol) { this.mol = new G.Molecule(this); this.dolly.add(this.mol.root); }
    var mol = this.mol;
    if (this.where !== 'protein') this.placeMolecule();
    mol.gene = gene; mol.status = 'Looking up ' + gene + ' in UniProt...'; mol.drawCard();
    var p;
    try { p = await pv.load(gene); } catch (err) { mol.status = err.message; mol.drawCard(); return; }
    if (mol.acc !== p.acc) await mol.load(gene, p);
    if (opts.focus && mol.variants) {
      var i = mol.variants.findIndex(function (v) { return v.source === 'sample'; });
      if (i >= 0) mol.focusVariant(i);
    }
  };

  // Held in front of the viewer, below eye level and to the right, so it does not hide the window.
  Atrium.prototype.placeMolecule = function () {
    var T = G.THREE, root = this.mol.root;
    if (root.parent !== this.dolly) this.dolly.attach(root);
    var head = this.camera.position.clone(), dir = new T.Vector3(0, 0, -1).applyQuaternion(this.camera.quaternion);
    dir.y = 0; dir.normalize();
    var right = new T.Vector3(-dir.z, 0, dir.x);
    root.position.copy(head).addScaledVector(dir, 0.8).addScaledVector(right, 0.35);
    root.position.y = head.y - 0.3;
    root.quaternion.identity(); root.lookAt(this.dolly.localToWorld(head.clone()));
    root.scale.setScalar(1);
    this.mol.cardHome();
  };

  // A small card of choices by a picked object: [{label, fn}], plus Cancel.
  Atrium.prototype.askChoice = function (point, title, options) {
    var T = G.THREE;
    this.closeChoice();
    options = options.concat([{ label: 'Cancel', fn: function () {} }]);
    var W = 900, ROW = 78, c = document.createElement('canvas'); c.width = W; c.height = 80 + ROW * options.length;
    var ctx = c.getContext('2d');
    ctx.fillStyle = 'rgba(14,14,22,0.95)'; ctx.fillRect(0, 0, W, c.height);
    ctx.strokeStyle = '#5ad2be'; ctx.lineWidth = 5; ctx.strokeRect(3, 3, W - 6, c.height - 6);
    ctx.fillStyle = '#5ad2be'; ctx.font = 'bold 38px Helvetica, Arial, sans-serif'; ctx.textBaseline = 'middle'; ctx.fillText(title, 28, 42);
    options.forEach(function (o, i) {
      var y = 80 + i * ROW;
      ctx.fillStyle = o.label === 'Cancel' ? 'rgba(90,40,40,0.9)' : 'rgba(40,44,58,0.95)'; ctx.fillRect(20, y + 6, W - 40, ROW - 12);
      ctx.fillStyle = '#fff'; ctx.font = '34px Helvetica, Arial, sans-serif'; ctx.fillText(o.label, 44, y + ROW / 2);
    });
    var tex = new T.CanvasTexture(c); tex.colorSpace = T.SRGBColorSpace;
    var w = 0.55, m = new T.Mesh(new T.PlaneGeometry(w, w * c.height / W), new T.MeshBasicMaterial({ map: tex, transparent: true, depthTest: false, side: T.DoubleSide }));
    m.renderOrder = 40;
    var cam = this.camera.getWorldPosition(new T.Vector3());
    m.position.copy(point).lerp(cam, 0.35); m.lookAt(cam);
    this.scene.add(m);
    this.choice = { mesh: m, options: options, rowFrac: ROW / c.height, headFrac: 80 / c.height };
  };
  Atrium.prototype.closeChoice = function () { if (this.choice) { this.scene.remove(this.choice.mesh); this.choice = null; } };
  // A trigger or click while a choice card is up: run the option hit, or close the card.
  Atrium.prototype.choiceHit = function (src) {
    if (!this.choice) return false;
    this.aim(src);
    var ch = this.choice, hit = this.raycaster.intersectObject(ch.mesh, false)[0];
    this.closeChoice();
    if (!hit || !hit.uv) return false; // a click elsewhere just closes it
    var row = Math.floor((1 - hit.uv.y - ch.headFrac) / ch.rowFrac);
    if (row >= 0 && row < ch.options.length) ch.options[row].fn();
    return true;
  };
  Atrium.prototype.closeMolecule = function () {
    if (!this.mol) return;
    this.mol.clearModel(); this.mol.root.parent.remove(this.mol.root); this.mol = null;
  };

  Atrium.prototype.panelActive = function (ms) {
    if (this.panel) this.panel.activeUntil = Math.max(this.panel.activeUntil || 0, performance.now() + (ms || 1500));
  };

  // Grip while pointing at the window: it follows the hand until the grip is let go.
  Atrium.prototype.grabPanel = function (c) {
    if (!this.panel || (!this.screenHit(c) && this.tabHit(c) < 0 && !this.portalHit(c))) return false;
    c.attach(this.panel.root); this.panel.heldBy = c;
    return true;
  };
  Atrium.prototype.releasePanel = function (c) {
    if (this.panel && this.panel.heldBy === c) { this.dolly.attach(this.panel.root); this.panel.heldBy = null; }
  };

  // The portal: render the other place from a slowly circling camera into the disc's texture.
  Atrium.prototype.renderPortal = function () {
    var p = this.panel;
    if (!p || !p.portalOn) return;
    if ((this.portalTick = (this.portalTick || 0) + 1) % 2) return; // every other frame is enough
    var t = performance.now() / 1000, cam = p.pcam, r = this.renderer;
    if (this.where === 'atrium') { // looking into the Landscape room
      cam.position.set(Math.sin(t * 0.15) * 5.2, 2.2, ROOM_Z + Math.cos(t * 0.15) * 5.2);
      cam.lookAt(0, 1.5, ROOM_Z);
    } else { // looking back at the Atrium from above its ring
      cam.position.set(Math.sin(t * 0.1) * 5.5, 3.6, Math.cos(t * 0.1) * 5.5);
      cam.lookAt(0, 1.0, 0);
    }
    var xr = r.xr.enabled, prev = r.getRenderTarget();
    r.xr.enabled = false; // render the portal view with the plain camera, not the headset's
    r.setRenderTarget(p.rt);
    r.render(this.scene, cam);
    r.xr.enabled = xr;
    r.setRenderTarget(prev);
  };

  // Move the viewer between the Atrium and the Landscape room.
  Atrium.prototype.goPlace = function (where) {
    if (where === this.where) return;
    if (where === 'protein' && !this.mol) return;
    this.closePanel(); this.closeChoice();
    var leaving = this.where;
    this.where = where;
    this.showTip(null);
    if (where === 'protein') this.enterProteinRoom(); else if (leaving === 'protein' && this.mol) this.placeMolecule();
    if (this.renderer.xr.isPresenting) this.resetPose();
    else { // desktop: move the orbit camera and its target
      if (where === 'room') { this.camera.position.set(0, 2.2, ROOM_Z + 6.5); this.controls.target.set(0, 1.5, ROOM_Z); }
      else if (where === 'protein') { this.camera.position.set(PROT_X, 2.0, 6.0); this.controls.target.set(PROT_X, 1.6, 0); }
      else this.homeView();
    }
    var back = document.getElementById('placeBack'); // the desktop's way back from a room
    if (back) back.hidden = where === 'atrium';
    var h = document.getElementById('help');
    if (h && G.app.view.mode === 'atrium') h.textContent = where === 'room'
      ? 'Landscape room: every dot is a genome window, placed by PCA; threads join neighbours on a chromosome. Click a dot to open it. The sign (or B/Y in VR) goes back.'
      : where === 'protein' ? 'Protein room: the protein at room size, normal and variant side by side at a site. Walk around and into it; the card is beside you. B/Y (or the card) goes back.'
      : 'Atrium: drag to orbit, wheel to zoom, click an object to open it in a window with view tabs';
  };

  // Esc on the desktop: close a choice card, else leave a room, else close the window.
  Atrium.prototype.escape = function () {
    if (this.choice) this.closeChoice();
    else if (this.where !== 'atrium') this.goPlace('atrium');
    else if (this.panel) this.closePanel();
  };

  // The Protein room: the protein at 9x in its own place, with a floor; its card stands
  // beside the viewer at normal size.
  Atrium.prototype.enterProteinRoom = function () {
    var T = G.THREE, root = this.mol.root;
    if (!this.protRoom) {
      this.protRoom = new T.Group(); this.protRoom.position.set(PROT_X, 0, 0); this.scene.add(this.protRoom);
      var floor = new T.Mesh(new T.CircleGeometry(8, 64), new T.MeshBasicMaterial({ color: 0x15151d }));
      floor.rotation.x = -Math.PI / 2; this.protRoom.add(floor);
    }
    this.protRoom.add(root);
    root.position.set(0, 1.7, 0); root.quaternion.identity(); root.scale.setScalar(PROT_SCALE);
    this.mol.cardAway(this.protRoom, new T.Vector3(-0.9, 1.05, 2.6));
  };

  // The Landscape room: the Landscape cloud at room size, around the viewer.
  Atrium.prototype.buildRoom = function (d) {
    var T = G.THREE, self = this;
    if (this.room) { this.scene.remove(this.room); this.pickables = this.pickables.filter(function (o) { return !o.userData.inRoom; }); }
    var room = this.room = new T.Group(); room.position.set(0, 0, ROOM_Z); this.scene.add(room);
    var floor = new T.Mesh(new T.CircleGeometry(7, 64), new T.MeshBasicMaterial({ color: 0x15151d }));
    floor.rotation.x = -Math.PI / 2; room.add(floor);
    var lm = G.landscape && G.landscape.model;
    if (!lm) return;
    var cs = d.genome.contigs, n = lm.windows.length;
    var cloud = this.cloud = new T.Group(); cloud.position.set(0, 1.5, 0); room.add(cloud);
    // Room coordinates per PCA axis: log-compressed about the axis median, so the dense core
    // opens up around the viewer, and scaled so 95% of windows fall within the room (3 m
    // sideways, 1.4 m up and down). Order along each axis is kept; distances are not.
    var reach = [3.0, 1.4, 3.0], fit = [0, 1, 2].map(function (k) {
      var v = lm.windows.map(function (w) { return Math.abs(w.p[k]); }).sort(function (x, y) { return x - y; });
      var med = v[Math.floor(0.5 * (v.length - 1))] || 1, q95 = v[Math.floor(0.95 * (v.length - 1))] || med;
      return function (x) { return Math.sign(x) * reach[k] * Math.log1p(Math.abs(x) / med) / Math.log1p(q95 / med); };
    });
    var pos = new Float32Array(3 * n);
    lm.windows.forEach(function (w, i) { pos[3 * i] = fit[0](w.p[0]); pos[3 * i + 1] = -fit[1](w.p[1]); pos[3 * i + 2] = fit[2](w.p[2]); });
    var inst = new T.InstancedMesh(new T.SphereGeometry(0.028, 10, 8), new T.MeshStandardMaterial({ roughness: 0.45 }), n), mtx = new T.Matrix4();
    lm.windows.forEach(function (w, i) { mtx.makeTranslation(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]); inst.setMatrixAt(i, mtx); inst.setColorAt(i, hsl(w.ci, cs.length)); });
    inst.userData.windows = lm.windows; inst.userData.inRoom = true;
    cloud.add(inst); this.pickables.push(inst);
    // threads: neighbouring windows on a chromosome, in its colour
    var tp = [], tc = [], cen = {};
    for (var i = 0; i < n; i++) {
      var w = lm.windows[i], c = cen[w.ci] || (cen[w.ci] = { x: 0, y: 0, z: 0, n: 0 });
      c.x += pos[3 * i]; c.y += pos[3 * i + 1]; c.z += pos[3 * i + 2]; c.n++;
      if (!i) continue;
      var w0 = lm.windows[i - 1];
      if (w0.ci !== w.ci || w.j !== w0.j + 1) continue;
      var col = hsl(w.ci, cs.length);
      tp.push(pos[3 * i - 3], pos[3 * i - 2], pos[3 * i - 1], pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]);
      tc.push(col.r, col.g, col.b, col.r, col.g, col.b);
    }
    var tg = new T.BufferGeometry(); tg.setAttribute('position', new T.Float32BufferAttribute(tp, 3)); tg.setAttribute('color', new T.Float32BufferAttribute(tc, 3));
    cloud.add(new T.LineSegments(tg, new T.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.22 })));
    // chromosome names at the centre of their windows
    Object.keys(cen).forEach(function (ci) {
      var c = cen[ci], lb = label(T, cs[ci].name.replace(/^chr/i, ''), '#' + hsl(+ci, cs.length).getHexString(), 0.16);
      lb.position.set(c.x / c.n, c.y / c.n + 0.12, c.z / c.n); cloud.add(lb);
    });
    // the way back
    var sign = card(T, ['Back to the Atrium', 'trigger here, or B / Y'], 0.7);
    sign.position.set(-0.8, 1.1, ROOM_START - 1.1); sign.rotation.y = 0.35;
    sign.userData.action = 'atrium'; sign.userData.info = 'Back to the Atrium'; sign.userData.inRoom = true;
    room.add(sign); this.pickables.push(sign);
  };

  // A card with the controls, shown when VR starts. It is placed in the world in front
  // of the viewer (head-locked panels feel like a smudge on the lens), and goes after
  // 20 seconds or the first trigger pull.
  Atrium.prototype.showHelp = function (on) {
    var T = G.THREE;
    if (this.help) { if (this.help.parent) this.help.parent.remove(this.help); this.help = null; }
    this.helpPending = !!on;
    if (!on) return;
    var c = document.createElement('canvas'); c.width = 1000; c.height = 470;
    var ctx = c.getContext('2d');
    ctx.fillStyle = 'rgba(10,10,20,0.9)'; ctx.fillRect(0, 0, c.width, c.height);
    ctx.fillStyle = '#fff'; ctx.font = '36px Helvetica, Arial, sans-serif';
    ['Right stick: move and turn.  Left stick: zoom, slide', 'Grip one hand: drag the world.  Both: pull to resize',
      'Trigger on anything: open it in a window', 'Grip on the window: grab it, pull it closer',
      'On the window: sticks zoom and pan, trigger clicks', 'Tabs switch views; Landscape is a portal',
      'B or Y: close the window.  A or X: start over', 'Left stick click: filters and comfort settings'].forEach(function (t, i) { ctx.fillText(t, 26, 50 + i * 54); });
    var tex = new T.CanvasTexture(c); tex.colorSpace = T.SRGBColorSpace;
    this.help = new T.Mesh(new T.PlaneGeometry(1.0, 0.47), new T.MeshBasicMaterial({ map: tex, transparent: true, depthTest: false }));
    this.help.renderOrder = 20;
  };
  Atrium.prototype.updateHelp = function () {
    if (!this.help) return;
    if (this.helpPending && this.renderer.xr.isPresenting) { // first frame with a head pose: place it ahead, below eye level
      var T = G.THREE, head = this.camera.position.clone(), dir = new T.Vector3(0, 0, -1).applyQuaternion(this.camera.quaternion);
      dir.y = 0; dir.normalize();
      this.help.position.copy(head).addScaledVector(dir, 1.4); this.help.position.y = head.y - 0.25;
      this.help.lookAt(this.dolly.localToWorld(head.clone()));
      this.dolly.add(this.help); this.helpPending = false; this.helpShown = performance.now();
    }
    if (this.helpShown && performance.now() - this.helpShown > 20000) { this.showHelp(false); this.helpShown = 0; }
  };

  Atrium.prototype.frame = function () {
    if (!this.running || !this.group) return;
    var dt = Math.min(0.1, this.clock.getDelta());
    this.applyLayers();
    if (this.renderer.xr.isPresenting) this.navigateXR(dt); else this.motion = 0;
    this.updateVignette();
    if (this.mol) {
      if (this.molDrag) { // hand movement turns the protein: 10 cm is about 35 degrees
        var T = G.THREE, wp = new T.Vector3(); this.molDrag.c.getWorldPosition(wp);
        var d = wp.clone().sub(this.molDrag.last), cq = this.camera.getWorldQuaternion(new T.Quaternion());
        this.mol.rotateBy(d.dot(new T.Vector3(1, 0, 0).applyQuaternion(cq)) * 6, -d.dot(new T.Vector3(0, 1, 0).applyQuaternion(cq)) * 6);
        this.molDrag.last = wp;
      }
      this.mol.update(dt);
    }
    this.updateHelp();
    this.updatePanel();
    if (!this.renderer.xr.isPresenting) {
      this.controls.update();
      var hit = this.screenHit(null) ? null : (this.aim(null), this.raycaster.intersectObjects(this.visiblePickables(), false)[0]);
      this.showTip(this.describe(hit), hit ? hit.point : null);
    }
    this.renderPortal();
    this.renderer.render(this.scene, this.camera);
  };

  G.Atrium = Atrium;
})(globalThis.G = globalThis.G || {});
