/*
 * Protein structures from the AlphaFold Protein Structure Database
 * (alphafold.ebi.ac.uk): the predicted model for a UniProt accession,
 * read from its mmCIF file.
 *
 *  - atoms: coordinates (Angstrom), element, atom name, residue number;
 *  - residues: name, number, CA/C/O atoms, pLDDT (AlphaFold's per-residue
 *    confidence, 0 to 100, stored in the B-factor column) and secondary
 *    structure (H helix, E strand, C other) from the file's struct_conf
 *    records (DSSP, written by AlphaFold);
 *  - bonds: atom pairs closer than a covalent cutoff, found with a grid.
 *
 * Only the UniProt accession is sent. AlphaFold DB has no full-length model
 * for human proteins over 2,700 residues; those report "no model".
 */
(function (G) {
  var API = 'https://alphafold.ebi.ac.uk/api/prediction/';

  // Rows of one mmCIF loop: { cols: {name: index}, rows: [[...]] }. Values are
  // split on whitespace, with 'quoted strings' kept whole.
  function cifLoop(lines, prefix) {
    var cols = {}, rows = [], i = 0, n = lines.length;
    for (; i < n; i++) if (lines[i].indexOf(prefix) === 0) break;
    if (i === n) return null;
    // a loop has one "_prefix.col" line per column; a single record is "_prefix.col value"
    var single = lines[i].trim().split(/\s+/).length > 1 && lines[i - 1].trim() !== 'loop_';
    if (single) {
      var row = [];
      for (; i < n && lines[i].indexOf(prefix) === 0; i++) {
        var sp = lines[i].trim().split(/\s+/);
        cols[sp[0].slice(prefix.length)] = row.length; row.push(sp.slice(1).join(' '));
      }
      return { cols: cols, rows: [row] };
    }
    for (var c = 0; i < n && lines[i].indexOf(prefix) === 0; i++) cols[lines[i].trim().slice(prefix.length)] = c++;
    for (; i < n; i++) {
      var l = lines[i];
      if (!l || l[0] === '#' || l[0] === '_' || l.indexOf('loop_') === 0) break;
      rows.push(l.indexOf("'") < 0 ? l.trim().split(/\s+/) : l.trim().match(/'[^']*'|\S+/g).map(function (t) { return t.replace(/^'|'$/g, ''); }));
    }
    return { cols: cols, rows: rows };
  }

  function parseCif(text) {
    var lines = text.split(/\r?\n/);
    var at = cifLoop(lines, '_atom_site.');
    if (!at) throw new Error('No atoms in the structure file');
    var C = at.cols, rows = at.rows.filter(function (r) { return r[C.group_PDB] === 'ATOM' && (C.pdbx_PDB_model_num === undefined || r[C.pdbx_PDB_model_num] === '1'); });
    var n = rows.length, x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n), b = new Float32Array(n), seq = new Int32Array(n), el = [], name = [];
    var residues = [], bySeq = {};
    rows.forEach(function (r, i) {
      x[i] = +r[C.Cartn_x]; y[i] = +r[C.Cartn_y]; z[i] = +r[C.Cartn_z]; b[i] = +r[C.B_iso_or_equiv];
      seq[i] = +r[C.label_seq_id]; el[i] = r[C.type_symbol]; name[i] = r[C.label_atom_id];
      var res = bySeq[seq[i]];
      if (!res) { res = bySeq[seq[i]] = { seq: seq[i], name: r[C.label_comp_id], plddt: b[i], ss: 'C', ca: -1, c: -1, o: -1, n: -1, atoms: [] }; residues.push(res); }
      res.atoms.push(i);
      if (name[i] === 'CA') res.ca = i; else if (name[i] === 'C') res.c = i; else if (name[i] === 'O') res.o = i; else if (name[i] === 'N') res.n = i;
    });
    var sc = cifLoop(lines, '_struct_conf.');
    if (sc) sc.rows.forEach(function (r) {
      var t = r[sc.cols.conf_type_id] || '', ss = /^HELX_RH/.test(t) ? 'H' : /^STRN/.test(t) ? 'E' : null;
      if (!ss) return;
      for (var s = +r[sc.cols.beg_label_seq_id]; s <= +r[sc.cols.end_label_seq_id]; s++) if (bySeq[s]) bySeq[s].ss = ss;
    });
    residues = residues.filter(function (r) { return r.ca >= 0; });
    return { atoms: { n: n, x: x, y: y, z: z, b: b, seq: seq, el: el, name: name }, residues: residues, bySeq: bySeq };
  }

  // Covalent bonds: pairs within 1.9 A (2.1 A when sulfur is involved), via a 2.2 A grid.
  function bonds(atoms) {
    var cell = 2.2, grid = new Map(), out = [], A = atoms;
    var key = function (i, j, k) { return i + ',' + j + ',' + k; };
    for (var a = 0; a < A.n; a++) {
      var k = key(Math.floor(A.x[a] / cell), Math.floor(A.y[a] / cell), Math.floor(A.z[a] / cell)), l = grid.get(k);
      if (l) l.push(a); else grid.set(k, [a]);
    }
    for (a = 0; a < A.n; a++) {
      var ci = Math.floor(A.x[a] / cell), cj = Math.floor(A.y[a] / cell), ck = Math.floor(A.z[a] / cell);
      for (var di = -1; di <= 1; di++) for (var dj = -1; dj <= 1; dj++) for (var dk = -1; dk <= 1; dk++) {
        var list = grid.get(key(ci + di, cj + dj, ck + dk));
        if (!list) continue;
        for (var q = 0; q < list.length; q++) {
          var o = list[q];
          if (o <= a) continue;
          var dx = A.x[a] - A.x[o], dy = A.y[a] - A.y[o], dz = A.z[a] - A.z[o], d2 = dx * dx + dy * dy + dz * dz;
          var cut = A.el[a] === 'S' || A.el[o] === 'S' ? 2.1 : 1.9;
          if (d2 < cut * cut && d2 > 0.16) out.push(a, o);
        }
      }
    }
    return out;
  }

  // AlphaFold's confidence colours: very high, confident, low, very low.
  function plddtColor(v) { return v >= 90 ? '#0053d6' : v >= 70 ? '#65cbf3' : v >= 50 ? '#ffdb13' : '#ff7d45'; }

  var cache = {};
  // The canonical model for a UniProt accession: { acc, entryId, model }.
  function fetchModel(acc) {
    if (cache[acc]) return cache[acc];
    cache[acc] = (async function () {
      var list = await G.net.fetchRetry(API + encodeURIComponent(acc)).then(function (r) {
        if (r.status === 404) return [];
        if (!r.ok) throw new Error('AlphaFold DB: HTTP ' + r.status);
        return r.json();
      });
      var e = (list || []).find(function (x) { return x.uniprotAccession === acc; }) || (list || [])[0];
      if (!e) throw new Error('No AlphaFold model for ' + acc + ' (AlphaFold DB has none for human proteins over 2,700 residues)');
      var text = await G.net.fetchRetry(e.cifUrl).then(function (r) { if (!r.ok) throw new Error('AlphaFold file: HTTP ' + r.status); return r.text(); });
      return { acc: acc, entryId: e.entryId, version: e.latestVersion, amUrl: e.amAnnotationsUrl || null, model: parseCif(text) };
    })();
    cache[acc].catch(function () { delete cache[acc]; });
    return cache[acc];
  }

  // ----- side chains from internal coordinates (for showing a missense change)
  //
  // Each atom is placed from three earlier atoms (a, b, c) by bond length c-d,
  // angle b-c-d and torsion a-b-c-d (the NeRF method). Torsions named chi1 to chi4
  // take the residue's rotamer; numbers are fixed; 'chiN+x' is chiN plus x.
  // Ideal geometry, close to Engh and Huber; rotamers are the most common ones
  // (Lovell et al. 2000 penultimate library). A modelled side chain, not a prediction.
  var CB = ['CB', 'C', 'C', 'N', 'CA', 1.53, 110.5, 'cb'];
  var SC = {
    ALA: [],
    SER: [['OG', 'O', 'N', 'CA', 'CB', 1.42, 111.1, 'chi1']],
    CYS: [['SG', 'S', 'N', 'CA', 'CB', 1.81, 114.0, 'chi1']],
    THR: [['OG1', 'O', 'N', 'CA', 'CB', 1.43, 109.2, 'chi1'], ['CG2', 'C', 'N', 'CA', 'CB', 1.52, 111.1, 'chi1+120']],
    VAL: [['CG1', 'C', 'N', 'CA', 'CB', 1.52, 110.7, 'chi1'], ['CG2', 'C', 'N', 'CA', 'CB', 1.52, 110.4, 'chi1-120']],
    ILE: [['CG1', 'C', 'N', 'CA', 'CB', 1.53, 110.4, 'chi1'], ['CG2', 'C', 'N', 'CA', 'CB', 1.52, 110.5, 'chi1+120'],
      ['CD1', 'C', 'CA', 'CB', 'CG1', 1.52, 113.8, 'chi2']],
    LEU: [['CG', 'C', 'N', 'CA', 'CB', 1.53, 116.1, 'chi1'], ['CD1', 'C', 'CA', 'CB', 'CG', 1.52, 110.3, 'chi2'],
      ['CD2', 'C', 'CA', 'CB', 'CG', 1.52, 110.6, 'chi2-120']],
    MET: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 114.0, 'chi1'], ['SD', 'S', 'CA', 'CB', 'CG', 1.81, 112.7, 'chi2'],
      ['CE', 'C', 'CB', 'CG', 'SD', 1.79, 100.6, 'chi3']],
    PHE: [['CG', 'C', 'N', 'CA', 'CB', 1.50, 113.9, 'chi1'], ['CD1', 'C', 'CA', 'CB', 'CG', 1.39, 120.7, 'chi2'],
      ['CD2', 'C', 'CA', 'CB', 'CG', 1.39, 120.7, 'chi2+180'], ['CE1', 'C', 'CB', 'CG', 'CD1', 1.39, 120.7, 180],
      ['CE2', 'C', 'CB', 'CG', 'CD2', 1.39, 120.7, 180], ['CZ', 'C', 'CG', 'CD1', 'CE1', 1.39, 120.0, 0]],
    TYR: [['CG', 'C', 'N', 'CA', 'CB', 1.51, 113.9, 'chi1'], ['CD1', 'C', 'CA', 'CB', 'CG', 1.39, 120.8, 'chi2'],
      ['CD2', 'C', 'CA', 'CB', 'CG', 1.39, 120.8, 'chi2+180'], ['CE1', 'C', 'CB', 'CG', 'CD1', 1.39, 121.2, 180],
      ['CE2', 'C', 'CB', 'CG', 'CD2', 1.39, 121.2, 180], ['CZ', 'C', 'CG', 'CD1', 'CE1', 1.38, 119.6, 0],
      ['OH', 'O', 'CD1', 'CE1', 'CZ', 1.38, 119.9, 180]],
    TRP: [['CG', 'C', 'N', 'CA', 'CB', 1.50, 114.1, 'chi1'], ['CD1', 'C', 'CA', 'CB', 'CG', 1.37, 127.1, 'chi2'],
      ['CD2', 'C', 'CA', 'CB', 'CG', 1.43, 126.6, 'chi2+180'], ['NE1', 'N', 'CB', 'CG', 'CD1', 1.38, 110.2, 180],
      ['CE2', 'C', 'CB', 'CG', 'CD2', 1.41, 107.2, 180], ['CE3', 'C', 'CB', 'CG', 'CD2', 1.40, 133.9, 0],
      ['CZ2', 'C', 'CG', 'CD2', 'CE2', 1.40, 122.4, 180], ['CZ3', 'C', 'CG', 'CD2', 'CE3', 1.39, 118.7, 180],
      ['CH2', 'C', 'CD2', 'CE2', 'CZ2', 1.37, 117.5, 0]],
    HIS: [['CG', 'C', 'N', 'CA', 'CB', 1.50, 113.7, 'chi1'], ['ND1', 'N', 'CA', 'CB', 'CG', 1.38, 122.7, 'chi2'],
      ['CD2', 'C', 'CA', 'CB', 'CG', 1.36, 131.0, 'chi2+180'], ['CE1', 'C', 'CB', 'CG', 'ND1', 1.32, 109.0, 180],
      ['NE2', 'N', 'CB', 'CG', 'CD2', 1.37, 107.0, 180]],
    ASP: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 113.0, 'chi1'], ['OD1', 'O', 'CA', 'CB', 'CG', 1.25, 119.2, 'chi2'],
      ['OD2', 'O', 'CA', 'CB', 'CG', 1.25, 118.2, 'chi2+180']],
    ASN: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 112.6, 'chi1'], ['OD1', 'O', 'CA', 'CB', 'CG', 1.23, 120.8, 'chi2'],
      ['ND2', 'N', 'CA', 'CB', 'CG', 1.33, 116.4, 'chi2+180']],
    GLU: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 114.1, 'chi1'], ['CD', 'C', 'CA', 'CB', 'CG', 1.52, 113.3, 'chi2'],
      ['OE1', 'O', 'CB', 'CG', 'CD', 1.25, 119.0, 'chi3'], ['OE2', 'O', 'CB', 'CG', 'CD', 1.25, 118.1, 'chi3+180']],
    GLN: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 114.1, 'chi1'], ['CD', 'C', 'CA', 'CB', 'CG', 1.52, 112.8, 'chi2'],
      ['OE1', 'O', 'CB', 'CG', 'CD', 1.23, 120.9, 'chi3'], ['NE2', 'N', 'CB', 'CG', 'CD', 1.33, 116.5, 'chi3+180']],
    LYS: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 114.1, 'chi1'], ['CD', 'C', 'CA', 'CB', 'CG', 1.52, 111.3, 'chi2'],
      ['CE', 'C', 'CB', 'CG', 'CD', 1.52, 111.3, 'chi3'], ['NZ', 'N', 'CG', 'CD', 'CE', 1.49, 111.9, 'chi4']],
    ARG: [['CG', 'C', 'N', 'CA', 'CB', 1.52, 114.1, 'chi1'], ['CD', 'C', 'CA', 'CB', 'CG', 1.52, 111.3, 'chi2'],
      ['NE', 'N', 'CB', 'CG', 'CD', 1.46, 112.0, 'chi3'], ['CZ', 'C', 'CG', 'CD', 'NE', 1.33, 124.2, 'chi4'],
      ['NH1', 'N', 'CD', 'NE', 'CZ', 1.33, 120.0, 0], ['NH2', 'N', 'CD', 'NE', 'CZ', 1.33, 119.6, 180]],
    PRO: [['CG', 'C', 'N', 'CA', 'CB', 1.50, 104.5, 'chi1'], ['CD', 'C', 'CA', 'CB', 'CG', 1.51, 105.5, 'chi2']]
  };
  // Most common rotamer per residue (degrees).
  var ROTAMER = { SER: [62], CYS: [-65], THR: [59], VAL: [175], ILE: [-65, 170], LEU: [-65, 175], MET: [-65, 180, -70],
    PHE: [-65, -85], TYR: [-65, -85], TRP: [-65, 95], HIS: [-65, -70], ASP: [-70, -15], ASN: [-65, -20],
    GLU: [-65, 180, -10], GLN: [-65, 180, -25], LYS: [-65, 180, 180, 180], ARG: [-65, 180, 180, 85], PRO: [30, -35] };
  var CB_TORSION = 122.6; // C-N-CA-CB for an L amino acid

  function v3(a) { return [a[0], a[1], a[2]]; }
  function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
  function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function norm(a) { var l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  // NeRF: d from a, b, c with |cd| = bond, angle bcd, torsion abcd (degrees).
  function place(a, b, c, bond, angle, torsion) {
    var R = Math.PI / 180, th = angle * R, ph = -torsion * R; // sign so that torsion(a, b, c, place(...)) === torsion
    var d2 = [-bond * Math.cos(th), bond * Math.sin(th) * Math.cos(ph), bond * Math.sin(th) * Math.sin(ph)];
    var bc = norm(sub(c, b)), n = norm(cross(sub(b, a), bc)), m = cross(n, bc);
    return [c[0] + bc[0] * d2[0] + m[0] * d2[1] + n[0] * d2[2], c[1] + bc[1] * d2[0] + m[1] * d2[1] + n[1] * d2[2], c[2] + bc[2] * d2[0] + m[2] * d2[1] + n[2] * d2[2]];
  }
  function torsion(a, b, c, d) {
    var b1 = sub(b, a), b2 = sub(c, b), b3 = sub(d, c), n1 = cross(b1, b2), n2 = cross(b2, b3);
    var m1 = cross(n1, norm(b2));
    return Math.atan2(dot(m1, n2), dot(n1, n2)) * 180 / Math.PI;
  }

  // Side chain atoms of residue type aa on backbone {N, CA, C} ([x,y,z] each), for
  // chi angles chis (default: the most common rotamer). Returns [{name, el, p}], CB first.
  function buildSideChain(aa, bb, chis) {
    if (aa === 'GLY' || !SC[aa]) return [];
    chis = chis || ROTAMER[aa] || [];
    var pos = { N: bb.N, CA: bb.CA, C: bb.C }, out = [];
    pos.CB = place(bb.C, bb.N, bb.CA, CB[5], CB[6], CB_TORSION);
    out.push({ name: 'CB', el: 'C', p: pos.CB });
    SC[aa].forEach(function (d) {
      var t = d[7];
      if (typeof t === 'string') { var m = /^chi(\d)([+-]\d+)?$/.exec(t); t = (chis[+m[1] - 1] || 0) + (m[2] ? +m[2] : 0); }
      pos[d[0]] = place(pos[d[2]], pos[d[3]], pos[d[4]], d[5], d[6], t);
      out.push({ name: d[0], el: d[1], p: pos[d[0]] });
    });
    return out;
  }

  // The chi angles of a residue as found in a model (for checks, and to keep the
  // normal residue's rotamer where the new one shares it).
  function measureChis(aa, at) {
    var defs = SC[aa] || [], chis = [];
    defs.forEach(function (d) {
      var m = typeof d[7] === 'string' && /^chi(\d)$/.exec(d[7]);
      if (!m || !at[d[0]] || !at[d[2]] || !at[d[3]] || !at[d[4]]) return;
      chis[+m[1] - 1] = torsion(at[d[2]], at[d[3]], at[d[4]], at[d[0]]);
    });
    return chis;
  }

  // Residue properties for describing a change (facts, not a judgement).
  var AA = {
    ALA: ['A', 'hydrophobic', 88.6, 1.8], ARG: ['R', 'positive', 173.4, -4.5], ASN: ['N', 'polar', 114.1, -3.5], ASP: ['D', 'negative', 111.1, -3.5],
    CYS: ['C', 'polar', 108.5, 2.5], GLN: ['Q', 'polar', 143.8, -3.5], GLU: ['E', 'negative', 138.4, -3.5], GLY: ['G', 'special (flexible)', 60.1, -0.4],
    HIS: ['H', 'positive (weak)', 153.2, -3.2], ILE: ['I', 'hydrophobic', 166.7, 4.5], LEU: ['L', 'hydrophobic', 166.7, 3.8], LYS: ['K', 'positive', 168.6, -3.9],
    MET: ['M', 'hydrophobic', 162.9, 1.9], PHE: ['F', 'hydrophobic, aromatic', 189.9, 2.8], PRO: ['P', 'special (rigid)', 112.7, -1.6], SER: ['S', 'polar', 89.0, -0.8],
    THR: ['T', 'polar', 116.1, -0.7], TRP: ['W', 'hydrophobic, aromatic', 227.8, -0.9], TYR: ['Y', 'polar, aromatic', 193.6, -1.3], VAL: ['V', 'hydrophobic', 140.0, 4.2]
  };
  var ONE = {}; Object.keys(AA).forEach(function (k) { ONE[AA[k][0]] = k; });

  // The change at one residue: properties, burial, clashes of the modelled side chain.
  // model: parseCif result; seq: residue number; alt: three-letter code.
  function mutate(model, seq, alt) {
    var res = model.bySeq[seq], A = model.atoms;
    if (!res) return null;
    var at = {};
    res.atoms.forEach(function (i) { at[A.name[i]] = [A.x[i], A.y[i], A.z[i]]; });
    var bb = { N: at.N, CA: at.CA, C: at.C };
    var refChis = measureChis(res.name, at), altChis = (ROTAMER[alt] || []).slice();
    if (refChis.length && altChis.length && SC[res.name] && SC[alt] && SC[res.name][0][0] === SC[alt][0][0]) altChis[0] = refChis[0]; // keep chi1 when the gamma atom is shared
    var side = buildSideChain(alt, bb, altChis);
    // clashes: modelled atoms (beyond CB) within 3.0 A of atoms of other residues, minus the neighbours' backbone bonded to it
    var clashes = [], near = 0, cb = side[0] ? side[0].p : bb.CA;
    for (var i = 0; i < A.n; i++) {
      if (A.seq[i] === seq) continue;
      var dx = A.x[i] - cb[0], dy = A.y[i] - cb[1], dz = A.z[i] - cb[2];
      if (dx * dx + dy * dy + dz * dz < 100) near++;
      side.forEach(function (s, k) {
        if (k === 0) return;
        var ex = A.x[i] - s.p[0], ey = A.y[i] - s.p[1], ez = A.z[i] - s.p[2], d = Math.sqrt(ex * ex + ey * ey + ez * ez);
        if (d < 3.0) clashes.push({ atom: s.name, with: i, d: d });
      });
    }
    var r = AA[res.name] || [], a = AA[alt] || [];
    return {
      seq: seq, ref: res.name, alt: alt, side: side, clashes: clashes, near: near,
      buried: near >= 120 ? 'buried' : near >= 70 ? 'partly buried' : 'on the surface',
      plddt: res.plddt, ss: res.ss,
      refClass: r[1], altClass: a[1], dVolume: a[2] - r[2], dHydropathy: a[3] - r[3]
    };
  }

  // AlphaMissense for one protein (CSV from AlphaFold DB: protein_variant, am_pathogenicity,
  // am_class): { byChange: {G551D: {score, cls}}, mean: {seq: mean over the substitutions} }.
  function parseAlphaMissense(text) {
    var byChange = {}, sum = {}, cnt = {};
    text.split(/\r?\n/).forEach(function (l, i) {
      if (!i || !l) return;
      var f = l.split(','), m = /^([A-Z])(\d+)([A-Z])$/.exec(f[0]);
      if (!m) return;
      var sc = +f[1], seq = +m[2];
      byChange[f[0]] = { score: sc, cls: f[2] };
      sum[seq] = (sum[seq] || 0) + sc; cnt[seq] = (cnt[seq] || 0) + 1;
    });
    var mean = {};
    Object.keys(sum).forEach(function (k) { mean[k] = sum[k] / cnt[k]; });
    return { byChange: byChange, mean: mean };
  }
  var amCache = {};
  function fetchAlphaMissense(url) {
    if (!amCache[url]) {
      amCache[url] = fetch(url).then(function (r) { if (!r.ok) throw new Error('AlphaMissense: HTTP ' + r.status); return r.text(); }).then(parseAlphaMissense);
      amCache[url].catch(function () { delete amCache[url]; });
    }
    return amCache[url];
  }

  G.structure = { parseAlphaMissense: parseAlphaMissense, fetchAlphaMissense: fetchAlphaMissense, parseCif: parseCif, bonds: bonds, plddtColor: plddtColor, fetchModel: fetchModel,
    buildSideChain: buildSideChain, measureChis: measureChis, torsion: torsion, mutate: mutate, AA: AA, ONE: ONE };
})(globalThis.G = globalThis.G || {});
