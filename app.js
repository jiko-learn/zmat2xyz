// Utility: degrees to radians
function degToRad(deg) {
  return (deg * Math.PI) / 180;
}

// Vector math helpers
function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function scale(v, s) { return [v[0] * s, v[1] * s, v[2] * s]; }
function dot(a, b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
function cross(a, b) { return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]; }
function norm(v) { return Math.sqrt(dot(v, v)); }
function normalize(v) {
  const n = norm(v);
  if (n === 0) return [0, 0, 0];
  return [v[0]/n, v[1]/n, v[2]/n];
}

function nearlyColinear(a, b) {
  const n = norm(cross(a, b));
  return n < 1e-8;
}

// Parse Z-matrix (simple, whitespace-separated tokens). Supports common format:
// 1) Elem
// 2) Elem j r
// 3) Elem j r k angle
// 4+) Elem j r k angle l dihedral
// Angles are degrees; indices are 1-based.
function parseZMatrix(text) {
  const lines = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length > 0 && !l.startsWith('#'));
  if (lines.length === 0) throw new Error('Empty input.');

  const atoms = [];
  for (let i = 0; i < lines.length; i++) {
    const tokens = lines[i].split(/\s+/);
    if (tokens.length < 1) throw new Error(`Invalid line ${i+1}`);

    const elem = tokens[0];
    const rec = { element: elem };

    if (i === 0) {
      // A1
      // nothing more
    } else if (i === 1) {
      // Elem j r
      if (tokens.length < 3) throw new Error(`Line ${i+1}: expected 3 tokens (Elem j r)`);
      rec.j = parseInt(tokens[1], 10) - 1;
      rec.r = parseFloat(tokens[2]);
    } else if (i === 2) {
      // Elem j r k angle
      if (tokens.length < 5) throw new Error(`Line ${i+1}: expected 5 tokens (Elem j r k angle)`);
      rec.j = parseInt(tokens[1], 10) - 1;
      rec.r = parseFloat(tokens[2]);
      rec.k = parseInt(tokens[3], 10) - 1;
      rec.theta = degToRad(parseFloat(tokens[4]));
    } else {
      // Elem j r k angle l dihedral
      if (tokens.length < 7) throw new Error(`Line ${i+1}: expected 7 tokens (Elem j r k angle l dihedral)`);
      rec.j = parseInt(tokens[1], 10) - 1;
      rec.r = parseFloat(tokens[2]);
      rec.k = parseInt(tokens[3], 10) - 1;
      rec.theta = degToRad(parseFloat(tokens[4]));
      rec.l = parseInt(tokens[5], 10) - 1;
      rec.phi = degToRad(parseFloat(tokens[6]));
    }
    atoms.push(rec);
  }
  return atoms;
}

// Convert parsed Z-matrix to Cartesian XYZ coordinates.
// Returns array of { element, x, y, z }
function zmatToXYZ(parsed) {
  const n = parsed.length;
  const coords = new Array(n);

  // 1) First atom at origin
  coords[0] = [0, 0, 0];

  if (n >= 2) {
    const a2 = parsed[1];
    if (a2.j == null || a2.r == null) throw new Error('Second atom requires j and r');
    coords[1] = [a2.r, 0, 0]; // place along +x
  }

  if (n >= 3) {
    const a3 = parsed[2];
    if (a3.j == null || a3.k == null || a3.r == null || a3.theta == null) throw new Error('Third atom requires j, r, k, theta');

    const j = a3.j; // bond reference
    const k = a3.k; // angle reference
    if (j < 0 || j >= 2 || k < 0 || k >= 2) throw new Error('Third atom references must be among first two atoms');

    const pj = coords[j];
    const pk = coords[k];

    // u: direction j -> k
    const u = normalize(sub(pk, pj));

    // Choose perpendicular e2 in a stable way in the plane z=0
    let arbitrary = [0, 0, 1];
    if (nearlyColinear(u, arbitrary)) arbitrary = [0, 1, 0];
    const e2 = normalize(cross(u, arbitrary));

    const r = a3.r;
    const theta = a3.theta;
    // position: pj + r * ( cos(theta) * u + sin(theta) * e2 )
    const dir = add(scale(u, Math.cos(theta)), scale(e2, Math.sin(theta)));
    coords[2] = add(pj, scale(dir, r));
  }

  for (let i = 3; i < n; i++) {
    const rec = parsed[i];
    const j = rec.j; // bond reference
    const k = rec.k; // angle reference
    const l = rec.l; // dihedral reference
    const r = rec.r;
    const theta = rec.theta; // angle at j-k-new
    const phi = rec.phi;     // dihedral with l-k-j-new

    if ([j, k, l].some(idx => idx == null)) throw new Error(`Atom ${i+1}: missing references j,k,l`);
    if (j < 0 || j >= i || k < 0 || k >= i || l < 0 || l >= i) throw new Error(`Atom ${i+1}: references must point to already placed atoms`);

    const pj = coords[j];
    const pk = coords[k];
    const pl = coords[l];

    // Build local orthonormal frame at j:
    // u = normalize(Rk - Rj)
    // a = normalize(Rl - Rj)
    // v = normalize(cross(u, a))
    // w = cross(v, u)
    const u = normalize(sub(pk, pj));
    let a = sub(pl, pj);
    a = normalize(a);
    let v = cross(u, a);
    v = normalize(v);

    // Handle degeneracy (l aligned with k relative to j)
    if (norm(v) === 0 || Number.isNaN(v[0])) {
      let arbitrary = [0, 0, 1];
      if (nearlyColinear(u, arbitrary)) arbitrary = [0, 1, 0];
      v = normalize(cross(u, arbitrary));
    }
    const w = cross(v, u); // already orthonormal if u,v are orthonormal

    // position: pj + r * ( cos(theta) * u + sin(theta) * ( cos(phi) * w + sin(phi) * v ) )
    const dir = add(
      scale(u, Math.cos(theta)),
      scale(add(scale(w, Math.cos(phi)), scale(v, Math.sin(phi))), Math.sin(theta))
    );
    coords[i] = add(pj, scale(dir, r));
  }

  // Package with elements
  return parsed.map((rec, i) => ({
    element: rec.element,
    x: coords[i][0],
    y: coords[i][1],
    z: coords[i][2]
  }));
}

function formatXYZ(atoms, comment = 'Converted from Z-Matrix') {
  const lines = [];
  lines.push(String(atoms.length));
  lines.push(comment);
  for (const a of atoms) {
    lines.push(`${a.element} ${a.x.toFixed(6)} ${a.y.toFixed(6)} ${a.z.toFixed(6)}`);
  }
  return lines.join('\n');
}

function centerAtoms(atoms) {
  const cx = atoms.reduce((s,a)=>s+a.x,0)/atoms.length;
  const cy = atoms.reduce((s,a)=>s+a.y,0)/atoms.length;
  const cz = atoms.reduce((s,a)=>s+a.z,0)/atoms.length;
  return atoms.map(a => ({...a, x: a.x - cx, y: a.y - cy, z: a.z - cz }));
}

// UI wiring
const els = {
  zmatInput: document.getElementById('zmatInput'),
  convertBtn: document.getElementById('convertBtn'),
  view3DBtn: document.getElementById('view3DBtn'),
  message: document.getElementById('message'),
  xyzOutput: document.getElementById('xyzOutput'),
  copyXYZBtn: document.getElementById('copyXYZBtn'),
  downloadXYZBtn: document.getElementById('downloadXYZBtn'),
  loadExampleBtn: document.getElementById('loadExampleBtn'),
  clearInputBtn: document.getElementById('clearInputBtn'),
  viewerDiv: document.getElementById('viewer'),
  resetViewBtn: document.getElementById('resetViewBtn'),
  centerCheckbox: document.getElementById('centerCheckbox'),
  showLabelsCheckbox: document.getElementById('showLabelsCheckbox'),
  showIndexLabelsCheckbox: document.getElementById('showIndexLabelsCheckbox')
};

let viewer = null; // 3Dmol viewer (fallback)
let jsmolApplet = null; // JSmol applet instance

function setMessage(text, type = '') {
  els.message.textContent = text;
  els.message.className = type ? type : '';
}

function convertNow() {
  try {
    setMessage('Converting...', '');
    const parsed = parseZMatrix(els.zmatInput.value);
    let atoms = zmatToXYZ(parsed);
    if (els.centerCheckbox.checked) atoms = centerAtoms(atoms);
    els.xyzOutput.value = formatXYZ(atoms);
    setMessage('Conversion successful.', 'success');
  } catch (err) {
    console.error(err);
    setMessage(err.message || 'Conversion failed.', 'error');
  }
}

// ensureViewer removed in favor of generic renderXYZInElement

function ensureJSmolApplet() {
  if (!window.Jmol) return null;
  if (jsmolApplet) return jsmolApplet;
  const Info = {
    width: els.viewerDiv.clientWidth || 600,
    height: els.viewerDiv.clientHeight || 400,
    use: 'HTML5',
    debug: false,
    color: '0x000000',
    // Use hosted assets to avoid bundling j2s/php locally
    j2sPath: 'https://chemapps.stolaf.edu/jmol/jsmol/j2s',
    serverURL: 'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php',
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true
  };
  jsmolApplet = Jmol.getApplet('jsmolApplet', Info);
  // Inject the applet HTML into the viewer container
  els.viewerDiv.innerHTML = Jmol.getAppletHtml(jsmolApplet);
  return jsmolApplet;
}

function renderWithJSmol(xyz) {
  const applet = ensureJSmolApplet();
  if (!applet) return false;
  // Ensure XYZ string ends with newline and is valid
  const data = xyz.endsWith('\n') ? xyz : xyz + '\n';
  
  let script = 'load DATA "model"\n' + data + 'END "model" {format xyz};\n' +
               'select *; spacefill 25%; wireframe 0.15; set zoomLarge false; zoomTo 0.0;';
  
  // Add labels if checkboxes are checked
  if (els.showLabelsCheckbox.checked && els.showIndexLabelsCheckbox.checked) {
    // Show both index and symbol (index first, no space)
    script += 'select *; label "%i%e"; set labelOffset 0 0; set labelSize 16; set labelColor black; set labelBgColor white;';
  } else if (els.showLabelsCheckbox.checked) {
    script += 'select *; label %e; set labelOffset 0 0; set labelSize 16; set labelColor black; set labelBgColor white;';
  } else if (els.showIndexLabelsCheckbox.checked) {
    script += 'select *; label %i; set labelOffset 0 0; set labelSize 16; set labelColor black; set labelBgColor white;';
  } else {
    script += 'select *; label OFF;';
  }
  
  Jmol.script(applet, script);
  return true;
}

function renderWith3DMol(xyz) {
  if (!window.$3Dmol) return false;
  if (!viewer) {
    viewer = $3Dmol.createViewer(els.viewerDiv, { backgroundColor: 'rgba(0,0,0,0)' });
  }
  viewer.clear();
  viewer.addModel(xyz, 'xyz');
  viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { scale: 0.28 } });
  
  // Add labels if checkboxes are checked
  if (els.showLabelsCheckbox.checked || els.showIndexLabelsCheckbox.checked) {
    // Parse the XYZ to get atom positions and elements for labeling
    const lines = xyz.trim().split('\n');
    const numAtoms = parseInt(lines[0]);
    
    // Create labels for each atom
    for (let i = 0; i < numAtoms && i + 2 < lines.length; i++) {
      const atomData = lines[i + 2].trim().split(/\s+/);
      if (atomData.length >= 4) {
        const element = atomData[0];
        const x = parseFloat(atomData[1]);
        const y = parseFloat(atomData[2]);
        const z = parseFloat(atomData[3]);
        
        // Determine label text based on checkboxes
        let labelText = '';
        if (els.showLabelsCheckbox.checked && els.showIndexLabelsCheckbox.checked) {
          labelText = `${i + 1}${element}`; // Show both index and symbol (index first, no space)
        } else if (els.showLabelsCheckbox.checked) {
          labelText = element; // Show only symbol
        } else if (els.showIndexLabelsCheckbox.checked) {
          labelText = String(i + 1); // Show only index
        }
        
        // Add label positioned slightly above the atom
        viewer.addLabel(labelText, {
          position: { x: x, y: y + 0.3, z: z }, // Offset y to place label above atom
          backgroundColor: 'rgba(255,255,255,0.9)',
          fontColor: 'black',
          backgroundOpacity: 0.9,
          fontSize: 14,
          fontFamily: 'Arial, sans-serif',
          alignment: 'center',
          showBackground: true,
          inFront: true,
          backgroundPadding: 1
        });
      }
    }
  }
  
  viewer.resize();
  viewer.zoomTo();
  viewer.render();
  return true;
}

function view3D() {
  if (!els.xyzOutput.value.trim()) {
    // If no XYZ yet, try converting automatically
    convertNow();
    if (!els.xyzOutput.value.trim()) return;
  }
  try {
    const xyz = els.xyzOutput.value;
    // Prefer JSmol when available
    const ok = renderWithJSmol(xyz) || renderWith3DMol(xyz);
    if (!ok) throw new Error('No supported 3D viewer found (JSmol or 3Dmol).');
    setMessage('3D view ready.', 'success');
  } catch (e) {
    console.error(e);
    setMessage(e.message || '3D viewer failed to initialize.', 'error');
  }
}

function copyXYZ() {
  if (!els.xyzOutput.value) return;
  navigator.clipboard.writeText(els.xyzOutput.value).then(() => {
    setMessage('XYZ copied to clipboard.', 'success');
  }).catch(() => setMessage('Copy failed.', 'error'));
}

function downloadXYZ() {
  if (!els.xyzOutput.value) return;
  const blob = new Blob([els.xyzOutput.value], { type: 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'molecule.xyz';
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function loadExample() {
  // Ammonia example: trigonal pyramid around nitrogen
  const example = [
    'N',
    'H 1 1.0100',
    'H 1 1.0100 2 107.80',
    'H 1 1.0100 2 107.80 3 120.00'
  ].join('\n');
  els.zmatInput.value = example;
  setMessage('Loaded example (NH3).', '');
}

function clearInput() {
  els.zmatInput.value = '';
  setMessage('');
}

function resetView() {
  if (jsmolApplet) {
    Jmol.script(jsmolApplet, 'zoomTo 0.0;');
    return;
  }
  if (viewer) {
    viewer.zoomTo();
    viewer.render();
  }
}

// Event listeners
els.convertBtn.addEventListener('click', convertNow);
els.view3DBtn.addEventListener('click', view3D);
els.copyXYZBtn.addEventListener('click', copyXYZ);
els.downloadXYZBtn.addEventListener('click', downloadXYZ);
els.loadExampleBtn.addEventListener('click', loadExample);
els.clearInputBtn.addEventListener('click', clearInput);
els.resetViewBtn.addEventListener('click', resetView);

// Re-render view when label checkboxes are toggled (if view is already displayed)
els.showLabelsCheckbox.addEventListener('change', () => {
  if (viewer || jsmolApplet) {
    view3D();
  }
});
els.showIndexLabelsCheckbox.addEventListener('change', () => {
  if (viewer || jsmolApplet) {
    view3D();
  }
});

// Auto-init example for convenience
loadExample();



