#!/usr/bin/env python3
"""CMFGEN populations as a SuperLite NLTE state (output.nltex_state format), for a
warm start with in_nltex_hold: SuperLite's transport and formal integral on CMFGEN's
own populations, which separates the radiative transfer from the rate equations.

The SuperLite model atoms must be those converted from the CMFGEN model's atomic
data (Tools/cmfgen2superlite_atom.py with the model's F_TO_S files): the levels of
each ion in the order of CMFGEN's POP files, summed into the super-levels of the
SuperLite atom; the ions below the lowest ion of an element in SuperLite go to its
ground state, those above its highest ion to the top state.  Each cell gets the
volume average of the state densities between its edges (log-linear in r between
CMFGEN's depths) over that of the element, CMFGEN's n_e alike, and the T_e of the
SuperLite structure.

usage: cmfgen2superlite_state.py <CMFGEN model dir> <SuperLite run dir> [out]
  (the run dir has input.str and input.par with in_nltex_dir, in_nltex_atoms;
  out: output.nltex_state there by default)"""
import os
import re
import sys
import numpy as np

POPFILE = {1: 'POPHYD', 2: 'POPHE', 6: 'POPCARB', 7: 'POPNIT', 8: 'POPOXY', 10: 'POPNEON',
           11: 'POPSOD', 12: 'POPMG', 13: 'POPALUM', 14: 'POPSIL', 15: 'POPPHOS', 16: 'POPSUL',
           17: 'POPCHL', 18: 'POPARG', 20: 'POPCAL', 24: 'POPCHRO', 25: 'POPMAN', 26: 'POPIRON',
           27: 'POPCOB', 28: 'POPNICK'}
# (other names CMFGEN uses for the same files)
POPALT = {12: 'POPMAG'}


def popfile(mod, z):
    for name in (POPFILE[z], POPALT.get(z)):
        if name and os.path.exists(os.path.join(mod, name)):
            return os.path.join(mod, name)
    return os.path.join(mod, POPFILE[z])


def numbers(s):
    return [float(x.replace('D', 'E')) for x in s.split()]


def read_pops(path, nd):
    """element density (nd) and per ion (name, levels (nd, nlev)) in file order, and the top ion (nd);
    depth 1 is the outer boundary"""
    lines = open(path, errors='replace').read().splitlines()
    i = 0
    while i < len(lines) and not re.match(r'^\s*[-+]?\d+\.\d+E', lines[i]):
        i += 1
    vals = []
    while len(vals) < nd:
        vals += numbers(lines[i])
        i += 1
    total = np.array(vals[:nd])
    ions, top = [], None
    while i < len(lines):
        m = re.match(r'\s*Number of (\S+) levels:\s+(\d+)', lines[i])
        if not m:
            i += 1
            continue
        name, nlev = m.group(1), int(m.group(2))
        i += 1
        while i < len(lines) and not re.match(r'^\s*[-+]?\d+\.\d+E', lines[i]):
            i += 1
        vals = []
        while i < len(lines) and len(vals) < (nlev + 1)*nd:
            if re.match(r'^\s*[A-Za-z]', lines[i]):
                break
            vals += numbers(lines[i])
            i += 1
        ions.append((name, np.array(vals[:nlev*nd]).reshape(nd, nlev)))
        if len(vals) >= (nlev + 1)*nd:
            top = np.array(vals[nlev*nd:(nlev + 1)*nd])
    return total, ions, top


def read_atom(path):
    """Z, ion stage (1: neutral), superlevel index of each level"""
    z = ion = nlev = None
    sl = []
    lines = open(path).read().splitlines()
    for k, s in enumerate(lines):
        p = s.split()
        if not p or p[0].startswith('#'):
            continue
        if p[0] == 'element':
            z = int(p[1])
        elif p[0] == 'ion':
            ion = int(p[1])
        elif p[0] == 'nlev':
            nlev = int(p[1])
            j = k + 1
            while len(sl) < nlev:
                q = lines[j].split()
                j += 1
                if q and not q[0].startswith('#'):
                    sl.append(int(q[-1]))
            break
    return z, ion, np.array(sl)


def namelist(path, key):
    m = re.search(r'\b' + key + r"\s*=\s*'([^']*)'", open(path).read(), re.I)
    return m.group(1) if m else None


def cell_average(rd, y, r1, r2):
    """volume average over [r1, r2] of y(r) log-linear in r between the depths rd (increasing)"""
    rr = np.linspace(r1, r2, 41)
    ly = np.interp(np.log(rr), np.log(rd), np.log(np.maximum(y, 1e-300)))
    w = rr**2
    return np.trapezoid(np.exp(ly)*w, rr)/np.trapezoid(w, rr)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    mod, run = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else os.path.join(run, 'output.nltex_state')
    par = os.path.join(run, 'input.par')
    adir = namelist(par, 'in_nltex_dir')
    atoms = namelist(par, 'in_nltex_atoms').split()
    if not os.path.isabs(adir):
        adir = os.path.join(run, adir)
# SuperLite elements: ions in order, the super-levels of each, and the top state
    els = {}
    for a in atoms:
        z, ion, sl = read_atom(os.path.join(adir, a))
        els.setdefault(z, []).append((ion, sl))
    zs = sorted(els)
    for z in zs:
        els[z].sort(key=lambda t: t[0])
# CMFGEN depths, T, n_e (RVTJ)
    rv = open(os.path.join(mod, 'RVTJ')).read().splitlines()
    nd = int(rv[3].split(':')[1])

    def section(title):
        k = [i for i, s in enumerate(rv) if s.strip() == title][0] + 1
        v = []
        while len(v) < nd:
            v += numbers(rv[k])
            k += 1
        return np.array(v)
    r = section('Radius (10^10 cm)')*1e10
    ne = section('Electron density')
    o = np.argsort(r)
# SuperLite cells
    rows = [l.split() for l in open(os.path.join(run, 'input.str')) if not l.startswith('#')]
    hdr = [l for l in open(os.path.join(run, 'input.str')) if l.startswith('#')]
    cols = hdr[-1].lstrip('#').split()
    ix = {c: k for k, c in enumerate(cols)}
    edges = np.array([[float(q[ix['x_left']]), float(q[ix['x_right']])] for q in rows])
    mass = np.array([float(q[ix['mass']]) for q in rows])
    temp = np.array([float(q[ix['temp']]) for q in rows])
    nc = len(rows)
    vol = 4*np.pi/3*(edges[:, 1]**3 - edges[:, 0]**3)
    mext = np.zeros(nc)
    mext[-1] = 0.5*mass[-1]
    for c in range(nc - 2, -1, -1):
        mext[c] = mext[c + 1] + 0.5*(mass[c + 1] + mass[c])
# states per element: densities per CMFGEN depth
    xs, nss = [], []
    for z in zs:
        total, ions, top = read_pops(popfile(mod, z), nd)
        sl_ions = els[z]
        ion0 = sl_ions[0][0]
        nst = sum(int(sl.max()) for _, sl in sl_ions) + 1
        st = np.zeros((nd, nst))
        # CMFGEN's ions in its file order are stages 1, 2, ... from its lowest ion; match by level count
        # from the lowest SuperLite ion upward
        cm_nlev = [p.shape[1] for _, p in ions]
        sl_nlev = [len(sl) for _, sl in sl_ions]
        k0 = None
        for k in range(len(ions) - len(sl_ions) + 1):
            if cm_nlev[k:k + len(sl_ions)] == sl_nlev:
                k0 = k
                break
        if k0 is None:
            sys.exit('Z=%d: CMFGEN ions %s do not match the SuperLite atoms %s' % (z, cm_nlev, sl_nlev))
        lumped = np.zeros(nd)
        for k in range(k0):
            lumped += ions[k][1].sum(axis=1)
        off = 0
        for j, (ion, sl) in enumerate(sl_ions):
            p = ions[k0 + j][1]
            for s in range(1, int(sl.max()) + 1):
                st[:, off + s - 1] = p[:, sl == s].sum(axis=1)
            if j == 0:
                st[:, 0] += lumped
            off += int(sl.max())
        above = np.zeros(nd)
        for k in range(k0 + len(sl_ions), len(ions)):
            above += ions[k][1].sum(axis=1)
        st[:, -1] = (top if top is not None else 0.0) + above
        dev = np.max(np.abs(st.sum(axis=1)/total - 1))
        print('Z=%2d: %d SuperLite states from CMFGEN ions %s (lowest %d lumped: max %.1e of the element); '
              'sum/total - 1 at most %.1e' % (z, nst, [n for n, _ in ions], k0, np.max(lumped/total), dev))
        xs.append(st)
        nss.append(nst)
# cell averages
    nst_all = sum(nss)
    X = np.zeros((nc, nst_all))
    NE = np.zeros(nc)
    for c in range(nc):
        r1, r2 = edges[c]
        NE[c] = cell_average(r[o], ne[o], r1, r2)
        off = 0
        for st in xs:
            dens = np.array([cell_average(r[o], st[o, s], r1, r2) for s in range(st.shape[1])])
            X[c, off:off + st.shape[1]] = dens/dens.sum()
            off += st.shape[1]
    with open(out, 'w') as f:
        f.write('# SuperLite NLTE state populations per cell (in_nltex_prev of a later model): '
                'CMFGEN model %s\n' % os.path.abspath(mod))
        f.write('nel nst ncell%9d%9d%9d\n' % (len(zs), nst_all, nc))
        f.write('elements Z ns' + ''.join(' %3d %5d' % (z, n) for z, n in zip(zs, nss)) + '\n')
        f.write('# per cell: exterior mass of the cell centre [g], mass [g], density [g/cm^3], '
                'T_e [K], n_e [cm^-3]; then the fractions of the states of each element\n')
        for c in range(nc):
            f.write(''.join('%16.7E' % v for v in (mext[c], mass[c], mass[c]/vol[c], temp[c], NE[c])) + '\n')
            for k in range(0, nst_all, 8):
                f.write(''.join('%16.7E' % v for v in X[c, k:k + 8]) + '\n')
    print('wrote', out, '(%d cells, %d states)' % (nc, nst_all))


if __name__ == '__main__':
    main()
