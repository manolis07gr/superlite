#!/usr/bin/env python3
"""SuperLite input.str of model A of Dessart et al. (2016, MNRAS 458, 2094) at one
epoch, from the vector curves of their Fig. 3 (T, V, rho and the radiative
luminosity of the 8 energy groups against radius, HERACLES), saved by the
extraction as modelA_<epoch>_raw.npz.  CMFGEN post-processed these snapshots
with T held fixed; here T is the temp column (in_tesolve = f) and the
luminosity L(r), the sum of the groups, the lum column (lumg source).
Composition: uniform, H 0.633 and He 0.36564 as their HERACLES runs, with
solar metals (Asplund et al. 2009) scaled into the remaining mass fraction.
usage: build_modelA.py raw.npz out_dir [--ncell N]"""
import sys, os, numpy as np
raw, out = sys.argv[1], sys.argv[2]
ncell = int(sys.argv[sys.argv.index('--ncell') + 1]) if '--ncell' in sys.argv else 120
d = np.load(raw)
def f(x, y):
    o = np.argsort(x); return x[o], y[o]
rT, T = f(d['rT'], d['T'] * 1e3)
rV, V = f(d['rV'], d['V'] * 1e8)               # cm/s
rR, lr = f(d['rR'], d['lrho'])
gr = [f(d['rL%d' % i], d['lL%d' % i]) for i in range(8)]
r0, r1 = max(rT.min(), rV.min(), rR.min(), 0.05), min(rT.max(), rV.max(), rR.max())
# cell edges: uniform in r, plus refinement where rho or V change fast
base = np.linspace(r0, r1, ncell // 2 + 1)
dl = np.abs(np.gradient(np.interp(base, rR, lr)))
dv = np.abs(np.gradient(np.interp(base, rV, V))) / V.max()
w = 1 + 40 * dl / dl.max() + 40 * dv / dv.max()
cw = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]))])
edges = np.interp(np.linspace(0, cw[-1], ncell + 1), cw, base) * 1e15
rc = 0.5 * (edges[1:] + edges[:-1])
x15 = lambda a: a / 1e15
rho = 10 ** np.interp(x15(rc), rR, lr)
vol = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
mass = rho * vol
tc = np.interp(x15(rc), rT, T)
vl, vr = np.interp(x15(edges[:-1]), rV, V), np.interp(x15(edges[1:]), rV, V)
Ltot = lambda x: sum(10 ** np.interp(x, g[0], g[1]) for g in gr)
lum = Ltot(x15(edges[:-1]))
lum = np.maximum.accumulate(lum)                # L(r) at the inner edges, non-decreasing
lbol = Ltot(x15(edges[-1:]))[0]
sol = {'c': 2.37e-3, 'n': 6.93e-4, 'o': 5.73e-3, 'ne': 1.26e-3, 'na': 2.9e-5, 'mg': 7.08e-4, 'si': 6.65e-4,
       's': 3.09e-4, 'ca': 6.4e-5, 'fe': 1.29e-3}
zmet = sum(sol.values())
xh, xhe = 0.633, 0.36564
sc = (1 - zmet) / (xh + xhe)
X = {'h': xh * sc, 'he': xhe * sc, **sol}
EL = ['h', 'he', 'c', 'n', 'o', 'ne', 'na', 'mg', 'si', 's', 'ca', 'fe']
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'input.str'), 'w') as fo:
    fo.write('# 1D spherical\n')
    fo.write('# %8d %10d 1 %d %d\n' % (ncell, 1, len(EL) + 8, len(EL)))
    fo.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(EL) + '\n')
    for i in range(ncell):
        row = [edges[i], edges[i + 1], vl[i], vr[i], mass[i], tc[i], tc[i], lum[i]] + [X[e] for e in EL]
        fo.write(' '.join('%.6e' % y for y in row) + '\n')
with open(os.path.join(out, 'structure.txt'), 'w') as fo:
    fo.write('model A of Dessart+2016 from %s: %d cells, %.3g-%.3g cm, mass %.4f Msun, L_bol %.4e erg/s\n'
             % (os.path.basename(raw), ncell, edges[0], edges[-1], mass.sum() / 1.989e33, lbol))
    fo.write('composition: ' + ' '.join('%s %.4g' % (e, X[e]) for e in EL) + '\n')
print(open(os.path.join(out, 'structure.txt')).read(), end='')
