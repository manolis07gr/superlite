#!/usr/bin/env python3
"""CMFGEN SN_HYDRO_DATA of a SuperLite input.str (model A of Dessart et al. 2016
at one epoch, Tools/dessart/build_modelA.py): the cell centres, outer point
first, with T (held fixed in CMFGEN: [FIX_T] T, [SN_T_OPT] USE_HYDRO), V
(non-monotonic: Sigma = dlnV/dlnR - 1 from the cell edges), rho and the mass
fractions of the cells, so that both codes run on the same structure.
usage: cmfgen_modelA.py input.str out_dir age_days"""
import sys, os, numpy as np
AMASS = {'h': 1.008, 'he': 4.003, 'c': 12.011, 'n': 14.007, 'o': 15.999, 'ne': 20.180, 'na': 22.990,
         'mg': 24.305, 'si': 28.085, 's': 32.06, 'ca': 40.078, 'fe': 55.845}
CMF = {'h': 'HYD', 'he': 'HE', 'c': 'CARB', 'n': 'NIT', 'o': 'OXY', 'ne': 'NEON', 'na': 'SOD', 'mg': 'MAG',
       'si': 'SIL', 's': 'SUL', 'ca': 'CAL', 'fe': 'IRON'}
src, out, age = sys.argv[1], sys.argv[2], float(sys.argv[3])
names = open(src).read().splitlines()[2].lstrip('#').split()
d = np.loadtxt(src, skiprows=3); col = {n: d[:, i] for i, n in enumerate(names)}
el = [n for n in names[8:]]
r = 0.5 * (col['x_left'] + col['x_right']); v = 0.5 * (col['vx_left'] + col['vx_right']) / 1e5
rho = col['mass'] / (4 / 3 * np.pi * (col['x_right'] ** 3 - col['x_left'] ** 3)); T = col['temp']
dlnv = (np.log(np.abs(col['vx_right'])) - np.log(np.abs(col['vx_left']))) * np.sign(col['vx_right'] - col['vx_left'])
sig = (col['vx_right'] - col['vx_left']) / (col['x_right'] - col['x_left']) * r / (v * 1e5) - 1
x = {e: col[e] / sum(col[k] for k in el) for e in el}
mu = 1 / sum(x[e] / AMASS[e] for e in el); natom = rho / (mu * 1.66053907e-24)
ne = natom * np.clip(T / 2e4, 0.05, 1.0); kap = 0.34 * np.clip(T / 1.5e4, 0.05, 1.0)
o = np.argsort(r)[::-1]
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, 'SN_HYDRO_DATA'), 'w') as f:
    f.write('! model A of Dessart et al. (2016) from %s (Tools/dessart/cmfgen_modelA.py)\n' % os.path.basename(os.path.dirname(os.path.abspath(src))))
    f.write('  Number of data points: %d\n  Number of mass fractions: %d\n  Number of isotopes: 1\n' % (len(r), len(el)))
    f.write('  Time(days) since explosion: %.6f\n\n' % age)
    def vec(t, a):
        f.write(' %s\n' % t)
        for i in range(0, len(a), 8): f.write(''.join('%16.8E' % y for y in a[i:i + 8]) + '\n')
        f.write('\n')
    vec('Radius grid (10^10 cm)', r[o] / 1e10); vec('Velocity (km/s)', v[o]); vec('Sigma (dlnV/dlnR-1)', sig[o])
    vec('Temperature (10^4 K)', T[o] / 1e4); vec('Density (gm/cm^3)', rho[o]); vec('Kappa (cm^2/g)', kap[o])
    vec('Electron density (/cm^3)', ne[o]); vec('Atom density (/cm^3)', natom[o])
    for e in el: vec('%s mass fraction' % CMF[e], x[e][o])
    vec('IRON  56  mass fraction', x['fe'][o])
print('points', len(r), 'r %.3e .. %.3e' % (r.min(), r.max()), 'v %.0f .. %.0f' % (v.min(), v.max()),
      'sigma min %.2f max %.2f' % (sig.min(), sig.max()), 'mass %.2f Msun' % (col['mass'].sum() / 1.989e33))
