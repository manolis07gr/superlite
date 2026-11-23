#!/usr/bin/env python3
"""A homologous twin of a SuperLite kit model, for SuperLite and CMFGEN alike.

A STELLA snapshot (Tools/stella2superlite_csm.py) is homologous in its
ejecta but not across the cold dense shell, where the velocity falls and
rises again with radius.  CMFGEN (VEL_LAW 11, SN_HYDRO_DATA) needs a
monotonic velocity field.  This script regrids the model onto a monotonic
velocity grid at one age t (v = r/t): each cell's mass, element masses and
added power are spread uniformly over the velocity range the cell spans
and summed per new cell, so the mass in each velocity interval (which sets
the line profiles) is kept, and the radii become r = v t.

The power the kit's lum column adds (STELLA's deposition and the shock
luminosity) is all deposited as heat in the twin: in SuperLite as the lum
column with in_nt_lumg = 0 and no shock cells, in CMFGEN as the
prescribed deposition table EDEP_PROFILE (v [km/s], power per volume
[erg/s/cm^3]), read by a local patch of get_edep_shock_power.f (see
Runs/iin_day20/C14H/README.md).  H, Ar and Ni are left out (no model atoms
in the comparison) and the other mass fractions renormalized.

usage: Tools/stella_twin_homologous.py <kit dir> <out dir> [--t days] [--ncell N] [--nhydro N]
  writes <out>/input.str, SN_HYDRO_DATA, EDEP_PROFILE, structure.txt
"""
import argparse
import os
import re

import numpy as np

MSUN = 1.989e33
AMU = 1.66054e-24
DAY = 86400.0
EL = ['he', 'c', 'n', 'o', 'ne', 'na', 'mg', 'si', 's', 'ca', 'fe']
CMF = {'he': 'HE', 'c': 'CARB', 'n': 'NIT', 'o': 'OXY', 'ne': 'NEON', 'na': 'SOD', 'mg': 'MAG',
       'si': 'SIL', 's': 'SUL', 'ca': 'CAL', 'fe': 'IRON'}
AMASS = {'he': 4.0026, 'c': 12.011, 'n': 14.007, 'o': 15.999, 'ne': 20.180, 'na': 22.990, 'mg': 24.305,
         'si': 28.085, 's': 32.06, 'ca': 40.078, 'fe': 55.845}


def read_kit(d):
    lines = open(os.path.join(d, 'input.str')).read().splitlines()
    labels = lines[2].lstrip('#').split()
    x = np.array([[float(v) for v in l.split()] for l in lines[3:] if l.strip()])
    c = {k: i for i, k in enumerate(labels)}
    lbol = None
    for l in open(os.path.join(d, 'input.par')):
        m = re.match(r'\s*in_L_bol\s*=\s*([0-9.eEdD+-]+)', l, re.I)
        if m:
            lbol = float(m.group(1).replace('d', 'e').replace('D', 'e'))
    return x, c, lbol


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('kit')
    ap.add_argument('out')
    ap.add_argument('--t', type=float, default=None, help='age [d] (default: median r/v of the cells)')
    ap.add_argument('--ncell', type=int, default=95, help='SuperLite cells')
    ap.add_argument('--nhydro', type=int, default=240, help='CMFGEN hydro points')
    a = ap.parse_args()
    x, c, lbol = read_kit(a.kit)
    r1, r2, v1, v2, m = (x[:, c[k]] for k in ('x_left', 'x_right', 'vx_left', 'vx_right', 'mass'))
    lum = x[:, c['lum']]
    temp = x[:, c['temp']]
    p_add = np.append(lum[1:], lbol) - lum               # power each cell adds [erg/s]
    l_in = lum[0]
    xe = np.array([x[:, c[e]] for e in EL]).T
    xe = xe / xe.sum(axis=1, keepdims=True)
    rc, vc = 0.5 * (r1 + r2), 0.5 * (v1 + v2)
    t = a.t * DAY if a.t else np.median(rc / vc)
    lo, hi = np.minimum(v1, v2) / 1e5, np.maximum(v1, v2) / 1e5      # km/s
    vmin, vmax = lo.min(), hi.max()

    def regrid(edges):
        """mass, element masses, power, mass-weighted T per velocity cell"""
        n = len(edges) - 1
        mm, me, pp, tt = np.zeros(n), np.zeros((n, len(EL))), np.zeros(n), np.zeros(n)
        for i in range(len(m)):
            w = hi[i] - lo[i]
            if w <= 0:
                k = min(max(np.searchsorted(edges, lo[i]) - 1, 0), n - 1)
                f = np.zeros(n)
                f[k] = 1.0
            else:
                f = np.clip(np.minimum(edges[1:], hi[i]) - np.maximum(edges[:-1], lo[i]), 0, None) / w
            mm += f * m[i]
            me += np.outer(f * m[i], xe[i])
            pp += f * p_add[i]
            tt += f * m[i] * temp[i]
        return mm, me, pp, tt / np.maximum(mm, 1e-300)

    # SuperLite cells: uniform in v
    ed = np.linspace(vmin, vmax, a.ncell + 1)
    mm, me, pp, tt = regrid(ed)
    ok = mm > 0
    if not ok.all():
        raise SystemExit('empty velocity cells: lower --ncell')
    lumc = l_in + np.concatenate([[0.0], np.cumsum(pp)[:-1]])
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'input.str'), 'w') as f:
        f.write('# 1D spherical\n')
        f.write('# %8d %10d 1 %d %d\n' % (a.ncell, 1, len(EL) + 8, len(EL)))
        f.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(EL) + '\n')
        for i in range(a.ncell):
            row = [ed[i] * 1e5 * t, ed[i + 1] * 1e5 * t, ed[i] * 1e5, ed[i + 1] * 1e5, mm[i], tt[i], tt[i],
                   lumc[i]] + list(me[i] / mm[i])
            f.write(' '.join('%.6e' % y for y in row) + '\n')
    # CMFGEN: finer grid, outer point first; density and composition per point
    eh = np.linspace(vmin, vmax, a.nhydro + 1)
    mh, meh, ph, th = regrid(eh)
    vh = 0.5 * (eh[1:] + eh[:-1])
    vol = 4.0 / 3.0 * np.pi * ((eh[1:] * 1e5 * t) ** 3 - (eh[:-1] * 1e5 * t) ** 3)
    rho = np.maximum(mh / vol, 1e-30)
    xh = meh / np.maximum(mh, 1e-300)[:, None]
    xh[mh <= 0] = xe[-1]
    mu = 1.0 / (xh / np.array([AMASS[e] for e in EL])).sum(axis=1)
    natom = rho / (mu * AMU)
    edep = ph / vol
    o = slice(None, None, -1)
    with open(os.path.join(a.out, 'SN_HYDRO_DATA'), 'w') as f:
        f.write('! homologous twin of %s (Tools/stella_twin_homologous.py)\n' % a.kit)
        f.write('  Number of data points: %d\n' % len(vh))
        f.write('  Number of mass fractions: %d\n' % len(EL))
        f.write('  Number of isotopes: 1\n')
        f.write('  Time(days) since explosion: %.6f\n\n' % (t / DAY))

        def vec(title, v):
            f.write(' %s\n' % title)
            for i in range(0, len(v), 8):
                f.write(''.join('%16.8E' % y for y in v[i:i + 8]) + '\n')
            f.write('\n')
        vec('Radius grid (10^10 cm)', (vh * 1e5 * t / 1e10)[o])
        vec('Velocity (km/s)', vh[o])
        vec('Sigma (dlnV/dlnR-1)', np.zeros(len(vh)))
        vec('Temperature (10^4 K)', (th / 1e4)[o])
        vec('Density (gm/cm^3)', rho[o])
        vec('Kappa (cm^2/g)', np.full(len(vh), 0.1))
        vec('Electron density (/cm^3)', (0.5 * natom)[o])
        vec('Atom density (/cm^3)', natom[o])
        for j, e in enumerate(EL):
            vec('%s mass fraction' % CMF[e], xh[o, j])
        vec('IRON  56  mass fraction', xh[o, EL.index('fe')])
    with open(os.path.join(a.out, 'EDEP_PROFILE'), 'w') as f:
        f.write('! deposited power per volume of the twin (Tools/stella_twin_homologous.py)\n')
        f.write('! total %.6e erg/s; age %.6f d\n' % (pp.sum(), t / DAY))
        f.write('%d\n' % len(vh))
        for v, e in zip(vh, edep):
            f.write('%14.6e %14.6e\n' % (v, e))
    with open(os.path.join(a.out, 'structure.txt'), 'w') as f:
        f.write('twin of %s: age %.4f d, %.0f-%.0f km/s, mass %.4f Msun (kit %.4f), %d SuperLite cells, '
                '%d CMFGEN points\n' % (a.kit, t / DAY, vmin, vmax, mm.sum() / MSUN, m.sum() / MSUN, a.ncell,
                                       a.nhydro))
        f.write('power added %.4e erg/s (all as heat), inner boundary %.4e erg/s, L_bol (in_L_bol) %.4e\n'
                % (pp.sum(), l_in, l_in + pp.sum()))
    print(open(os.path.join(a.out, 'structure.txt')).read(), end='')


if __name__ == '__main__':
    main()
