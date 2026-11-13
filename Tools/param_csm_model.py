#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""Parameterized interaction model (SN ejecta running into a wind-like CSM)
as a SuperLite input.str with the 'lum' column of Tools/stella2superlite_csm.py
(in_srctype = 'lumg'), for H-poor (Type Ibn: He-rich CSM; Type Icn: C/O-rich
CSM) or H-rich interaction when no hydrodynamical model is at hand.

Layers, inside out, at time t after explosion (homologous ejecta v = r/t):
  ejecta      from r_in to the cold dense shell, rho ~ r^-n, mass --mej
  CDS         the cold dense shell at R_s = v_s t: the CSM swept up to R_s
              (Mdot/v_w (R_s - v_w t)), thickness --dcds R_s, speed v_s
  shock       --nshock cells over --dshock R_s outside the CDS in which the
              speed falls linearly from v_s to v_w and the density from the
              shell's to the wind's (as STELLA's smeared shock, so that
              SuperLite's shock cells, in_shock_dvmin, find it); the shock
              luminosity L_sh = eps/2 Mdot/v_w (v_s - v_w)^3 is injected
              across them (the lum column rises by L_sh)
  CSM         rho = Mdot/(4 pi r^2 v_w) at speed v_w out to --rout R_s
The luminosity entering at r_in (from the deeper ejecta) is --lin (in L_sh).
Temperatures are first guesses (the run finds its own: thermal balance).

Compositions (mass fractions): ibn (He-rich, CNO-processed: N up, C and O
down; H free unless --xh), icn (C/O-rich, He-poor), iin (solar, H-rich);
metals heavier than O at solar values (Asplund et al. 2009).  Ejecta take
the same composition unless --ejcomp is given.

Writes input.str, input.par (from the A4 template of Runs/iin_day20 with the
grid size, L_bol, name and no warm start) into the output directory; the
group grid files (input.wlgrid, input.fluxwl) are copied from the A4 kit.

usage: param_csm_model.py -o <dir> [--comp ibn|icn|iin] [--t 20] [--vs 4000]
         [--vw 1000] [--mdot 0.3] [--mej 0.5] [--eps 0.5] ...
"""
import argparse
import os
import shutil
import sys

import numpy as np

MSUN, YR, DAY = 1.98841e33, 3.15576e7, 86400.0
EL = ['h', 'he', 'c', 'n', 'o', 'ne', 'na', 'mg', 'si', 's', 'ar', 'ca', 'fe', 'ni']
# solar mass fractions (Asplund et al. 2009, X = 0.7381, Y = 0.2485, Z = 0.0134)
SOLAR = dict(h=0.7381, he=0.2485, c=2.37e-3, n=6.93e-4, o=5.73e-3, ne=1.26e-3, na=2.93e-5,
             mg=7.08e-4, si=6.66e-4, s=3.10e-4, ar=7.33e-5, ca=6.41e-5, fe=1.29e-3, ni=7.10e-5)
HEAVY = ['ne', 'na', 'mg', 'si', 's', 'ar', 'ca', 'fe', 'ni']


def composition(name, xh=0.0):
    x = {e: 0.0 for e in EL}
    for e in HEAVY:
        x[e] = SOLAR[e]
    if name == 'iin':
        return dict(SOLAR)
    if name == 'ibn':
        # CNO-processed He-rich material: C + N + O conserved, mostly N
        cno = SOLAR['c'] + SOLAR['n'] + SOLAR['o']
        x['c'], x['n'], x['o'] = 0.02 * cno, 0.93 * cno, 0.05 * cno
        x['h'] = xh
        x['he'] = 1.0 - xh - sum(x[e] for e in EL if e not in ('h', 'he'))
    elif name == 'icn':
        # He-burning products (C/O-rich, some He, Ne from alpha capture)
        x['ne'] = 0.03
        x['h'] = xh
        x['he'] = 0.15
        rest = 1.0 - xh - x['he'] - sum(x[e] for e in HEAVY)
        x['c'], x['o'], x['n'] = 0.45 * rest, 0.55 * rest, 0.0
    else:
        sys.exit('composition: ibn, icn or iin')
    return x


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--name', default=None)
    ap.add_argument('--comp', default='ibn', help='CSM composition: ibn, icn, iin')
    ap.add_argument('--ejcomp', default=None, help='ejecta composition (default: as the CSM)')
    ap.add_argument('--xh', type=float, default=0.0, help='H mass fraction of ibn/icn material')
    ap.add_argument('--t', type=float, default=20.0, help='time since explosion [d]')
    ap.add_argument('--vs', type=float, default=4000.0, help='shock (CDS) speed [km/s]')
    ap.add_argument('--vw', type=float, default=1000.0, help='CSM wind speed [km/s]')
    ap.add_argument('--mdot', type=float, default=0.3, help='CSM mass-loss rate [Msun/yr]')
    ap.add_argument('--eps', type=float, default=0.5, help='radiative efficiency of the shock')
    ap.add_argument('--mej', type=float, default=0.5, help='ejecta mass in the grid [Msun]')
    ap.add_argument('--nej', type=float, default=8.0, help='ejecta density exponent')
    ap.add_argument('--rin', type=float, default=0.6, help='inner boundary [R_s]')
    ap.add_argument('--rout', type=float, default=4.0, help='outer boundary [R_s]')
    ap.add_argument('--dcds', type=float, default=0.01, help='CDS thickness [R_s]')
    ap.add_argument('--dshock', type=float, default=0.15, help='width of the smeared shock [R_s]')
    ap.add_argument('--nej-cells', type=int, default=16)
    ap.add_argument('--ncds', type=int, default=4)
    ap.add_argument('--nshock', type=int, default=8)
    ap.add_argument('--ncsm', type=int, default=50)
    ap.add_argument('--lin', type=float, default=0.3, help='luminosity entering at r_in [L_sh]')
    ap.add_argument('--kit', default=None, help='template directory (default: Runs/iin_day20/A4)')
    a = ap.parse_args()

    src = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    kit = a.kit or os.path.join(src, 'Runs', 'iin_day20', 'A4')
    t = a.t * DAY
    vs, vw = a.vs * 1e5, a.vw * 1e5
    mdot = a.mdot * MSUN / YR
    rs = vs * t
    q = mdot / vw                                   # g/cm: rho = q/(4 pi r^2)
    lsh = 0.5 * a.eps * q * (vs - vw) ** 3
    lin = a.lin * lsh
    mcds = q * max(rs - vw * t, 0.0)                # the CSM swept up to R_s
    xcsm = composition(a.comp, a.xh)
    xej = composition(a.ejcomp or a.comp, a.xh)

    cells = []   # (r1, r2, v1, v2, mass, T, L_inner, x)
    # ejecta
    r_in, r_cds = a.rin * rs, rs * (1.0 - a.dcds)
    re = np.linspace(r_in, r_cds, a.nej_cells + 1)
    shells = (re[1:] ** 3 - re[:-1] ** 3) * (0.5 * (re[1:] + re[:-1])) ** (-a.nej)
    me = a.mej * MSUN * shells / shells.sum()
    for k in range(a.nej_cells):
        tk = max(1.2e4 * (re[k] / r_cds) ** -0.5, 6e3)
        cells.append((re[k], re[k + 1], re[k] / t, re[k + 1] / t, me[k], tk, lin, xej))
    # cold dense shell (moving at v_s; the reverse-shocked ejecta join it)
    rc = np.linspace(r_cds, rs, a.ncds + 1)
    for k in range(a.ncds):
        cells.append((rc[k], rc[k + 1], vs, vs, mcds / a.ncds, 7e3, lin, xcsm))
    # smeared shock: speed v_s -> v_w, density from the shell's to the wind's
    rh = np.linspace(rs, rs * (1.0 + a.dshock), a.nshock + 1)
    vh = np.linspace(vs, vw, a.nshock + 1)
    rho_cds = mcds / (4.0 / 3.0 * np.pi * (rs ** 3 - r_cds ** 3))
    for k in range(a.nshock):
        rm = 0.5 * (rh[k] + rh[k + 1])
        rho_w = q / (4 * np.pi * rm ** 2)
        f = (k + 0.5) / a.nshock
        rho = np.exp((1 - f) * np.log(rho_cds) + f * np.log(rho_w)) if f < 0.5 else rho_w * (1 + 4 * (1 - f))
        m = rho * 4.0 / 3.0 * np.pi * (rh[k + 1] ** 3 - rh[k] ** 3)
        cells.append((rh[k], rh[k + 1], vh[k], vh[k + 1], m, 1e4, lin + lsh * k / a.nshock, xcsm))
    # unshocked CSM wind
    rw = np.geomspace(rh[-1], a.rout * rs, a.ncsm + 1)
    for k in range(a.ncsm):
        m = q * (rw[k + 1] - rw[k])
        tk = max(2e4 * (rw[k] / rh[-1]) ** -0.5, 6e3)
        cells.append((rw[k], rw[k + 1], vw, vw, m, tk, lin + lsh, xcsm))

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, 'input.str'), 'w') as fh:
        fh.write('# 1D spherical\n')
        fh.write('# %8d %10d 1 %d %d\n' % (len(cells), 1, len(EL) + 8, len(EL)))
        fh.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(EL) + '\n')
        for r1, r2, v1, v2, m, tk, lk, x in cells:
            fh.write(' '.join('%.6e' % v for v in [r1, r2, v1, v2, m, tk, tk, lk] + [x[e] for e in EL]) + '\n')
    ltot = lin + lsh
    name = a.name or '%s_d%g' % (a.comp, a.t)
    par = open(os.path.join(kit, 'input.par')).read().splitlines()
    out = []
    for l in par:
        s = l.strip().lower()
        if s.startswith('in_name'):
            l = " in_name = '%s'" % name
        elif s.startswith('in_comment'):
            l = (" in_comment = 'parameterized %s interaction (Tools/param_csm_model.py): t %g d, v_s %g, "
                 "v_w %g km/s, Mdot %g Msun/yr, eps %g, M_ej %g Msun'" % (a.comp, a.t, a.vs, a.vw, a.mdot, a.eps, a.mej))
        elif s.startswith('in_ndim'):
            l = ' in_ndim = %d, 1, 1' % len(cells)
        elif s.startswith('in_l_bol'):
            l = ' in_L_bol = %.4e' % ltot
        elif s.startswith('in_nltex_init'):
            l = " in_nltex_init = 'none'"
        out.append(l)
    open(os.path.join(a.out, 'input.par'), 'w').write('\n'.join(out) + '\n')
    for f in ('input.wlgrid', 'input.fluxwl'):
        shutil.copy(os.path.join(kit, f), a.out)
    mtot = sum(c[4] for c in cells) / MSUN
    print('%s: %d cells (%d ejecta, %d CDS, %d shock, %d CSM), R_s %.3e cm, M_CDS %.3f Msun, '
          'M_CSM(grid) %.3f Msun, total %.3f Msun' % (name, len(cells), a.nej_cells, a.ncds, a.nshock, a.ncsm,
                                                      rs, mcds / MSUN, q * (a.rout * rs - rh[-1]) / MSUN, mtot))
    print('L_sh %.3e erg/s (eps %g), L_in %.3e, L_bol %.3e; CSM density at R_s %.3e g/cm^3, n_e(He++) ~ %.2e'
          % (lsh, a.eps, lin, ltot, q / (4 * np.pi * rs ** 2), q / (4 * np.pi * rs ** 2) / 1.6726e-24 / 2.0))


if __name__ == '__main__':
    main()
