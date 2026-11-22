#!/usr/bin/env python3
"""A homologous Type Ibn dense-shell model for SuperLite and CMFGEN alike.

One continuous structure, written as CMFGEN's SN_HYDRO_DATA (VEL_LAW 11)
and as SuperLite's input.str, so that both codes compute the same model:
homologous expansion at time t (v = r/t); a dense shell (CDS) whose density
is a Gaussian in velocity (centre v_cds, sigma s_cds) holding M_cds of the
he4mix composition, as the CDS models of Dessart, Hillier & Kuncarayakti
2022 (A&A 658, A130); inner ejecta (rho ~ v^-n_ej, M_ej, he4mix) below the
shell; a He-rich CSM outside (rho = Mdot/(4 pi r^2 v_w), the density of a
wind of speed v_w at the radius it now occupies, in homologous motion).
Power P deposited with CMFGEN's shock-power profile exp(-((v - v_cds)/dv)^2)
per unit volume (dv = sqrt(2) s_cds by default), in SuperLite as the lum
column (lumg source, in_nt_lumg = 1: all of it to fast electrons, as
CMFGEN's non-thermal solution of its shock power); no core luminosity
(CMFGEN IB_METH = ZERO_FLUX, SuperLite a reflecting inner boundary).
Ar and Ni are left out of both (no model atoms in the comparison).

usage: Tools/ibn_cds_reference.py <out dir> [options]
  writes <out>/SN_HYDRO_DATA, <out>/input.str, <out>/structure.txt
"""
import argparse
import os

import numpy as np

MSUN = 1.989e33
AMU = 1.66054e-24
DAY = 86400.0
EL = ['he', 'c', 'n', 'o', 'ne', 'na', 'mg', 'si', 's', 'ca', 'fe']
CMF = {'he': 'HE', 'c': 'CARB', 'n': 'NIT', 'o': 'OXY', 'ne': 'NEON', 'na': 'SOD', 'mg': 'MAG',
       'si': 'SIL', 's': 'SUL', 'ca': 'CAL', 'fe': 'IRON'}
AMASS = {'he': 4.0026, 'c': 12.011, 'n': 14.007, 'o': 15.999, 'ne': 20.180, 'na': 22.990, 'mg': 24.305,
         'si': 28.085, 's': 32.06, 'ca': 40.078, 'fe': 55.845}
ZNUM = {'he': 2, 'c': 6, 'n': 7, 'o': 8, 'ne': 10, 'na': 11, 'mg': 12, 'si': 14, 's': 16, 'ca': 20, 'fe': 26}
HE4MIX = {'he': 0.50, 'c': 0.04, 'n': 0.005, 'o': 0.216, 'ne': 0.06, 'na': 0.002, 'mg': 0.04,
          'si': 0.07, 's': 0.03, 'ca': 0.012, 'fe': 0.016}
CSM = {'he': 0.98674, 'c': 1.7586e-4, 'n': 8.1775e-3, 'o': 4.3965e-4, 'ne': 1.2e-3, 'na': 3.3e-5,
       'mg': 6.4e-4, 'si': 7.0e-4, 's': 3.3e-4, 'ca': 6.4e-5, 'fe': 1.3e-3}


def norm(x):
    s = sum(x.values())
    return {k: v / s for k, v in x.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('--t', type=float, default=1e7 / DAY, help='age [d] (default: 2000 km/s at 2e15 cm)')
    ap.add_argument('--vcds', type=float, default=2000.0)
    ap.add_argument('--scds', type=float, default=70.0, help='Gaussian sigma of the shell [km/s]')
    ap.add_argument('--mcds', type=float, default=1.5)
    ap.add_argument('--mej', type=float, default=0.2)
    ap.add_argument('--nej', type=float, default=8.0)
    ap.add_argument('--mdot', type=float, default=0.05, help='CSM [Msun/yr]')
    ap.add_argument('--vw', type=float, default=600.0, help='CSM wind speed giving its density [km/s]')
    ap.add_argument('--vmin', type=float, default=1200.0)
    ap.add_argument('--vmax', type=float, default=6000.0)
    ap.add_argument('--power', type=float, default=6e42)
    ap.add_argument('--dvpow', type=float, default=None, help="CMFGEN DVLOC_SHOCK_POWER [km/s]")
    ap.add_argument('--ncell', type=int, default=80, help='SuperLite cells')
    a = ap.parse_args()
    t = a.t * DAY
    dvp = a.dvpow if a.dvpow else np.sqrt(2.0) * a.scds
    x_sh, x_csm = norm(HE4MIX), norm(CSM)

    def rho_parts(v):
        """density [g/cm^3] of shell, ejecta, CSM at velocity v [km/s] (unnormalized shell/ejecta)"""
        r = v * 1e5 * t
        sh = np.exp(-0.5 * ((v - a.vcds) / a.scds) ** 2)
        ej = np.where(v < a.vcds, (v / a.vcds) ** (-a.nej), 0.0)
        cs = np.where(v > a.vcds, a.mdot * MSUN / 3.15576e7 / (4 * np.pi * r ** 2 * a.vw * 1e5), 0.0)
        return sh, ej, cs

    # normalize the shell and ejecta masses on a fine grid
    vf = np.linspace(a.vmin, a.vmax, 200001)
    rf = vf * 1e5 * t
    sh, ej, cs = rho_parts(vf)
    dm = 4 * np.pi * rf ** 2 * np.gradient(rf)
    ash = a.mcds * MSUN / np.sum(sh * dm)
    aej = a.mej * MSUN / np.sum(ej * dm)

    def rho_x(v):
        sh, ej, cs = rho_parts(v)
        sh, ej = ash * sh, aej * ej
        rho = sh + ej + cs
        x = {e: (sh * x_sh[e] + ej * x_sh[e] + cs * x_csm[e]) / rho for e in EL}
        return rho, x

    pw = lambda v: np.exp(-((v - a.vcds) / dvp) ** 2)
    os.makedirs(a.out, exist_ok=True)

    # --- CMFGEN SN_HYDRO_DATA: fine in the shell, outer point first
    vs = np.unique(np.concatenate([np.geomspace(a.vmin, a.vcds - 6 * a.scds, 60),
                                   np.linspace(a.vcds - 6 * a.scds, a.vcds + 6 * a.scds, 121),
                                   np.geomspace(a.vcds + 6 * a.scds, a.vmax, 80)]))[::-1]
    rs = vs * 1e5 * t
    rho, x = rho_x(vs)
    mu = 1.0 / sum(x[e] / AMASS[e] for e in EL)
    natom = rho / (mu * AMU)
    tguess = np.where(np.abs(vs - a.vcds) < 3 * a.scds, 0.8, 1.0)     # 10^4 K
    ne = natom * np.where(np.abs(vs - a.vcds) < 3 * a.scds, 0.5, 1.0)
    kap = 0.2 * np.where(np.abs(vs - a.vcds) < 3 * a.scds, 0.15, 0.5)   # cm^2/g, for the grid only
    with open(os.path.join(a.out, 'SN_HYDRO_DATA'), 'w') as f:
        f.write('! Ibn dense-shell reference model (Tools/ibn_cds_reference.py), homologous\n')
        f.write('  Number of data points: %d\n' % len(vs))
        f.write('  Number of mass fractions: %d\n' % len(EL))
        f.write('  Number of isotopes: 1\n')
        f.write('  Time(days) since explosion: %.6f\n\n' % a.t)

        def vec(title, v):
            f.write(' %s\n' % title)
            for i in range(0, len(v), 8):
                f.write(''.join('%16.8E' % y for y in v[i:i + 8]) + '\n')
            f.write('\n')
        vec('Radius grid (10^10 cm)', rs / 1e10)
        vec('Velocity (km/s)', vs)
        vec('Sigma (dlnV/dlnR-1)', np.zeros(len(vs)))
        vec('Temperature (10^4 K)', tguess)
        vec('Density (gm/cm^3)', rho)
        vec('Kappa (cm^2/g)', kap)
        vec('Electron density (/cm^3)', ne)
        vec('Atom density (/cm^3)', natom)
        for e in EL:
            vec('%s mass fraction' % CMF[e], x[e])
        vec('IRON  56  mass fraction', x['fe'])

    # --- SuperLite input.str: cells uniform in v away from the shell, fine across it
    nsh = a.ncell // 2
    edges = np.unique(np.concatenate([np.linspace(a.vmin, a.vcds - 5 * a.scds, (a.ncell - nsh) // 3 + 1),
                                      np.linspace(a.vcds - 5 * a.scds, a.vcds + 5 * a.scds, nsh + 1),
                                      np.geomspace(a.vcds + 5 * a.scds, a.vmax,
                                                   a.ncell - nsh - (a.ncell - nsh) // 3 + 1)]))
    cells = []
    psum = 0.0
    pcell = []
    for k in range(len(edges) - 1):
        vv = np.linspace(edges[k], edges[k + 1], 401)
        rr = vv * 1e5 * t
        rh, xx = rho_x(vv)
        w = 4 * np.pi * rr ** 2 * np.gradient(rr)
        m = np.sum(rh * w)
        xm = {e: np.sum(rh * xx[e] * w) / m for e in EL}
        p = np.sum(pw(vv) * w)
        cells.append((rr[0], rr[-1], vv[0] * 1e5, vv[-1] * 1e5, m, xm))
        pcell.append(p)
    pcell = np.array(pcell) * a.power / np.sum(pcell)
    lum = np.concatenate([[0.0], np.cumsum(pcell)[:-1]])
    with open(os.path.join(a.out, 'input.str'), 'w') as f:
        f.write('# 1D spherical\n')
        f.write('# %8d %10d 1 %d %d\n' % (len(cells), 1, len(EL) + 8, len(EL)))
        f.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(EL) + '\n')
        for i, (r1, r2, v1, v2, m, xm) in enumerate(cells):
            vm = 0.5 * (v1 + v2) / 1e5
            tk = 8000.0 if abs(vm - a.vcds) < 3 * a.scds else 1.0e4
            f.write(' '.join('%.6e' % y for y in [r1, r2, v1, v2, m, tk, tk, lum[i]] + [xm[e] for e in EL]) + '\n')
    mtot = sum(c[4] for c in cells) / MSUN
    with open(os.path.join(a.out, 'structure.txt'), 'w') as f:
        f.write('t %.6f d (%.4e s); shell %g Msun at %g km/s, sigma %g km/s; ejecta %g Msun v^-%g; '
                'CSM Mdot %g Msun/yr / v_w %g km/s; %g-%g km/s; power %.3e erg/s, profile '
                'exp(-((v-%g)/%.2f)^2)\n' % (a.t, t, a.mcds, a.vcds, a.scds, a.mej, a.nej, a.mdot, a.vw,
                                             a.vmin, a.vmax, a.power, a.vcds, dvp))
        f.write('SuperLite cells %d, mass %.4f Msun; CMFGEN hydro points %d\n' % (len(cells), mtot, len(vs)))
        f.write('L_bol (in_L_bol) %.4e erg/s\n' % a.power)
    print(open(os.path.join(a.out, 'structure.txt')).read(), end='')


if __name__ == '__main__':
    main()
