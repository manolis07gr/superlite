"""Write input.str and input.par for the Sobolev wind benchmarks.

A pure-scattering H-alpha line (in_nlteh_tau0) in a spherical wind around an
absorbing photosphere of radius RC that emits all the energy (lum constant,
in_srctype='lumg').  No continuum opacity.  Velocity laws:
  const : v = VINF                       (dv/dr = 0: tau(mu) = tau0/(1-mu^2))
  beta  : v = VINF*(1 - 0.9*RC/r)        (accelerating, dv/dr > v/r near RC)
  homol : v = VINF*r/ROUT                (homologous: isotropic tau, S = W I_c)
  shell : v = VINF*r/RS inside RS = 3 RC (homologous ejecta), falling linearly
          to VCSM = 0.1 VINF at 1.1 RS, VCSM outside (slow circumstellar gas):
          non-monotonic, Q < 0 in the decelerating zone, and the ejecta and
          the CSM coupled where |projected velocity| < VCSM (sobolev_nl.py)
  accdec: v = VINF*(0.2 + 0.8 sin(pi/2 (r-RC)/(RM-RC))) up to VINF at RM = 3 RC,
          then a linear fall to 0.4 VINF at ROUT (accelerating, then
          decelerating: several resonances on many rays)
  rh78  : v = VINF (RC/r)^0.5 out to 100 RC (Rybicki & Hummer 1978, alpha =
          1/2: their tau(mu) = tau0/|1 - 1.5 mu^2| is tau0 (v/r)/|Q(mu)|,
          the benchmark depth): a decelerating wind, every point coupled to
          the faster gas inside it
The Sobolev optical depth is tau(mu) = tau0*(v/r)_cell/|Q(mu)|,
Q(mu) = mu^2 dv/dr + (1-mu^2) v/r, as in nltehtrmod (benchmark mode).

usage: python3 make_bench.py {const,beta,homol,shell,accdec,rh78} [--tau0 3] [--n 100] [--n2s 20]
"""
import argparse

import numpy as np

RC, ROUT, VINF = 1e15, 1e16, 2e8

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('law', choices=['const', 'beta', 'homol', 'shell', 'accdec', 'rh78'])
ap.add_argument('--tau0', type=float, default=3.0)
ap.add_argument('--n', type=int, default=100)
ap.add_argument('--n2s', type=int, default=20)
a = ap.parse_args()

r = RC * (ROUT / RC) ** (np.arange(a.n + 1) / a.n)
RS, RM = 3 * RC, 3 * RC


def shell(r):
    v = VINF * r / RS
    f = (r - RS) / (0.1 * RS)
    v = np.where(r > RS, VINF + (0.1 * VINF - VINF) * np.clip(f, 0, 1), v)
    return v


def accdec(r):
    up = VINF * (0.2 + 0.8 * np.sin(0.5 * np.pi * np.clip((r - RC) / (RM - RC), 0, 1)))
    return np.where(r > RM, VINF + (0.4 * VINF - VINF) * (r - RM) / (ROUT - RM), up)


if a.law == 'shell':
    # cell edges at RS and 1.1 RS, so that v is linear in each cell
    n1 = a.n // 2
    r = np.r_[RC * (RS / RC) ** (np.arange(n1) / n1), RS + 0.1 * RS * np.arange(10) / 10,
              1.1 * RS * (ROUT / (1.1 * RS)) ** (np.arange(a.n - n1 - 10 + 1) / (a.n - n1 - 10))]
if a.law == 'rh78':
    ROUT = 100 * RC
    r = RC * (ROUT / RC) ** (np.arange(a.n + 1) / a.n)
v = {'const': lambda r: VINF * np.ones_like(r), 'beta': lambda r: VINF * (1 - 0.9 * RC / r),
     'homol': lambda r: VINF * r / ROUT, 'shell': shell, 'accdec': accdec,
     'rh78': lambda r: VINF * (RC / r) ** 0.5}[a.law](r)
rho, L = 1e-16, 1e42
rows = ['# 1D spherical', '#%9d%11d%2d%3d%3d' % (a.n, 1, 1, 9, 1),
        '# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  h']
for i in range(a.n):
    m = rho * 4 / 3 * np.pi * (r[i + 1] ** 3 - r[i] ** 3)
    rows.append('  '.join('%.8e' % q for q in (r[i], r[i + 1], v[i], v[i + 1], m, 1e4, 1e4, L, 1.0)))
open('input.str', 'w').write('\n'.join(rows) + '\n')
open('input.par', 'w').write(f"""&inputpars
in_name = 'wind_{a.law}'
in_comment = 'pure-scattering Sobolev H-alpha in a {a.law} wind'
 in_nomp = 1
 in_niter = 1
 in_ndim = {a.n}, 1, 1
 in_grp_ng = 0
 in_flx_ndim = -1, 1, 1
 in_src_n2s = {a.n2s}
 in_L_bol = 1d42
 in_srctype = 'lumg'
 in_srcmax = 1d4
 in_surfsrcmu = 'isot'
 in_trn_errorfatal = f
 in_nobbopac = t
 in_nobfopac = t
 in_noffopac = t
 in_nothmson = t
 in_nlte = t
 in_nlte_nelem = 1
 in_nlteh = t
 in_nlteh_tau0 = {a.tau0}d0
 in_innerbc = 'absb'
/
""")
print(f'{a.law} wind: {a.n} cells, RC={RC:.1e} ROUT={ROUT:.1e} cm, v={v[0]/1e5:.0f}..{v[-1]/1e5:.0f} km/s, tau0={a.tau0}')
