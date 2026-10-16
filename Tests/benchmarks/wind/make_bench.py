"""Write input.str and input.par for the Sobolev wind benchmarks.

A pure-scattering H-alpha line (in_nlteh_tau0) in a spherical wind around an
absorbing photosphere of radius RC that emits all the energy (lum constant,
in_srctype='lumg').  No continuum opacity.  Velocity laws:
  const : v = VINF                       (dv/dr = 0: tau(mu) = tau0/(1-mu^2))
  beta  : v = VINF*(1 - 0.9*RC/r)        (accelerating, dv/dr > v/r near RC)
  homol : v = VINF*r/ROUT                (homologous: isotropic tau, S = W I_c)
The Sobolev optical depth is tau(mu) = tau0*(v/r)_cell/|Q(mu)|,
Q(mu) = mu^2 dv/dr + (1-mu^2) v/r, as in nltehtrmod (benchmark mode).

usage: python3 make_bench.py {const,beta,homol} [--tau0 3] [--n 100] [--n2s 20]
"""
import argparse

import numpy as np

RC, ROUT, VINF = 1e15, 1e16, 2e8

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('law', choices=['const', 'beta', 'homol'])
ap.add_argument('--tau0', type=float, default=3.0)
ap.add_argument('--n', type=int, default=100)
ap.add_argument('--n2s', type=int, default=20)
a = ap.parse_args()

r = RC * (ROUT / RC) ** (np.arange(a.n + 1) / a.n)
v = {'const': VINF * np.ones_like(r), 'beta': VINF * (1 - 0.9 * RC / r), 'homol': VINF * r / ROUT}[a.law]
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
