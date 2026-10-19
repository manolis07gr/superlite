"""Compare SuperLite's resonance estimator of the H-alpha incident intensity
per cell (output.nlteh, Jres_Ha/B) with the generalized Sobolev solution
(sobolev_nl.py), in units of the core intensity I_c = L B_nu/(4 pi R_c^2
sigma T^4) of the inner-boundary Planck spectrum (T = in_srcmax, input.par).
The estimator sum e c/(nu |Q|)/(4 pi V) is the plain angle mean <I_inc> of
the incident intensity; with in_nltex_jbbeta its terms carry the weight
beta(mu)/beta_bar, and it is <beta I_inc>/beta_bar, which the rates take
(times beta_bar) and which for the pure-scattering line is S.

usage: python3 compare_s.py run_dir nl.S [--every N] [--beta]
"""
import argparse
import os
import re

import numpy as np

SIG = 5.670374419e-5


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run'); ap.add_argument('ref')
    ap.add_argument('--every', type=int, default=5)
    ap.add_argument('--beta', action='store_true', help='the estimator is beta-weighted (in_nltex_jbbeta)')
    ap.add_argument('--rates', action='store_true',
                    help='the J_b of the rates (Jb_Ha/B: the estimator, or the weighted continuum where it has '
                         'fewer than 30 events) against S, all cells')
    a = ap.parse_args()
    par = open(os.path.join(a.run, 'input.par')).read()
    lbol = float(re.search(r'in_L_bol\s*=\s*([\d.eEdD+-]+)', par).group(1).replace('d', 'e'))
    tsrc = float(re.search(r'in_srcmax\s*=\s*([\d.eEdD+-]+)', par).group(1).replace('d', 'e'))
    s = np.loadtxt(os.path.join(a.run, 'input.str'), comments='#')
    rc = s[0, 0]
    h = np.loadtxt(os.path.join(a.run, 'output.nlteh'), usecols=range(22))
    jres, nres = h[:, 21], h[:, 20]
    if a.rates:
        jres = h[:, 17]
        a.beta = True
    # J_b/B with B at the cell temperature: in units of I_c (both Planck at T_src = T_cell)
    tcell = s[:, 6]
    if np.any(np.abs(tcell/tsrc - 1) > 1e-6):
        raise SystemExit('cell temperature differs from the source temperature: scale by B_nu yourself')
    smc = jres*4*np.pi*rc**2*SIG*tsrc**4/lbol
    ref = np.loadtxt(a.ref)
    rr = ref[:, 0]
    sref, sloc, jref = np.empty(len(s)), np.empty(len(s)), np.empty(len(s))
    for i, (r0, r1) in enumerate(s[:, :2]):
        x = np.linspace(r0, r1, 201)
        w = x**2
        for out, col in ((sref, 1), (sloc, 3), (jref, 4)):
            out[i] = np.trapezoid(np.interp(x, rr, ref[:, col])*w, x)/np.trapezoid(w, x)
    target = sref if a.beta else jref
    # noise of the estimator from the number of resonance events
    err = smc/np.sqrt(np.maximum(nres, 1))
    ok = nres >= 100
    if a.rates:
        ok = np.ones(len(nres), bool)
        err = np.full(len(nres), np.inf)
    chi = (smc - target)/np.where(err > 0, err, 1)
    what = '<beta I_inc>/beta_bar' if a.beta else '<I_inc>'
    print(f'{a.run}: estimator against {what}, {ok.sum()} cells with >= 100 events: chi2/dof = '
          f'{np.mean(chi[ok]**2):.2f}, max |dev| {np.abs(chi[ok]).max():.1f} sigma, max |MC/ref - 1| = '
          f'{np.abs(smc/target - 1)[ok].max():.3f}; local Sobolev S off by up to {np.abs(sloc/sref - 1).max():.3f}; '
          f'<I_inc>/S from {np.min(jref/sref):.3f} to {np.max(jref/sref):.3f}')
    print('%5s %10s %9s %9s %9s %9s %8s %7s' % ('cell', 'r', 'J_b MC', '<I_inc>', 'S', 'S_local', 'MC/ref', 'events'))
    for i in range(0, len(s), a.every):
        print('%5d %10.3e %9.4f %9.4f %9.4f %9.4f %8.4f %7d' % (i + 1, 0.5*(s[i, 0] + s[i, 1]), smc[i], jref[i], sref[i],
                                                               sloc[i], smc[i]/target[i], nres[i]))


if __name__ == '__main__':
    main()
