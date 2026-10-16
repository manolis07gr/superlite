#!/usr/bin/env python3
"""Write the SuperLite spectrum of one iteration as wavelength [A] vs L_lambda [erg/s/A].

usage: superlite_spectrum.py <run dir> [--it N] [-o spectrum.dat]
(N: 1-based iteration, default the last; the flux grid is read from output.flx_grid)
"""
import argparse
import os

import numpy as np


def read(run, it=-1):
    with open(os.path.join(run, 'output.flx_grid')) as fh:
        ng = int(fh.readline().split()[1])
        edges = np.array(fh.readline().split(), float)
    rows = [l for l in open(os.path.join(run, 'output.flx_luminos')) if l.strip()]
    lum = np.array(rows[it].split(), float)
    if len(lum) != ng or len(edges) != ng + 1:
        raise SystemExit('inconsistent flux grid and luminosity rows')
    lam = 0.5 * (edges[1:] + edges[:-1]) * 1e8
    return lam, lum / (np.diff(edges) * 1e8), len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run')
    ap.add_argument('--it', type=int, default=0, help='iteration (1-based), default last')
    ap.add_argument('-o', '--out', default=None)
    a = ap.parse_args()
    lam, ll, n = read(a.run, a.it - 1 if a.it > 0 else -1)
    out = a.out or os.path.join(a.run, 'spectrum.dat')
    np.savetxt(out, np.c_[lam, ll], fmt='%.6e', header=f'wavelength[A] L_lambda[erg/s/A]  (iteration {a.it or n} of {n})')
    print(f'wrote {out}: {len(lam)} bins, L = {np.sum(ll * np.gradient(lam)):.4e} erg/s')


if __name__ == '__main__':
    main()
