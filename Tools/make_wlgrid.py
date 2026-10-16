#!/usr/bin/env python3
"""Write input.wlgrid (opacity groups) and input.fluxwl (spectrum bins) with fine
velocity resolution in the optical and coarse resolution elsewhere.

Each segment is logarithmic (constant velocity width).  Use with
  in_grp_ng = 0              (groups from input.wlgrid)
  in_flx_ndim = -1, 1, 1     (flux bins from input.fluxwl)

Example (100 km/s from 3000 to 10000 A, 1000 km/s from 1 A to 3 um):
  python3 make_wlgrid.py -o run_dir --fine 3000 10000 100 --coarse 1 30000 1000
"""
import argparse
import os

import numpy as np

C_KMS = 299792.458


def segment(lo, hi, dv):
    n = max(1, int(np.ceil(np.log(hi / lo) / np.log1p(dv / C_KMS))))
    return lo * (hi / lo) ** (np.arange(n + 1) / n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--out', default='.')
    ap.add_argument('--fine', nargs=3, type=float, default=[3000.0, 10000.0, 100.0],
                    metavar=('LO_A', 'HI_A', 'DV_KMS'))
    ap.add_argument('--coarse', nargs=3, type=float, default=[1.0, 30000.0, 1000.0],
                    metavar=('LO_A', 'HI_A', 'DV_KMS'))
    a = ap.parse_args()
    flo, fhi, fdv = a.fine
    clo, chi, cdv = a.coarse
    if not (clo < flo < fhi < chi):
        raise SystemExit('need coarse_lo < fine_lo < fine_hi < coarse_hi')
    edges = np.concatenate([segment(clo, flo, cdv)[:-1], segment(flo, fhi, fdv)[:-1], segment(fhi, chi, cdv)])
    edges_cm = edges * 1e-8
    ng = len(edges) - 1
    os.makedirs(a.out, exist_ok=True)
    body = ' '.join(f'{x:.8e}' for x in edges_cm)
    with open(os.path.join(a.out, 'input.wlgrid'), 'w') as fh:       # read: ng, then ng+1 bounds
        fh.write(f'{ng}\n{body}\n')
    with open(os.path.join(a.out, 'input.fluxwl'), 'w') as fh:       # read: 3 header lines, 1 skipped, n, bounds
        fh.write('#wavelength bin boundaries for flux. units: [cm]\n# nbin\n# wlbounds\n#\n')
        fh.write(f'{ng}\n{body}\n')
    print(f'{ng} groups: {clo:g}-{flo:g} A and {fhi:g}-{chi:g} A at {cdv:g} km/s, '
          f'{flo:g}-{fhi:g} A at {fdv:g} km/s')


if __name__ == '__main__':
    main()
