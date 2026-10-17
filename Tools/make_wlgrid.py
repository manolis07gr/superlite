#!/usr/bin/env python3
"""Write input.wlgrid (opacity groups) and input.fluxwl (spectrum bins) with fine
velocity resolution in the optical and coarse resolution elsewhere.

Each segment is logarithmic (constant velocity width).  Use with
  in_grp_ng = 0              (groups from input.wlgrid)
  in_flx_ndim = -1, 1, 1     (flux bins from input.fluxwl)

Example (100 km/s from 3000 to 10000 A, 1000 km/s from 1 A to 3 um, 30 000 km/s
from 3 to 300 um):
  python3 make_wlgrid.py -o run_dir --fine 3000 10000 100 --coarse 1 30000 1000 \
      --far 3000000 30000
The far-IR segment holds the fine-structure and forbidden lines of the NLTE model
atoms ([Fe II] 26 um, [Si II] 35 um, [O I] 63 um, [C II] 158 um, [N II] 205 um);
a line emitted outside the grid escapes and is reported as luminosity outside
the flux bins.  --far 0 0 omits it.

Bound-free edges are group boundaries (--edges, default EDGES): a group across
an edge has the cross-section averaged over it, so the photons of the lines
that converge on the edge from the red (the Lyman series above n = 20 within
2 A of 911.8 A) would photoionize.  The boundary sits 30 km/s to the red of
the edge, so that the group to its red has no part above the threshold of the
model atom (energies differ from the table by up to ~1e-5).  The boundary
nearest an edge moves to it when closer than 25 % of the group width, else the
edge splits the group.  --from adds the edges to an existing input.wlgrid.
"""
import argparse
import os

import numpy as np

C_KMS = 299792.458
HC_EV_A = 12398.4198
EDGE_SHIFT = 1e-4             # boundary to the red of the edge (30 km/s)
# ionization energies [eV] of the ground and the main excited levels whose
# edges shape the spectra of winds and SNe: vacuum wavelengths hc/E
EDGES = {'H I 1': 13.598434, 'H I 2': 3.399609, 'H I 3': 1.510937, 'H I 4': 0.849902,
         'He I 1s2': 24.587389, 'He I 2s 3S': 4.767598, 'He I 2s 1S': 3.971548,
         'He I 2p 3P': 3.623, 'He I 2p 1P': 3.369, 'He II 1': 54.417765, 'He II 2': 13.604441,
         'He II 3': 6.046418, 'C I': 11.26030, 'C II': 24.38332, 'C III': 47.8878, 'C IV': 64.4939,
         'N I': 14.53414, 'N II': 29.6013, 'N III': 47.44924, 'N IV': 77.4735, 'O I': 13.61806,
         'O II': 35.11730, 'O III': 54.9355, 'O IV': 77.41353, 'Ne I': 21.56454, 'Ne II': 40.96296,
         'Na I': 5.139076, 'Mg I': 7.646235, 'Mg II': 15.03527, 'Al I': 5.985768, 'Si I': 8.151683,
         'Si II': 16.34585, 'Si III': 33.49300, 'Si IV': 45.14179, 'S I': 10.36001, 'S II': 23.33788,
         'Ca I': 6.113155, 'Ca II': 11.87172, 'Ca III': 50.91316, 'Fe I': 7.902468, 'Fe II': 16.19921,
         'Fe III': 30.651, 'Fe IV': 54.91}


def segment(lo, hi, dv):
    n = max(1, int(np.ceil(np.log(hi / lo) / np.log1p(dv / C_KMS))))
    return lo * (hi / lo) ** (np.arange(n + 1) / n)


def add_edges(edges, wl_edges, frac=0.25):
    """group boundaries [A] with the edges wl_edges [A] as boundaries (a
    boundary that is already an edge stays)"""
    w = list(np.array(edges, dtype=float))
    fixed = {0, len(w) - 1}
    for e in sorted(wl_edges):
        if e <= w[0] or e >= w[-1]:
            continue
        i = int(np.searchsorted(w, e)) - 1      # w[i] <= e < w[i+1]
        dw = w[i + 1] - w[i]
        if e - w[i] < frac * dw and i not in fixed:
            w[i] = e
            fixed.add(i)
        elif w[i + 1] - e < frac * dw and i + 1 not in fixed:
            w[i + 1] = e
            fixed.add(i + 1)
        elif e != w[i]:
            w.insert(i + 1, e)
            fixed = {k + 1 if k > i else k for k in fixed} | {i + 1}
    return np.array(w)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--out', default='.')
    ap.add_argument('--fine', nargs=3, type=float, default=[3000.0, 10000.0, 100.0],
                    metavar=('LO_A', 'HI_A', 'DV_KMS'))
    ap.add_argument('--coarse', nargs=3, type=float, default=[1.0, 30000.0, 1000.0],
                    metavar=('LO_A', 'HI_A', 'DV_KMS'))
    ap.add_argument('--far', nargs=2, type=float, default=[3.0e6, 30000.0],
                    metavar=('HI_A', 'DV_KMS'), help='far-IR segment from the coarse upper end')
    ap.add_argument('--edges', default='default',
                    help="bound-free edges as group boundaries: 'default' (EDGES), 'none', or a "
                         "comma-separated list of wavelengths [A] added to the default")
    ap.add_argument('--from', dest='src', default=None,
                    help='add the edges to this input.wlgrid instead of building the segments')
    a = ap.parse_args()
    flo, fhi, fdv = a.fine
    clo, chi, cdv = a.coarse
    if not (clo < flo < fhi < chi):
        raise SystemExit('need coarse_lo < fine_lo < fine_hi < coarse_hi')
    if a.src:
        tok = open(a.src).read().split()
        edges = np.array([float(x) for x in tok[1:1 + int(tok[0]) + 1]]) * 1e8
    else:
        edges = np.concatenate([segment(clo, flo, cdv)[:-1], segment(flo, fhi, fdv)[:-1], segment(fhi, chi, cdv)])
        if a.far[0] > chi:
            edges = np.concatenate([edges[:-1], segment(chi, a.far[0], a.far[1])])
    if a.edges != 'none':
        wl_edges = [HC_EV_A / ev * (1 + EDGE_SHIFT) for ev in EDGES.values()]
        if a.edges != 'default':
            wl_edges += [float(x) for x in a.edges.split(',')]
        n0 = len(edges)
        edges = add_edges(edges, wl_edges)
        print(f'{len(wl_edges)} bound-free edges as boundaries ({len(edges) - n0} groups added)')
    edges_cm = edges * 1e-8
    ng = len(edges) - 1
    os.makedirs(a.out, exist_ok=True)
    body = ' '.join(f'{x:.8e}' for x in edges_cm)
    with open(os.path.join(a.out, 'input.wlgrid'), 'w') as fh:       # read: ng, then ng+1 bounds
        fh.write(f'{ng}\n{body}\n')
    with open(os.path.join(a.out, 'input.fluxwl'), 'w') as fh:       # read: 3 header lines, 1 skipped, n, bounds
        fh.write('#wavelength bin boundaries for flux. units: [cm]\n# nbin\n# wlbounds\n#\n')
        fh.write(f'{ng}\n{body}\n')
    if a.src:
        print(f'{ng} groups from {a.src}')
    else:
        print(f'{ng} groups: {clo:g}-{flo:g} A and {fhi:g}-{chi:g} A at {cdv:g} km/s, '
              f'{flo:g}-{fhi:g} A at {fdv:g} km/s' +
              (f', {chi:g}-{a.far[0]:g} A at {a.far[1]:g} km/s' if a.far[0] > chi else ''))


if __name__ == '__main__':
    main()
