#!/usr/bin/env python3
"""STELLA profile -> SuperLite input.str (+ input.par) for the CSM-physics branch.

Differences from stella2superlite.py:
  * writes the luminosity at each cell's inner edge ('lum' column), used by
    in_srctype='lumg' to inject energy where L(r) grows outward (shock
    dissipation, cooling of the ejecta) instead of blackbody 'surfaces';
  * keeps the whole unshocked CSM out to the outer boundary;
  * merges zones that are geometrically thinner than --dr-min (relative) or
    optically thick (STELLA tau > --tau-merge), but never across a velocity
    jump larger than --dv-keep, so shocks stay resolved;
  * truncates the interior at STELLA optical depth --tau-max;
  * writes an input.par with the new physics switched on.

The STELLA columns are read by name from the header line.

Example:
  python3 stella2superlite_csm.py mesa.day020_post_Lbol_max.data -o run_dir --tau-max 50
"""
import argparse
import os
import re
import sys

import numpy as np

MSUN = 1.98841e33

ELEMENTS = ['h', 'he', 'c', 'n', 'o', 'ne', 'na', 'mg', 'si', 's', 'ar', 'ca', 'ti', 'cr', 'fe', 'co', 'ni']

HEADER_MAP = {
    'mass of cell (g)': 'm_cell', 'cell center m (g)': 'm_c', 'cell center R (cm)': 'r_c',
    'cell center v (cm/s)': 'v_c', 'avg density': 'rho', 'radiation pressure': 'prad',
    'avg temperature': 'temp', 'radiation temperature': 'trad', 'avg opacity': 'opac',
    'tau': 'tau', 'outer edge m (g)': 'm_out', 'outer edge r (cm)': 'r_out',
    'luminosity': 'lum', 'n_bar': 'n_bar', 'n_e': 'n_e'}


def read_stella(path):
    with open(path) as fh:
        lines = fh.readlines()
    hdr = {}
    for l in lines[:6]:
        m = re.match(r'\s*(days post max Lbol|zones)\s+([0-9.Ee+-]+)', l)
        if m:
            hdr[m.group(1)] = float(m.group(2))
    # the column-name line: names separated by >= 2 spaces
    names = None
    for i, l in enumerate(lines[:12]):
        if 'mass of cell' in l and 'cell center' in l:
            names = [s.strip() for s in re.split(r'\s{2,}', l.strip()) if s.strip()]
            start = i + 1
            break
    if names is None:
        sys.exit('cannot find the column header in ' + path)
    data = np.loadtxt(lines[start:])
    if data.shape[1] == len(names) + 1:       # leading zone index
        data = data[:, 1:]
    if data.shape[1] != len(names):
        sys.exit(f'{data.shape[1]} columns but {len(names)} names in {path}')
    col = {}
    for n, c in zip(names, data.T):
        col[HEADER_MAP.get(n, n)] = c
    return hdr, col


def element_fractions(col, n):
    """Mass fractions per element (isotopes summed; 56Ni kept as Ni)."""
    out = {e: np.zeros(n) for e in ELEMENTS}
    for k, v in col.items():
        m = re.fullmatch(r'([a-z]+)(\d+)', k)
        if m and m.group(1) in out:
            out[m.group(1)] += v
    tot = sum(out.values())
    for e in ELEMENTS:
        out[e] = out[e] / np.where(tot > 0, tot, 1.0)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('profile')
    ap.add_argument('-o', '--out', default='.')
    ap.add_argument('--tau-max', type=float, default=50.0, help='truncate the interior at this STELLA tau')
    ap.add_argument('--tau-merge', type=float, default=10.0, help='merge optically thick zones above this tau')
    ap.add_argument('--merge-thick', type=int, default=3, help='zones per cell where tau > tau-merge')
    ap.add_argument('--dr-min', type=float, default=2e-3, help='minimum cell width dr/r')
    ap.add_argument('--dv-keep', type=float, default=100.0, help='never merge across a velocity drop [km/s]')
    ap.add_argument('--keep-csm-dl', action='store_true',
                    help="keep the rise of L(r) through the unshocked CSM (by default it is assigned to the "
                         "forward shock: cold CSM cannot radiate it, it is a time-delay effect of the hydro L(r))")
    ap.add_argument('--dv-shock', type=float, default=300.0, help='velocity drop that defines a shock [km/s]')
    ap.add_argument('--n2s', type=int, default=18)
    ap.add_argument('--niter', type=int, default=6)
    a = ap.parse_args()

    hdr, col = read_stella(a.profile)
    n0 = len(col['r_c'])
    r_out = col['r_out']
    r_in = np.r_[col['r_c'][0] - (r_out[0] - col['r_c'][0]), r_out[:-1]]
    # velocities at the edges: linear in r through the cell-centre values
    v_out = np.interp(r_out, col['r_c'], col['v_c'])
    v_in = np.r_[np.interp(r_in[0], col['r_c'], col['v_c']), v_out[:-1]]
    l_out = col['lum']
    l_in = np.r_[l_out[0], l_out[:-1]]
    abund = element_fractions(col, n0)

    # truncate the optically thick interior
    k0 = int(np.argmax(col['tau'] < a.tau_max)) if col['tau'][0] >= a.tau_max else 0
    # group zones into cells
    groups, cur = [], [k0]
    for k in range(k0 + 1, n0):
        dv_drop = (v_out[k - 1] - v_out[k]) / 1e5          # velocity decrease across the boundary
        width = r_out[cur[-1]] - r_in[cur[0]]
        thick = col['tau'][k] > a.tau_merge and len(cur) < a.merge_thick
        thin = width / r_out[cur[-1]] < a.dr_min
        if dv_drop < a.dv_keep and (thick or thin):
            cur.append(k)
        else:
            groups.append(cur)
            cur = [k]
    groups.append(cur)

    m = col['m_cell']
    rows = []
    for g in groups:
        mg = m[g].sum()
        w = m[g] / mg
        row = dict(x_left=r_in[g[0]], x_right=r_out[g[-1]], vx_left=v_in[g[0]], vx_right=v_out[g[-1]],
                   mass=mg, rad_temp=(w * col['trad'][g] ** 4).sum() ** 0.25,
                   temp=(w * col['temp'][g]).sum(), lum=l_in[g[0]])
        for e in ELEMENTS:
            row[e] = (w * abund[e][g]).sum()
        rows.append(row)
    # the forward shock: the outermost contiguous run of cells with v dropping
    # outward by >= dv_shock in total; L beyond it is assigned to the run
    if not a.keep_csm_dl:
        vl = np.array([r['vx_left'] for r in rows]); vr = np.array([r['vx_right'] for r in rows])
        # shock cells: velocity drop > 50 km/s and > 2% (a gently decelerating wind is not a shock)
        dec = (vl - vr > 5e6) & (vl - vr > 0.02 * vl)
        runs, i = [], 0
        while i < len(rows):
            if dec[i]:
                j = i
                while j + 1 < len(rows) and dec[j + 1]:
                    j += 1
                if (vl[i] - vr[j]) / 1e5 >= a.dv_shock:
                    runs.append((i, j))
                i = j + 1
            else:
                i += 1
        if runs:
            i_s, i_e = runs[-1]
            lum = np.array([r['lum'] for r in rows] + [l_out[-1]])
            l0, l1 = lum[i_s], lum[i_e + 1]
            if l1 > l0:
                for k in range(i_s, i_e + 1):
                    rows[k]['lum'] = l0 + (lum[k] - l0) * (l_out[-1] - l0) / (l1 - l0)
            for k in range(i_e + 1, len(rows)):
                rows[k]['lum'] = l_out[-1]
            print(f'forward shock: cells {i_s + 1}-{i_e + 1} (v {vl[i_s]/1e5:.0f} -> {vr[i_e]/1e5:.0f} km/s); '
                  f'L rise through the CSM ({l_out[-1] - l1:.3e} erg/s) assigned to it')
    vx_left = np.array([r['vx_left'] for r in rows])
    if np.any(vx_left < 0):
        sys.exit('negative velocities (infall) are not supported')

    os.makedirs(a.out, exist_ok=True)
    cols = ['x_left', 'x_right', 'vx_left', 'vx_right', 'mass', 'rad_temp', 'temp', 'lum'] + \
           [e for e in ELEMENTS if any(r[e] > 0 for r in rows)]
    nab = len(cols) - 8
    with open(os.path.join(a.out, 'input.str'), 'w') as fh:
        fh.write('# 1D spherical\n')
        fh.write(f'#{len(rows):9d}{1:11d}{1:2d}{len(cols):3d}{nab:3d}\n')
        fh.write('# ' + '  '.join(cols) + '\n')
        for r in rows:
            fh.write('  '.join(f'{r[c]:.6e}' for c in cols) + '\n')
    lbol = l_out[-1]
    par = f"""&inputpars

in_name = '{os.path.basename(os.path.abspath(a.out))}'
in_comment = 'CSM physics branch, {os.path.basename(a.profile)}'

!-- misc
 in_nomp = 1
 in_niter = {a.niter}
 in_io_grabstdout = f
 in_io_opacdump = 'off'
 in_io_profdump = t

!-- grid settings
 in_ndim = {len(rows)}, 1, 1

!-- flux and group settings (1 A - 3 um, 6000 logarithmic groups = 515 km/s)
 in_flx_ndim = 6000, 1, 1
 in_flx_wlmin = 1d-8
 in_flx_wlmax = 30d-5
 in_grp_ng = 6000
 in_grp_wlmin = 1d-8
 in_grp_wlmax = 30d-5

!-- particle numbers and source
 in_src_n2s = {a.n2s}
 in_L_bol = {lbol:.4e}
 in_srctype = 'lumg'
 in_surfsrcmu = 'isot'

!-- transport
 in_trn_tauddmc = 5d0
 in_taulump = 5d0
 in_trn_errorfatal = f
 in_thm_thermal = t

!-- opacities and NLTE hydrogen
 in_nobbopac = f
 in_nobfopac = f
 in_noffopac = f
 in_nothmson = f
 in_nlte = t
 in_nlte_nelem = 1
 in_nlteh = t
 in_nlteh_nlev = 30
 in_nlteh_damp = 0.5d0

/
""".replace('e+', 'd+').replace('e-', 'd-')
    # keep the d-exponent only for in_L_bol
    par = re.sub(r"in_L_bol = ([0-9.]+)d\+(\d+)", r"in_L_bol = \1d+\2", par)
    with open(os.path.join(a.out, 'input.par'), 'w') as fh:
        fh.write(par)
    tau_in = col['tau'][k0]
    print(f'{a.profile}: day {hdr.get("days post max Lbol")}, {n0} STELLA zones -> {len(rows)} cells '
          f'(interior cut at tau={tau_in:.1f}, r={r_in[k0]:.3e} cm); L_in={l_in[k0]:.3e}, L_bol={lbol:.3e} erg/s; '
          f'M={sum(r["mass"] for r in rows)/MSUN:.3f} Msun')


if __name__ == '__main__':
    main()
