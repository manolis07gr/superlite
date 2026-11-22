#!/usr/bin/env python3
"""Export a SuperLite 1D model (input.str) as a plain table for other codes.

For a reference calculation of the same model with another code (e.g.
CMFGEN): per cell the radii, velocity, density (volume mean and, with
microclumping, in the clumps), filling factor, temperature, the power the
'lumg' source puts into the cell (the rise of the lum column) and where it
goes, and the mass fractions.  Optionally the electron temperature and
density of a SuperLite run (output.nltex) for comparison.

usage: Tools/export_model_table.py <model dir> <out file> [output.nltex] [--shock-dvmin km/s]

The lum column holds the luminosity entering each cell from below; the
power a cell adds is lum(i+1) - lum(i) (the last cell's is L_bol - lum).
Cells whose velocity drops across them by more than --shock-dvmin
(default 100 km/s, in_shock_dvmin) radiate it as a CIE shock spectrum in
SuperLite; the others deposit it in the gas (with in_nt_sf, in_nt_lumg as
fast electrons).
"""
import re
import sys

import numpy as np

MSUN = 1.989e33


def read_str(path):
    lines = open(path).read().splitlines()
    hdr = lines[1].split()
    ncol, nab = int(hdr[4]), int(hdr[5])
    labels = lines[2].lstrip('#').split()
    data = np.array([[float(v) for v in l.split()] for l in lines[3:] if l.strip()])
    assert data.shape[1] == ncol == len(labels), 'column count'
    return labels, data, labels[ncol - nab:]


def read_lbol(par):
    for l in open(par):
        m = re.match(r'\s*in_L_bol\s*=\s*([0-9.eEdD+-]+)', l, re.I)
        if m:
            return float(m.group(1).replace('d', 'e').replace('D', 'e'))
    return None


def read_nltex_te(path):
    t, ne = [], []
    with open(path) as f:
        f.readline()
        f.readline()
        for l in f:
            if l.startswith('#'):
                break
            v = l.split()
            t.append(float(v[1]))
            ne.append(float(v[2]))
    return np.array(t), np.array(ne)


def main():
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    if len(a) < 2:
        raise SystemExit(__doc__)
    dvmin = 100.0
    if '--shock-dvmin' in sys.argv:
        dvmin = float(sys.argv[sys.argv.index('--shock-dvmin') + 1])
    d, out = a[0], a[1]
    labels, x, ab = read_str(d + '/input.str')
    c = {k: i for i, k in enumerate(labels)}
    r1, r2, v1, v2, m = x[:, c['x_left']], x[:, c['x_right']], x[:, c['vx_left']], x[:, c['vx_right']], x[:, c['mass']]
    temp = x[:, c['temp']]
    fcl = x[:, c['fcl']] if 'fcl' in c else np.ones(len(m))
    lum = x[:, c['lum']] if 'lum' in c else np.zeros(len(m))
    lbol = read_lbol(d + '/input.par')
    lnext = np.append(lum[1:], lbol if lbol else lum[-1])
    dep = lnext - lum
    vol = 4.0 / 3.0 * np.pi * (r2 ** 3 - r1 ** 3)
    rho = m / vol
    shock = (v1 - v2) > dvmin * 1e5
    te = ne = None
    if len(a) > 2:
        te, ne = read_nltex_te(a[2])
    with open(out, 'w') as f:
        f.write('# %s exported by Tools/export_model_table.py from %s/input.str\n' % (out.split('/')[-1], d))
        f.write('# homologous? v/r of first and last cell [1/s]: %.4e %.4e (not homologous: the shell and the '
                'shock have their own speeds)\n' % (0.5 * (v1[0] + v2[0]) / (0.5 * (r1[0] + r2[0])),
                                                   0.5 * (v1[-1] + v2[-1]) / (0.5 * (r1[-1] + r2[-1]))))
        f.write('# L_bol %.4e erg/s; luminosity entering at the inner boundary %.4e erg/s\n' % (lbol or 0, lum[0]))
        f.write('# power added in the gas (deposited): %.4e erg/s; in shock cells (radiated as CIE spectra): %.4e '
                'erg/s\n' % (dep[~shock].sum(), dep[shock].sum()))
        f.write('# total mass %.4f Msun; cells %d\n' % (m.sum() / MSUN, len(m)))
        cols = ['cell', 'r_in[cm]', 'r_out[cm]', 'v_in[km/s]', 'v_out[km/s]', 'mass[g]', 'rho_mean[g/cm3]',
                'fcl', 'rho_clump[g/cm3]', 'T_start[K]', 'P_added[erg/s]', 'P_per_clumpvol[erg/s/cm3]', 'kind']
        if te is not None:
            cols += ['T_e(SuperLite)', 'n_e(SuperLite)']
        cols += ['X_' + e for e in ab]
        f.write('# ' + '  '.join(cols) + '\n')
        for i in range(len(m)):
            row = ['%4d' % (i + 1)] + ['%.6e' % v for v in (r1[i], r2[i])] + \
                  ['%10.2f' % (v1[i] / 1e5), '%10.2f' % (v2[i] / 1e5)] + \
                  ['%.6e' % v for v in (m[i], rho[i], fcl[i], rho[i] / fcl[i], temp[i], dep[i],
                                        dep[i] / (vol[i] * fcl[i]))] + \
                  ['shock' if shock[i] else ('gas' if dep[i] > 0 else '-')]
            if te is not None:
                row += ['%.4e' % te[i], '%.4e' % ne[i]]
            row += ['%.4e' % x[i, c[e]] for e in ab]
            f.write(' '.join(row) + '\n')
    print('wrote %s: %d cells, deposited %.3e erg/s, in shocks %.3e erg/s' % (out, len(m), dep[~shock].sum(),
                                                                          dep[shock].sum()))


if __name__ == '__main__':
    main()
