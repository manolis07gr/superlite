#!/usr/bin/env python3
"""StaNdaRT SN Ia toy model (Blondin et al. 2022) -> SuperLite input.str.

The toy model (snia_toy01_2d.dat: zone velocities, masses, composition) is
expanded homologously to time t, 56Ni -> 56Co -> 56Fe decayed to t, and
rebinned (mass conserving) to n cells uniform in velocity.  The 'lum'
column (luminosity at the inner edge of each cell, for in_srctype='lumg')
is the cumulative energy deposition of a reference code (the StaNdaRT
edep_<model>_<code>.txt file, interpolated to t), scaled so that the
emergent luminosity is that code's L_bol at t: deposition plus release of
stored radiation, distributed like the deposition.  temp and rad_temp
start from that code's temperature (phys_<model>_<code>.txt).

usage: standart2superlite.py <model .dat> <data dir> <time [d]> <out dir>
       [--code cmfgen] [--n 100] [--vmin 1000] [--vmax 30000]
Writes input.str and a summary (L_bol, total deposition, masses) to stdout.
"""
import argparse
import numpy as np

MSUN = 1.989e33
DAY = 86400.0
T_NI = 6.075 / np.log(2.0)   # mean lives [d], as in mk_snia_toy_model.py
T_CO = 77.236 / np.log(2.0)


def read_model(fname):
    rows = [l.split() for l in open(fname) if l.strip() and not l.startswith('#')]
    a = np.array(rows, dtype=float)
    cols = dict(vel=1, dm=2, xige0=4, xni0=5, xime=6, xti=7, xco_cb=8,
                xca=16, xs=17, xsi=18, xo=19, xc=20)
    return {k: a[:, j] for k, j in cols.items()}


def read_blocks(fname):
    """phys_/ionfrac_ files: blocks '#TIME: t' with a header line."""
    out, t, rows, hdr = {}, None, [], None
    for l in open(fname):
        if l.startswith('#TIME:'):
            if t is not None:
                out[t] = (hdr, np.array(rows, dtype=float))
            t, rows = float(l.split()[1]), []
        elif l.startswith('#'):
            if t is not None and 'vel_mid' in l:
                hdr = l[1:].split()
        elif l.strip() and t is not None:
            rows.append(l.split())
    if t is not None:
        out[t] = (hdr, np.array(rows, dtype=float))
    return out


def read_edep(fname):
    times, rows = None, []
    for l in open(fname):
        if l.startswith('#TIMES'):
            times = np.array(l.split(':')[1].split(), dtype=float)
        elif not l.startswith('#') and l.strip():
            rows.append(l.split())
    a = np.array(rows, dtype=float)
    return times, a[:, 0], a[:, 1:]


def read_lbol(fname):
    a = np.loadtxt(fname, comments='#')
    return a[:, 0], a[:, 1], a[:, 2]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('model')
    p.add_argument('data')
    p.add_argument('time', type=float)
    p.add_argument('out')
    p.add_argument('--code', default='cmfgen')
    p.add_argument('--n', type=int, default=100)
    p.add_argument('--vmin', type=float, default=1000.0)
    p.add_argument('--vmax', type=float, default=30000.0)
    a = p.parse_args()
    name = 'toy01' if 'toy01' in a.model else 'toy06'
    m = read_model(a.model)
    t = a.time
# composition at t (per zone): 56Ni chain from X_56Ni0; stable Fe/Ni from X_IGE0
# with the toy model's relative (Ni, Fe) = (0.1, 0.9)
    x0 = m['xni0']
    xni = x0 * np.exp(-t / T_NI) + 0.1 * m['xige0']
    xco = x0 * T_CO / (T_CO - T_NI) * (np.exp(-t / T_CO) - np.exp(-t / T_NI))
    xfe = x0 - x0 * np.exp(-t / T_NI) - xco + 0.9 * m['xige0']
    comp = dict(c=m['xc'], o=m['xo'], si=m['xsi'], s=m['xs'], ca=m['xca'],
                ti=m['xti'], fe=xfe, co=xco, ni=xni)
# zone edges in velocity [km/s]: uniform zones of width 2*vel(1)
    dv = 2.0 * m['vel'][0]
    vl = m['vel'] - 0.5 * dv
    vr = m['vel'] + 0.5 * dv
    edges = np.linspace(a.vmin, a.vmax, a.n + 1)
    mass = np.zeros(a.n)
    mfr = {k: np.zeros(a.n) for k in comp}
    dm = m['dm'] * MSUN
    for j in range(len(dm)):
# share of zone j in each cell: by volume (velocity^3), uniform density in a zone
        lo = np.clip(edges[:-1], vl[j], vr[j])
        hi = np.clip(edges[1:], vl[j], vr[j])
        w = (hi**3 - lo**3) / (vr[j]**3 - vl[j]**3)
        mass += w * dm[j]
        for k in comp:
            mfr[k] += w * dm[j] * comp[k][j]
    for k in comp:
        mfr[k] /= np.maximum(mass, 1e-300)
    tot = sum(mfr[k] for k in comp)
    for k in comp:
        mfr[k] /= np.maximum(tot, 1e-300)
# reference code: deposition and L_bol at t, temperature
    te, ve, ed = read_edep(f'{a.data}/edep_{name}_{a.code}.txt')
    it = np.searchsorted(te, t)
    w1 = (t - te[it - 1]) / (te[it] - te[it - 1])
    edt = np.exp((1 - w1) * np.log(np.maximum(ed[:, it - 1], 1e-300)) + w1 * np.log(np.maximum(ed[:, it], 1e-300)))
    tl, lb, eb = read_lbol(f'{a.data}/lbol_edep_{name}_{a.code}.txt')
    lbol = np.interp(t, tl, lb)
    edtot = np.interp(t, tl, eb)
# cumulative deposition at the cell edges: edep [erg/s/cm^3] on the
# reference velocity grid, integrated over r = v t
    r = ve * 1e5 * t * DAY
    rmid = 0.5 * (r[1:] + r[:-1])
    re = np.concatenate([[0.0], rmid, [r[-1] + 0.5 * (r[-1] - r[-2])]])
    shell = 4.0 / 3.0 * np.pi * (re[1:]**3 - re[:-1]**3) * edt
    cum = np.concatenate([[0.0], np.cumsum(shell)])
    redge = edges * 1e5 * t * DAY
    lum = np.interp(redge, re, cum)
    ltot_grid = lum[-1]
    lum = lum * lbol / max(ltot_grid, 1e-300)
    ph = read_blocks(f'{a.data}/phys_{name}_{a.code}.txt')
    tk = min(ph, key=lambda x: abs(x - t))
    hdr, arr = ph[tk]
    vmid = 0.5 * (edges[1:] + edges[:-1])
    temp = np.interp(vmid, arr[:, 0], arr[:, 1])
    tt = t * DAY
    labels = ['c', 'o', 'si', 's', 'ca', 'ti', 'fe', 'co', 'ni']
    with open(f'{a.out}/input.str', 'w') as f:
        f.write('# 1D spherical\n')
        f.write('# %8d %10d %d %d %d\n' % (a.n, 1, 1, 8 + len(labels), len(labels)))
        f.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(labels) + '\n')
        for i in range(a.n):
            v0, v1 = edges[i] * 1e5, edges[i + 1] * 1e5
            row = [v0 * tt, v1 * tt, v0, v1, mass[i], temp[i], temp[i], lum[i]] + [mfr[k][i] for k in labels]
            f.write('  '.join('%.6e' % x for x in row) + '\n')
    linner = np.interp(a.vmin * 1e5 * tt, re, cum)
    print(f'{name} at t = {t} d ({a.code}): L_bol = {lbol:.7e} erg/s, deposition = {edtot:.4e} erg/s '
          f'(on its grid {cum[-1]:.4e}, inside {a.vmin} km/s {linner:.4e}, inside {a.vmax} km/s {ltot_grid:.4e})')
    print(f'mass {mass.sum()/MSUN:.4f} Msun in the grid (model {dm.sum()/MSUN:.4f}); '
          f'56Ni-chain Fe+Co+Ni {sum((mfr[k]*mass).sum() for k in ("fe","co","ni"))/MSUN:.4f} Msun; '
          f'temperature from {hdr[1]} at t = {tk} d')


if __name__ == '__main__':
    main()
