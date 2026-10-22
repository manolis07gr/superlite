#!/usr/bin/env python3
"""StaNdaRT SN Ia toy model (Blondin et al. 2022) -> SuperLite input.str.

Also the CMFGEN delayed-detonation models of StaNdaRT (DDC10_0.976d.dat,
DDC25_1.300d.dat: zone centres, masses, elemental and isotopic mass
fractions at t0 ~ 1 d): all their elements, with the 56Ni -> 56Co -> 56Fe
chain decayed from t0 to t (the other chains are left as at t0), zones
bounded halfway between their centres.

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
       [--code cmfgen] [--n 100] [--vmin 1000] [--vmax 30000] [--edep mc|gray]
--edep mc (gray): the deposition from the decay of the model's 56Ni computed
here (Tools/gamma_deposition.py: Monte Carlo gamma-ray transport with
Compton scattering and photoabsorption, or gray absorption) instead of the
reference code's; L_bol is then the deposition (no scaling).
Writes input.str and a summary (L_bol, total deposition, masses) to stdout.
"""
import os
import sys
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


def read_ddc(fname):
    """CMFGEN model of StaNdaRT (DDC10, DDC25): zone edges [km/s], masses [g],
    t0 [d] and the mass fractions of each element and of 56Ni, 56Co at t0."""
    t0, hdr, rows = None, None, []
    for l in open(fname):
        if l.startswith('# TIME SINCE EXPLOSION'):
            t0 = float(l.split('=')[1].split()[0])
        elif l.startswith('#vel'):
            hdr = l[1:].split()
        elif l.strip() and not l.startswith('#'):
            rows.append(l.split())
    a = np.array(rows, dtype=float)
    col = {h.split('[')[0].lower(): j for j, h in enumerate(hdr)}
    v = a[:, col['vel']]
    vm = 0.5 * (v[1:] + v[:-1])
    edges = np.concatenate([[max(v[0] - (vm[0] - v[0]), 0.0)], vm, [v[-1] + (v[-1] - vm[-1])]])
    els = [h for h in hdr[6:] if not any(ch.isdigit() for ch in h)]
    return dict(t0=t0, edges=edges, dm=a[:, col['dmass']],
                x={e: a[:, col[e]] for e in els}, ni56=a[:, col['ni56']], co56=a[:, col['co56']])


def ddc_comp(m, t):
    """elemental mass fractions at t [d]: the 56Ni chain decayed from t0."""
    dt = t - m['t0']
    n0, c0 = m['ni56'], m['co56']
    ni = n0 * np.exp(-dt / T_NI)
    co = c0 * np.exp(-dt / T_CO) + n0 * T_CO / (T_CO - T_NI) * (np.exp(-dt / T_CO) - np.exp(-dt / T_NI))
    x = dict(m['x'])
    x['ni'] = x['ni'] - n0 + ni
    x['co'] = x['co'] - c0 + co
    x['fe'] = x['fe'] + (n0 + c0) - (ni + co)
    return x


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
    p.add_argument('--edep', choices=['code', 'mc', 'gray'], default='code')
    p.add_argument('--npkt', type=int, default=400000)
    a = p.parse_args()
    base = os.path.basename(a.model).lower()
    name = next((k for k in ('toy01', 'toy06', 'ddc10', 'ddc25') if k in base), 'toy06')
    t = a.time
    if name.startswith('ddc'):
        if a.edep != 'code':
            sys.exit('--edep mc|gray: toy models only')
        m = read_ddc(a.model)
        comp = ddc_comp(m, t)
        vl, vr = m['edges'][:-1], m['edges'][1:]
        dm = m['dm']
    else:
        m = read_model(a.model)
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
        dm = m['dm'] * MSUN
# (no cells inside the model: SuperLite needs mass in every cell)
    if a.vmin < vl[0]:
        print(f'--vmin {a.vmin} is inside the model\'s inner edge; using {vl[0]:.1f} km/s', file=sys.stderr)
        a.vmin = float(vl[0])
    edges = np.linspace(a.vmin, a.vmax, a.n + 1)
    mass = np.zeros(a.n)
    mfr = {k: np.zeros(a.n) for k in comp}
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
# (or the deposition computed here from the model's 56Ni: on the model's
# zones, the reference's velocity points replaced by them)
    if a.edep != 'code':
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import gamma_deposition as gd
        gm = gd.read_model(a.model)
        vz = np.concatenate([vl, [vr[-1]]]) * 1e5
        vb, mb, xb, zb = gd.rebin(vz, gm['dm'], gm['xni0'], 200, gm['zmix'])
        rb = mb / (4.0 / 3.0 * np.pi * ((vb[1:] * t * DAY)**3 - (vb[:-1] * t * DAY)**3))
        if a.edep == 'mc':
            dg, dp = gd.deposition_mc(vb, rb, xb, t, npkt=a.npkt, zmix=zb)
        else:
            dg, dp, _, _ = gd.deposition(vb, rb, xb, t)
        ve = 0.5 * (vb[1:] + vb[:-1]) / 1e5
        edt = dg + dp
        lbol = None
    r = ve * 1e5 * t * DAY
    rmid = 0.5 * (r[1:] + r[:-1])
    re = np.concatenate([[0.0], rmid, [r[-1] + 0.5 * (r[-1] - r[-2])]])
    shell = 4.0 / 3.0 * np.pi * (re[1:]**3 - re[:-1]**3) * edt
    cum = np.concatenate([[0.0], np.cumsum(shell)])
    redge = edges * 1e5 * t * DAY
    lum = np.interp(redge, re, cum)
    ltot_grid = lum[-1]
    if lbol is None:
        lbol = ltot_grid
        edtot = cum[-1]
    lum = lum * lbol / max(ltot_grid, 1e-300)
    ph = read_blocks(f'{a.data}/phys_{name}_{a.code}.txt')
    tk = min(ph, key=lambda x: abs(x - t))
    hdr, arr = ph[tk]
    vmid = 0.5 * (edges[1:] + edges[:-1])
    temp = np.interp(vmid, arr[:, 0], arr[:, 1])
    tt = t * DAY
    labels = [k for k in comp if mfr[k].max() > 0.0]
    with open(f'{a.out}/input.str', 'w') as f:
        f.write('# 1D spherical\n')
        f.write('# %8d %10d %d %d %d\n' % (a.n, 1, 1, 8 + len(labels), len(labels)))
        f.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(labels) + '\n')
        for i in range(a.n):
            v0, v1 = edges[i] * 1e5, edges[i + 1] * 1e5
            row = [v0 * tt, v1 * tt, v0, v1, mass[i], temp[i], temp[i], lum[i]] + [mfr[k][i] for k in labels]
# (lum in full precision: the injected energy of a cell is the difference of
# neighbouring values, 1e-8 of the total and less in the outer ejecta)
            f.write('  '.join(('%.15e' if j == 7 else '%.6e') % x for j, x in enumerate(row)) + '\n')
    linner = np.interp(a.vmin * 1e5 * tt, re, cum)
    print(f'{name} at t = {t} d ({a.code}): L_bol = {lbol:.7e} erg/s, deposition = {edtot:.4e} erg/s '
          f'(on its grid {cum[-1]:.4e}, inside {a.vmin} km/s {linner:.4e}, inside {a.vmax} km/s {ltot_grid:.4e})')
    print(f'mass {mass.sum()/MSUN:.4f} Msun in the grid (model {dm.sum()/MSUN:.4f}); '
          f'56Ni-chain Fe+Co+Ni {sum((mfr[k]*mass).sum() for k in ("fe","co","ni"))/MSUN:.4f} Msun; '
          f'temperature from {hdr[1]} at t = {tk} d')


if __name__ == '__main__':
    main()
