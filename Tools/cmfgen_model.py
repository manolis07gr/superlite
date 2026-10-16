#!/usr/bin/env python3
"""Read a converged CMFGEN model (RVTJ, POP<element> files, MOD_SUM) and write
the SuperLite input of the same structure, or the CMFGEN ion fractions for a
comparison.

  cmfgen_model.py input <model dir> <out dir> [--ncell N] [--tau-in 30]
      input.str (cells between CMFGEN radius points, merged to about N cells
      uniform in log r; mass from rho; temp, rad_temp from T; lum = L* at
      every edge: the inner boundary emits L*), and the abundances by mass.
  cmfgen_model.py ions <model dir> [element ...]
      table of r, v, T, n_e and the ion fractions of each element (POP files).

Clumping (volume filling factor f < 1) is not representable in SuperLite;
the converter stops unless --allow-clumping is given (then the mean density
is used).  Run with python3 -I on downloaded data.
"""
import argparse
import os
import re
import sys

import numpy as np

MSUN = 1.989e33
LSUN = 3.828e33
AMU = 1.66053907e-24
# CMFGEN species names -> (symbol, Z, mass [amu])
ELEMENTS = {'HYD': ('h', 1, 1.008), 'HE': ('he', 2, 4.003), 'CARB': ('c', 6, 12.011),
            'NIT': ('n', 7, 14.007), 'OXY': ('o', 8, 15.999), 'NEON': ('ne', 10, 20.180),
            'SOD': ('na', 11, 22.990), 'MAG': ('mg', 12, 24.305), 'ALUM': ('al', 13, 26.982),
            'SIL': ('si', 14, 28.085), 'PHOS': ('p', 15, 30.974), 'SUL': ('s', 16, 32.06),
            'CHL': ('cl', 17, 35.45), 'ARG': ('ar', 18, 39.948), 'POT': ('k', 19, 39.098),
            'CAL': ('ca', 20, 40.078), 'SCAN': ('sc', 21, 44.956), 'TIT': ('ti', 22, 47.867),
            'VAN': ('v', 23, 50.942), 'CHRO': ('cr', 24, 51.996), 'MAN': ('mn', 25, 54.938),
            'IRON': ('fe', 26, 55.845), 'COB': ('co', 27, 58.933), 'NICK': ('ni', 28, 58.693),
            'BAR': ('ba', 56, 137.33)}
POPFILE = {'HYD': 'POPHYD', 'HE': 'POPHE', 'CARB': 'POPCARB', 'NIT': 'POPNIT', 'OXY': 'POPOXY',
           'NEON': 'POPNEON', 'SOD': 'POPSOD', 'MAG': 'POPMAG', 'ALUM': 'POPALUM', 'SIL': 'POPSIL',
           'PHOS': 'POPPHOS', 'SUL': 'POPSUL', 'CHL': 'POPCHL', 'ARG': 'POPARG', 'CAL': 'POPCAL',
           'CHRO': 'POPCHRO', 'MAN': 'POPMAN', 'IRON': 'POPIRON', 'COB': 'POPCOB', 'NICK': 'POPNICK'}


def numbers(s):
    return [float(x.replace('D', 'E')) for x in re.findall(r'[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?', s)]


def read_rvtj(path):
    """Sections of ND values keyed by their title line."""
    lines = open(path, errors='replace').read().splitlines()
    hdr, sec = {}, {}
    nd = None
    i = 0
    while i < len(lines):
        s = lines[i]
        m = re.match(r'\s*(ND|NC|NP|NCF|Mdot\(Msun/yr\)|L\(Lsun\)|H/He abundance):\s+(\S+)', s)
        if m:
            hdr[m.group(1)] = float(m.group(2).replace('D', 'E'))
            if m.group(1) == 'ND':
                nd = int(hdr['ND'])
            i += 1
            continue
        if nd and re.match(r'^\s*[A-Za-z]', s) and not re.match(r'^\s*[-+]?\d', s):
            title = s.strip()
            vals = []
            j = i + 1
            while j < len(lines) and len(vals) < nd:
                t = lines[j]
                if re.match(r'^\s*[A-Za-z]', t):
                    break
                vals += numbers(t)
                j += 1
            if len(vals) >= nd:
                sec[title] = np.array(vals[:nd])
                i = j
                continue
        i += 1
    return hdr, sec


def sec_get(sec, key):
    for k, v in sec.items():
        if k.lower().startswith(key.lower()):
            return v
    raise KeyError(key)


def read_pop(path, nd):
    """Element density and the population of each ion per depth (sum of its
    levels); after the levels of each ion follow ND values of the next ion,
    the last of which is the top ion."""
    lines = open(path, errors='replace').read().splitlines()
    i = 0
    while i < len(lines) and not re.match(r'^\s*[-+]?\d+\.\d+E', lines[i]):
        i += 1
    vals = []
    while len(vals) < nd:
        vals += numbers(lines[i])
        i += 1
    total = np.array(vals[:nd])
    ions = []
    top = None
    while i < len(lines):
        m = re.match(r'\s*Number of (\S+) levels:\s+(\d+)', lines[i])
        if not m:
            i += 1
            continue
        name, nlev = m.group(1), int(m.group(2))
        i += 1
        while i < len(lines) and not re.match(r'^\s*[-+]?\d+\.\d+E', lines[i]):
            i += 1
        vals = []
        while i < len(lines) and len(vals) < (nlev + 1) * nd:
            if re.match(r'^\s*[A-Za-z]', lines[i]):
                break
            vals += numbers(lines[i])
            i += 1
        pops = np.array(vals[:nlev * nd]).reshape(nd, nlev)
        ions.append((name, pops.sum(axis=1)))
        if len(vals) >= (nlev + 1) * nd:
            top = np.array(vals[nlev * nd:(nlev + 1) * nd])
    if top is not None:
        ions.append(('top', top))
    return total, ions


def abundances_vadat(moddir):
    """Mass fractions from VADAT ([XXX/X] keys: positive = number abundance,
    negative = mass fraction; lines starting with '!' are inactive)."""
    num, mass = {}, {}
    for s in open(os.path.join(moddir, 'VADAT'), errors='replace'):
        if s.lstrip().startswith('!'):
            continue
        m = re.match(r'\s*([-+]?[\d.]+(?:[EeDd][-+]?\d+)?)\s+\[(\w+)/X\]', s)
        if not m or m.group(2) not in ELEMENTS:
            continue
        v = float(m.group(1).replace('D', 'E').replace('d', 'e'))
        sym, z, amu = ELEMENTS[m.group(2)]
        if v < 0:
            mass[sym] = -v
        elif v > 0:
            num[sym] = v * amu
    rest = 1.0 - sum(mass.values())
    tot = sum(num.values())
    x = dict(mass)
    for k, v in num.items():
        x[k] = rest * v / tot
    return x


def abundances(moddir):
    """Mass fractions from MOD_SUM ('SPECIES Rel.# Fraction Mass Fraction ...'),
    else from VADAT."""
    if not os.path.exists(os.path.join(moddir, 'MOD_SUM')):
        return abundances_vadat(moddir)
    x = {}
    started = False
    for s in open(os.path.join(moddir, 'MOD_SUM'), errors='replace'):
        if 'SPECIES' in s and 'Mass Fraction' in s:
            started = True
            continue
        if started:
            f = s.split()
            if len(f) >= 3 and f[0] in ELEMENTS:
                try:
                    x[ELEMENTS[f[0]][0]] = float(f[2].replace('D', 'E'))
                except ValueError:
                    pass
            elif f and not f[0] in ELEMENTS and x:
                break
    return x


def cmd_input(a):
    hdr, sec = read_rvtj(os.path.join(a.model, 'RVTJ'))
    r = sec_get(sec, 'Radius') * 1e10
    v = sec_get(sec, 'Velocity') * 1e5
    t = sec_get(sec, 'Temperature') * 1e4
    rho = sec_get(sec, 'Mass Density')
    try:
        clump = sec_get(sec, 'Clumping Factor')
    except KeyError:
        clump = np.ones_like(r)
    if np.any(clump < 0.999) and not a.allow_clumping:
        raise SystemExit('clumped model (min filling factor %.3g): not representable' % clump.min())
    o = np.argsort(r)
    r, v, t, rho = r[o], v[o], t[o], rho[o]
    lstar = hdr['L(Lsun)'] * LSUN
    xab = abundances(a.model)
    labels = sorted(xab, key=lambda k: [e[1] for e in ELEMENTS.values() if e[0] == k][0])
    tot = sum(xab.values())
    # cell edges: the CMFGEN radii, merged (largest log r ratio kept) to at
    # most ncell cells
    edges = r.copy()
    while len(edges) - 1 > a.ncell:
        k = 1 + int(np.argmin(np.diff(np.log(edges))[1:] + np.diff(np.log(edges))[:-1]))
        edges = np.delete(edges, k)
    vedges = np.interp(edges, r, v)
    # mass by integrating rho over the shells (log-log interpolation)
    rr = np.geomspace(edges[0], edges[-1], 4000)
    rhor = np.exp(np.interp(np.log(rr), np.log(r), np.log(rho)))
    dm = 4 * np.pi * rr**2 * rhor
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (dm[1:] + dm[:-1]) * np.diff(rr))])
    mass = np.diff(np.interp(edges, rr, cum))
    rmid = np.sqrt(edges[1:] * edges[:-1])
    tmid = np.interp(rmid, r, t)
    os.makedirs(a.out, exist_ok=True)
    n = len(edges) - 1
    with open(os.path.join(a.out, 'input.str'), 'w') as f:
        f.write('# 1D spherical\n')
        f.write('# %8d %10d %d %d %d\n' % (n, 1, 1, 8 + len(labels), len(labels)))
        f.write('# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  ' + '  '.join(labels) + '\n')
        for i in range(n):
            row = [edges[i], edges[i + 1], vedges[i], vedges[i + 1], mass[i], tmid[i], tmid[i], lstar] + \
                  [xab[k] / tot for k in labels]
            f.write('  '.join('%.6e' % x for x in row) + '\n')
    print(f'{n} cells, r = {edges[0]:.4e} - {edges[-1]:.4e} cm, v = {vedges[0]/1e5:.3g} - {vedges[-1]/1e5:.4g} km/s, '
          f'M = {mass.sum()/MSUN:.4e} Msun, L* = {lstar:.4e} erg/s, T(R*) = {t[0]:.0f} K')
    print('mass fractions: ' + ', '.join(f'{k} {xab[k]/tot:.3e}' for k in labels))
    print(f'in_L_bol = {lstar:.4e};  rad_temp of cell 1 (inner boundary) {tmid[0]:.0f} K')


def cmd_ions(a):
    hdr, sec = read_rvtj(os.path.join(a.model, 'RVTJ'))
    nd = int(hdr['ND'])
    r = sec_get(sec, 'Radius') * 1e10
    v = sec_get(sec, 'Velocity')
    t = sec_get(sec, 'Temperature') * 1e4
    ne = sec_get(sec, 'Electron density')
    els = a.elements or [k for k in POPFILE if os.path.exists(os.path.join(a.model, POPFILE[k]))]
    print('# r[cm] v[km/s] T[K] ne ' + ' '.join(els))
    for el in els:
        total, ions = read_pop(os.path.join(a.model, POPFILE[el]), nd)
        print(f'# {el}: ' + ' '.join(n for n, _ in ions))
        for d in range(nd):
            fr = [p[d] / total[d] for _, p in ions]
            print(f'{r[d]:.4e} {v[d]:.4e} {t[d]:.4e} {ne[d]:.4e} ' + ' '.join(f'{x:.4e}' for x in fr))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('input')
    p.add_argument('model'); p.add_argument('out')
    p.add_argument('--ncell', type=int, default=200)
    p.add_argument('--allow-clumping', action='store_true')
    p = sp.add_parser('ions')
    p.add_argument('model'); p.add_argument('elements', nargs='*')
    a = ap.parse_args()
    {'input': cmd_input, 'ions': cmd_ions}[a.cmd](a)


if __name__ == '__main__':
    main()
