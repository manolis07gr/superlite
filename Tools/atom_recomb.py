#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""
Total recombination coefficient of SuperLite model atoms: the sum over the
levels of the Milne integrals of their tabulated cross-sections (linear in
u = nu/nu0, (u_l/u)^3 beyond the table), as SuperLite's milne() in
nlteatommod.f90.  Dielectronic resonances of the cross-sections
(make_nlte_atoms.sh, DIE) are included: the difference between atoms made
with and without them is the dielectronic rate of the matched transitions.

usage: python3 -I atom_recomb.py [--temps 1e4,2e4,3e4] file.atom ...
"""
import argparse

import numpy as np

H, KB, C, ME = 6.62607015e-27, 1.380649e-16, 2.99792458e10, 9.1093837015e-28
SAHAC = (H * H / (2 * np.pi * ME * KB)) ** 1.5


def read_atom(path):
    """chi, g of the next ion, levels (name, g, E) and cross-section tables"""
    rec = (s.split('!')[0].split() for s in open(path))
    rec = (w for w in rec if w and not w[0].startswith('#'))
    kv, levels, tables = {}, [], []
    for w in rec:
        if w[0] in ('chi', 'gion'):
            kv[w[0]] = float(w[1])
        elif w[0] == 'nlev':
            for _ in range(int(w[1])):
                x = next(rec)
                levels.append((x[1], float(x[2]), float(x[3])))
        elif w[0] == 'nphot':
            for _ in range(int(w[1])):
                n = int(next(rec)[1])
                a = np.array([[float(v) for v in next(rec)[:2]] for _ in range(n)])
                tables.append((a[:, 0], a[:, 1]))
    return kv, levels, tables


def alpha(kv, level, table, temp):
    """Milne recombination coefficient [cm^3/s] of one level"""
    name, g, e = level
    nu0 = (kv['chi'] - e) * C
    u, s = table
    ykt = H * nu0 / (KB * temp)
    y = np.unique(np.concatenate([(u - 1) * ykt, np.linspace(0.0, 50.0, 20001)]))
    y = y[(y >= 0) & (y <= 50)]
    uu = 1 + y / ykt
    sg = np.interp(uu, u, s)
    tl = uu > u[-1]
    sg[tl] = s[-1] * (u[-1] / uu[tl]) ** 3
    f = sg * (nu0 * uu) ** 2 * np.exp(-y)
    integ = np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(y))
    return 4 * np.pi * g / (2 * kv['gion']) * SAHAC * temp ** -1.5 * 2 / C ** 2 * (KB * temp / H) * integ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--temps', default='1e4,2e4,3e4')
    ap.add_argument('atoms', nargs='+')
    a = ap.parse_args()
    temps = [float(t) for t in a.temps.split(',')]
    print('%-40s %s' % ('atom', ' '.join('%10s' % ('T=%.3g' % t) for t in temps)))
    for p in a.atoms:
        kv, levs, tabs = read_atom(p)
        tot = [sum(alpha(kv, lv, tb, t) for lv, tb in zip(levs, tabs)) for t in temps]
        print('%-40s %s' % (p, ' '.join('%10.4e' % x for x in tot)))


if __name__ == '__main__':
    main()
