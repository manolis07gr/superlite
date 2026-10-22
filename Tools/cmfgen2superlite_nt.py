#!/usr/bin/env python3
"""Electron-impact ionization of the subshells of each ion for the
non-thermal electrons of SuperLite (nltentmod): the fits of Arnaud &
Rothenflug (1985),

    sigma = 1e-14/(u I^2) [A (1 - 1/u) + B (1 - 1/u)^2 + C ln u + D ln(u)/u]  cm^2,
    u = E/I,

per subshell of the ground configuration, from CMFGEN's atomic data
(misc/arnaud_rothenflug_*, the most recent unless given), written as
nt_ar.dat in the model-atom directory.  Each line: Z, electrons of the ion,
n and l of the subshell, its threshold I [eV], A, B, C, D (the final state of
CMFGEN's file is dropped).

usage: cmfgen2superlite_nt.py <CMFGEN atomic dir> <out file> [--file name]
"""
import argparse
import glob
import os
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('atomic', help='CMFGEN atomic data directory (containing misc/)')
    ap.add_argument('out', help='output file (nt_ar.dat in the model-atom directory)')
    ap.add_argument('--file', help='the arnaud_rothenflug file of misc/ (default: the most recent)')
    a = ap.parse_args()
    if a.file:
        src = os.path.join(a.atomic, 'misc', a.file)
    else:
        cand = glob.glob(os.path.join(a.atomic, 'misc', 'arnaud_rothenflug_*'))
        if not cand:
            sys.exit('no misc/arnaud_rothenflug_* in ' + a.atomic)
        months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

        def date(f):
            s = os.path.basename(f).split('_')[-1]
            d, m, y = s[:-5], s[-5:-2], s[-2:]
            try:
                return (int(y), months.index(m.lower()), int(d))
            except ValueError:
                return (0, 0, 0)
        src = max(cand, key=date)
    rows = []
    for line in open(src):
        f = line.split()
        if len(f) < 9 or line.lstrip().startswith('!'):
            continue
        try:
            z, ne, n, l = (int(x) for x in f[:4])
            p = [float(x) for x in f[4:9]]
        except ValueError:
            continue
        rows.append((z, ne, n, l, p))
    with open(a.out, 'w') as o:
        o.write('# electron-impact ionization of the subshells (Arnaud & Rothenflug 1985), from %s\n'
                % os.path.basename(src))
        o.write('#  Z  Ne  n  l     I[eV]         A          B          C          D\n')
        for z, ne, n, l, p in rows:
            o.write('%4d %3d %2d %2d %10.4f %10.4f %10.4f %10.4f %10.4f\n' % (z, ne, n, l, *p))
    print('%s: %d subshells of %d ions' % (a.out, len(rows), len({(r[0], r[1]) for r in rows})))


if __name__ == '__main__':
    main()
