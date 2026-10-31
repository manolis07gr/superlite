#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""A saved NLTE state (output.nltex_state) for a larger set of model atoms:
the warm start (in_nltex_init) needs the atoms of the run that wrote it.
Every ion of the old set must be in the new one with the same atom (the
same super-levels); the new set may add ions above an element's old top
ion.  The states of the old ions are copied; the old top state (the ion
above the last old one) goes to the ground state of the first added ion;
the added ions above it and the new top start empty.  The statistical
equilibrium of the first iteration then moves them.  The Monte Carlo
estimators (output.nltex_jmc) belong to the old lines and are not
carried over: the first iteration builds its tables with J = B(T_R).

usage: extend_nltex_state.py <old state> <old atom dir> <new atom dir> <new state>
  each atom dir holds the .atom files of its set (those of in_nltex_atoms,
  and no others); the elements are in the order of Z, as in the run
"""
import glob
import os
import sys

import numpy as np


def atoms(d):
    """per element Z: the ions in order, each (stage, number of states, file)"""
    el = {}
    for n in sorted(glob.glob(os.path.join(d, '*.atom'))):
        z = stage = None
        nsl = 0
        for line in open(n):
            w = line.split()
            if not w or w[0].startswith('#'):
                continue
            if w[0] == 'element':
                z = int(w[1])
            elif w[0] == 'ion':
                stage = int(w[1])
            elif w[0] == 'nlev':
                nlev = int(w[1])
                rows = []
                break
        for line in open(n):
            w = line.split()
            if len(w) >= 5 and w[0].isdigit() and len(rows) < nlev:
                rows.append(int(w[-1]))
        nsl = max(rows)
        el.setdefault(z, []).append((stage, nsl, os.path.basename(n)))
    for z in el:
        el[z].sort()
    return el


def read_state(p):
    lines = [l for l in open(p) if not l.startswith('#')]
    nel, nst, ncell = (int(x) for x in lines[0].split()[3:6])
    zs = [int(x) for x in lines[1].split()[3:]]
    zn = list(zip(zs[0::2], zs[1::2]))
    cells = []
    vals = ' '.join(lines[2:]).split()
    k = 0
    for _ in range(ncell):
        head = [float(x) for x in vals[k:k + 5]]; k += 5
        x = np.array([float(v) for v in vals[k:k + nst]]); k += nst
        cells.append((head, x))
    return zn, cells


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    old, dold, dnew, out = sys.argv[1:]
    ao, an = atoms(dold), atoms(dnew)
    zn, cells = read_state(old)
    if [z for z, _ in zn] != sorted(an) or sorted(ao) != sorted(an):
        sys.exit('the elements of the state, the old and the new atoms differ')
    plan = []
    for z, ns in zn:
        io, inew = ao[z], an[z]
        if sum(i[1] for i in io) + 1 != ns:
            sys.exit(f'Z={z}: the state has {ns} states, the old atoms {sum(i[1] for i in io) + 1}')
        if [i[:2] for i in inew[:len(io)]] != [i[:2] for i in io]:
            sys.exit(f'Z={z}: the new atoms do not start with the old ions (same super-levels)')
        add = inew[len(io):]
        if add and add[0][0] != io[-1][0] + 1:
            sys.exit(f'Z={z}: the added ions do not follow the old top ion')
        plan.append((z, ns, sum(i[1] for i in inew) + 1, add))
    nst = sum(p[2] for p in plan)
    with open(out, 'w') as f:
        f.write('# SuperLite NLTE state populations per cell (in_nltex_prev of a later model)\n')
        f.write(f'nel nst ncell{len(plan):9d}{nst:9d}{len(cells):9d}\n')
        f.write('elements Z ns' + ''.join(f' {z:3d} {n:5d}' for z, _, n, _ in plan) + '\n')
        f.write('# per cell: exterior mass of the cell centre [g], mass [g], density [g/cm^3], '
                'T_e [K], n_e [cm^-3]; then the fractions of the states of each element\n')
        for head, x in cells:
            y = []
            o = 0
            for z, ns, nn, add in plan:
                xe = x[o:o + ns]; o += ns
                if add:
                    ye = np.zeros(nn)
                    ye[:ns - 1] = xe[:-1]
                    ye[ns - 1] = xe[-1]          # old top -> ground state of the first added ion
                else:
                    ye = xe
                y.append(ye)
            y = np.concatenate(y)
            f.write(''.join(f'{v:16.7E}' for v in head) + '\n')
            for i in range(0, len(y), 8):
                f.write(''.join(f'{v:16.7E}' for v in y[i:i + 8]) + '\n')
    added = {z: [i[2] for i in add] for z, _, _, add in plan if add}
    print(f'{len(cells)} cells, {nst} states (was {sum(n for _, n in zn)}); added: {added}')


if __name__ == '__main__':
    main()
