#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""CMFGEN's comoving mean intensity J(nu) (EDDFACTOR) at the radii of cells of
a SuperLite run, for Tools/cell_jtest.f90: one file per cell, cmfj_<cell>.dat
(nu [Hz], J [erg/s/cm^2/Hz/sr]), interpolated in log r between CMFGEN's
depths; prints CMFGEN's fraction of H in the ground state there (POPHYD).

usage: cmfgen_j_cells.py <cmfgen model dir> <superlite run dir> <cell,cell,...> <out dir>

The J of a CMFGEN model with incoherent electron scattering (COH_ES = F) can
differ from the J of its final populations by ~1 % (every Lyman line has the
same J/S - 1 at a given radius); use a model with coherent scattering for
tests that need J/S of thick transitions."""
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cmfgen_model as CM


def numbers(s):
    return [float(x) for x in re.findall(r'[-+]?\d+\.\d+E[-+]\d+', s)]


def ground_fraction(path, nd):
    """n(H I ground)/n(H) per depth from POPHYD (first block: total H)."""
    L = open(path, errors='replace').read().splitlines()
    i = 0
    while not re.match(r'^\s*[-+]?\d+\.\d+E', L[i]):
        i += 1
    v = []
    while len(v) < nd:
        v += numbers(L[i])
        i += 1
    tot = np.array(v[:nd])
    while not re.match(r'\s*Number of (\S+) levels:\s+(\d+)', L[i]):
        i += 1
    nlev = int(re.match(r'\s*Number of (\S+) levels:\s+(\d+)', L[i]).group(2))
    i += 1
    while not re.match(r'^\s*[-+]?\d+\.\d+E', L[i]):
        i += 1
    v = []
    while len(v) < nlev*nd:
        v += numbers(L[i])
        i += 1
    return np.array(v[:nlev*nd]).reshape(nd, nlev)[:, 0]/tot


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    mod, run, out = sys.argv[1], sys.argv[2], sys.argv[4]
    cells = [int(x) for x in sys.argv[3].split(',')]
    os.makedirs(out, exist_ok=True)
    info = open(os.path.join(mod, 'EDDFACTOR_INFO')).read().split('\n')
    nd = int(info[2].split()[0])
    a = np.fromfile(os.path.join(mod, 'EDDFACTOR'), dtype='<f8')
    rec = a[:a.size//(nd+1)*(nd+1)].reshape(-1, nd+1)
    ok = (rec[:, nd] > 1e-4) & (rec[:, nd] < 1e3) & np.all(np.isfinite(rec[:, :nd]), axis=1) & (rec[:, nd-1] > 0)
    rec = rec[ok]
    nu, J = rec[:, nd]*1e15, rec[:, :nd]
    o = np.argsort(nu)
    nu, J = nu[o], J[o]
    keep = np.concatenate([[True], np.diff(nu) > 0])
    nu, J = nu[keep], J[keep]
    _, rv = CM.read_rvtj(os.path.join(mod, 'RVTJ'))
    rc = rv['Radius (10^10 cm)']*1e10
    o2 = np.argsort(rc)
    lr = np.log(rc[o2])
    x1 = ground_fraction(os.path.join(mod, 'POPHYD'), len(rc))
    r = np.loadtxt(os.path.join(run, 'output.nlteh'))[:, 0]
    lJ = np.log(np.maximum(J[:, o2], 1e-300))
    for c in cells:
        r0 = r[c-1]
        Jc = np.exp(np.array([np.interp(np.log(r0), lr, row) for row in lJ]))
        np.savetxt(os.path.join(out, 'cmfj_%d.dat' % c), np.column_stack([nu, Jc]), fmt='%.8e')
        print('cell %4d  r %.4e cm  CMFGEN n(H I ground)/n(H) %.4e' %
              (c, r0, np.exp(np.interp(np.log(r0), lr, np.log(x1[o2])))))


if __name__ == '__main__':
    main()
