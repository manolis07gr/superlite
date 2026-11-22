#!/usr/bin/env python3
"""CMFGEN input of the Ibn dense-shell reference model (Tools/ibn_cds_reference.py).

Writes into <model dir>: MODEL_SPEC (the ions of a SuperLite atom table, with
the same number of levels and the same super-level files), batch_ins.sh (the
atomic-data links: osc, col, photoionization routes A-D, autoionization,
F_TO_S) and a list of the SuperLite atom files of the same ions, so that both
codes run on the same atomic data.  VADAT is written separately
(Runs/ibn_cmfgen/VADAT).

usage: Tools/cmfgen_ibn_setup.py <atom table> <CMFGEN atomic dir> <model dir> ion [ion ...]
  ions as named in the table (he1 he2 c1 ...)
"""
import os
import sys

ABR = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Sk', 16: 'S',
       18: 'Ar', 20: 'Ca', 26: 'Fe', 28: 'Nk'}
ROMAN = ['I', '2', 'III', 'IV', 'V', 'SIX', 'SEV', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI']


def nsuper(path, n):
    """number of super-levels of the first n levels of an F_TO_S file"""
    sl = {}
    k = None
    for line in open(path, errors='replace'):
        if '!Entry number of link' in line:
            k = int(line.split()[0])       # column of the super-level (the name is column 1)
            continue
        if k is None:
            continue
        p = line.split()
        if len(p) < k + 1:
            continue
        try:
            s = int(p[k - 1])
        except ValueError:
            continue
        sl[len(sl) + 1] = s        # levels in order (older files have no level number)
    if not sl:
        raise SystemExit('no levels read from ' + path)
    return len({sl[k] for k in sl if k <= n})


def main():
    if len(sys.argv) < 5:
        raise SystemExit(__doc__)
    table, atomic, out = sys.argv[1:4]
    want = sys.argv[4:]
    rows = {}
    for line in open(table):
        p = line.split()
        if len(p) >= 5 and not line.startswith('#'):
            rows[p[0]] = p
    os.makedirs(out, exist_ok=True)
    spec, links = [], []
    for name in want:
        p = rows[name]
        z, stage, d = int(p[1]), int(p[2]), p[3]
        nl = p[4]
        f2s = None
        if ':' in nl:
            nl, f2s = nl.split(':')
        nf = int(nl)
        ion = ABR[z] + ROMAN[stage - 1]
        dd = os.path.join(atomic, d)
        osc, phot, col = (p[5], p[6], p[7]) if len(p) >= 8 else ('osc_data', 'phot_data_A', 'col_data')
        ns = nsuper(os.path.join(dd, f2s), nf) if f2s else nf
        spec.append('%d,%d       [%s_NSF]' % (ns, nf, ion))
        links.append(' ln -sf  %s/%s  %s_F_OSCDAT' % (dd, osc, ion))
        links.append(' ln -sf  %s/%s  %s_COL_DATA' % (dd, col, ion))
        links.append(' ln -sf  %s/%s  PHOT%s_A' % (dd, phot, ion))
        if len(p) < 8:
            for r in 'BCD':
                if os.path.exists(os.path.join(dd, 'phot_data_' + r)):
                    links.append(' ln -sf  %s/phot_data_%s  PHOT%s_%s' % (dd, r, ion, r))
            if os.path.exists(os.path.join(dd, 'auto_data')):
                links.append(' ln -sf  %s/auto_data  %s_AUTO_DATA' % (dd, ion))
        if f2s:
            links.append(' ln -sf  %s/%s  %s_F_TO_S' % (dd, f2s, ion))
        links.append('')
        print('%-5s %-6s NS %4d NF %4d  %s' % (name, ion, ns, nf, d))
    with open(os.path.join(out, 'MODEL_SPEC'), 'w') as f:
        f.write('80          [ND]                !Number of depth points\n')
        f.write('10          [NC]                !Number of core rays\n')
        f.write('90          [NP]                !Number of impact parameters [ND+NC]\n')
        f.write('3           [NUM_BNDS]          !Use 3 for tridiagonal matrix\n')
        f.write('500000      [NCF_MAX]           !Total number of frequencies\n')
        f.write('500         [MAX_SIM]           !Maximum # of lines whose Doppler profile overlap\n')
        f.write('300000      [NLINE_MAX]         !Maximum # of bound-bound transitions\n')
        f.write('21          [NLF]               !For CMF option only\n\n')
        f.write('\n'.join(spec) + '\n')
    with open(os.path.join(out, 'batch_ins.sh'), 'w') as f:
        f.write('#!/bin/sh\n# atomic data links (Tools/cmfgen_ibn_setup.py)\n')
        A = os.path.join(atomic, 'misc')
        f.write(' ln -sf  %s/two_phot_data.dat  TWO_PHOT_DATA\n' % A)
        f.write(' ln -sf  %s/xray_phot_fits.dat  XRAY_PHOT_FITS\n' % A)
        f.write(' ln -sf  %s/HYD/I/5dec96/hyd_l_data.dat  HYD_L_DATA\n' % atomic)
        f.write(' ln -sf  %s/HYD/I/5dec96/gbf_n_data.dat  GBF_N_DATA\n' % atomic)
        f.write('\n'.join(links) + '\n')
    with open(os.path.join(out, 'superlite_atoms.txt'), 'w') as f:
        f.write(' '.join(n + '.atom' for n in want) + '\n')


if __name__ == '__main__':
    main()
