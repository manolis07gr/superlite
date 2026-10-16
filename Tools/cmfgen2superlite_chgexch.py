#!/usr/bin/env python3
"""CMFGEN charge-exchange data (atomic/misc/chg_exch_data.*) -> SuperLite.

Reaction  X(n+1) + Y(m) -> X(n) + Y(m+1)  with rate coefficient
  type 1: a T4^b
  type 2: a T4^b (1 + c exp(d T4))        (Kingdon & Ferland 1996)
  type 3: a T4^b exp(c T4)
[cm^3/s], T4 = T/1e4 K clipped to [tlo, thi] (in 1e4 K).  The reverse rate
follows from detailed balance in SuperLite.

Output, one reaction per line:
  zA sA namesA  zB sB namesB  zC sC namesC  zD sD namesD  type a b c d tlo thi
z: atomic number, s: ionization stage (1 = neutral), names: the level or
term name with its alternatives separated by '/', normalized (no '_',
no J), 'ION' for the bare ion state.

usage: cmfgen2superlite_chgexch.py <chg_exch_data file> <output file>
"""
import re
import sys

ELEM = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sk': 21, 'Tk': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Nk': 28,
        'Ba': 56}
STAGE = {'I': 1, '2': 2, 'III': 3, 'IV': 4, 'V': 5, 'SIX': 6, 'SEV': 7, 'VIII': 8, 'IX': 9,
         'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16, 'XVII': 17}


def species(s):
    """CMFGEN species name (e.g. FeIII, Fe2, HI, Nk2) -> (Z, stage)."""
    for n in (2, 1):
        e, rest = s[:n], s[n:]
        if e in ELEM and rest in STAGE:
            return ELEM[e], STAGE[rest]
    raise ValueError('unknown species ' + s)


def norm(name):
    if name.upper() == 'ION':
        return 'ION'
    return '/'.join(re.sub(r'\[.*?\]', '', a).replace('_', '') for a in name.split('/'))


def tokens(fname):
    for line in open(fname):
        line = line.split('!')[0].strip()
        if line:
            yield line


def main():
    src, out = sys.argv[1], sys.argv[2]
    lines = list(tokens(src))
# header: the line with the number of reactions is the first purely numeric one
    i0 = next(i for i, l in enumerate(lines) if re.fullmatch(r'\d+', l.split()[0]) and len(l.split()) == 1)
    nreac = int(lines[i0])
    rest = [l for l in lines[i0 + 1:] if not re.fullmatch(r'\d{1,2}-\w{3}-\d{4}', l.split()[0])]
    k = 0
    reacs = []
    while len(reacs) < nreac:
        f = rest[k].split()
        k += 1
        if len(f) != 8:
            raise ValueError('reaction line expected: ' + rest[k - 1])
        parts = [(species(f[2 * j]), norm(f[2 * j + 1])) for j in range(4)]
        num = []
        while len(num) < 2 or len(num) < 2 + int(num[1]) + 2:
            num += [float(x.replace('D', 'E').replace('d', 'e')) for x in rest[k].split()]
            k += 1
        typ, nc = int(num[0]), int(num[1])
        coef = num[2:2 + nc] + [0.0] * (4 - nc)
        tlo, thi = num[2 + nc], num[3 + nc]
        reacs.append((parts, typ, coef, tlo, thi))
    with open(out, 'w') as fo:
        fo.write('# SuperLite charge-exchange reactions, converted from CMFGEN %s\n' % src.split('/')[-1])
        fo.write('# X(n+1) + Y(m) -> X(n) + Y(m+1): zA sA A  zB sB B  zC sC C  zD sD D  type a b c d tlo thi\n')
        for parts, typ, coef, tlo, thi in reacs:
            s = '  '.join('%d %d %s' % (z, st, nm) for (z, st), nm in parts)
            fo.write('%s  %d %.4e %.4e %.4e %.4e %.3e %.3e\n' % (s, typ, *coef, tlo, thi))
    print('%d reactions -> %s' % (len(reacs), out))


if __name__ == '__main__':
    main()
