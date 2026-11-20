#!/usr/bin/env python3
"""Name the element-wide NLTE states of SuperLite diagnostics.

States are numbered over all ions of an element, lowest ion first, each
ion's super-levels in order, the bare top ion last (nltexmod: xoff + soff
+ super-level).  Given the model-atom directory and an element, prints for
each requested state its ion, super-level, member levels (name, g, E in
cm^-1) and the lines among the requested states (f, wavelength).

usage: Tools/atom_states.py <atom dir> <element symbol, e.g. si> <state> [state ...]
"""
import os
import sys


def records(path):
    with open(path) as f:
        for line in f:
            line = line.split('!')[0].strip()
            if not line or line.startswith('#'):
                continue
            yield line


def read_ion(path):
    it = records(path)
    hdr = {}
    for key in ('element', 'ion', 'zion', 'chi', 'gion', 'nlev'):
        r = next(it)
        if not r.startswith(key):
            raise SystemExit(f'{path}: expected {key}, found {r}')
        hdr[key] = r[len(key):].split()
    nlev = int(hdr['nlev'][0])
    levs = []
    for i in range(nlev):
        p = next(it).split()
        name = p[1]
        g, e = float(p[2]), float(p[3])
        sl = int(p[4]) if len(p) > 4 else i + 1
        levs.append((i + 1, name, g, e, sl))
    r = next(it)
    nlin = int(r.split()[-1])
    lines = []
    for _ in range(nlin):
        p = next(it).split()
        lines.append((int(p[0]), int(p[1]), float(p[2]), float(p[4])))
    return int(hdr['ion'][0]), levs, lines


def main():
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    d, el = sys.argv[1], sys.argv[2].lower()
    want = sorted({int(a) for a in sys.argv[3:]})
    ions = []
    for stage in range(1, 40):
        p = os.path.join(d, f'{el}{stage}.atom')
        if os.path.exists(p):
            ions.append(read_ion(p))
    if not ions:
        raise SystemExit(f'no {el}*.atom in {d}')
    off = 0
    smap = {}
    print(f'{el}: ions', ' '.join(str(s) for s, _, _ in ions))
    for stage, levs, lines in ions:
        nsl = max(l[4] for l in levs)
        print(f'  ion {stage}: states {off+1}-{off+nsl} ({len(levs)} levels)')
        for s in range(1, nsl + 1):
            smap[off + s] = (stage, s, [l for l in levs if l[4] == s], lines, levs)
        off += nsl
    print(f'  top ion: state {off+1}')
    for w in want:
        if w not in smap:
            print(f'state {w}: top ion' if w == off + 1 else f'state {w}: none')
            continue
        stage, s, mem, _, _ = smap[w]
        print(f'state {w}: ion {stage}, super-level {s}, {len(mem)} level(s)')
        for i, name, g, e, _ in mem[:12]:
            print(f'    level {i:4d} {name:24s} g={g:5.0f} E={e:12.2f} cm-1')
        if len(mem) > 12:
            print(f'    ... {len(mem)-12} more')
    for a in want:
        for b in want:
            if a >= b or a not in smap or b not in smap or smap[a][0] != smap[b][0]:
                continue
            stage, sa, _, lines, levs = smap[a]
            sb = smap[b][1]
            sl = {l[0]: l[4] for l in levs}
            sel = [x for x in lines if {sl[x[0]], sl[x[1]]} == {sa, sb}]
            sel.sort(key=lambda x: -x[2])
            print(f'lines between states {a} and {b}: {len(sel)}')
            for l, u, f, lam in sel[:8]:
                print(f'    {l:4d} -> {u:4d}  f={f:10.3e}  lambda={lam:12.3f}')


if __name__ == '__main__':
    main()
