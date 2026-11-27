import sys, re, glob, numpy as np
atoms = sys.argv[1]; files = sys.argv[2:]
def rd(f):
    L = open(f).read().splitlines()
    tok = L[2].split()[3:]; Z = [int(t) for t in tok[0::2]]; NS = [int(t) for t in tok[1::2]]
    nst = sum(NS); vals = ' '.join(L[4:]).split(); cells = []; i = 0
    while i < len(vals):
        hdr = [float(v) for v in vals[i:i+5]]; i += 5
        x = np.array(vals[i:i+nst], float); i += nst
        cells.append((hdr, x))
    return Z, NS, cells
sym = {1:'h',2:'he',6:'c',7:'n',8:'o',10:'ne',11:'na',12:'mg',14:'si',16:'s',20:'ca',26:'fe'}
Z, NS, _ = rd(files[0])
groups = []   # (label, slice) of ion stages per element
o = 0
for z, ns in zip(Z, NS):
    k, s = 1, o
    while True:
        f = '%s/%s%d.atom' % (atoms, sym[z], k)
        try: txt = open(f).read().splitlines()
        except FileNotFoundError: break
        i0 = next(i for i, l in enumerate(txt) if l.startswith('nlev')); nl = int(txt[i0].split()[1])
        sl, i = [], i0 + 1
        while len(sl) < nl:
            t = txt[i].split(); i += 1
            if t and t[0][0].isdigit() and not t[0].startswith('#'): sl.append(int(t[-1]))
        n = max(sl)
        groups.append(('%s%d' % (sym[z].capitalize(), k), slice(s, s+n))); s += n; k += 1
    groups.append(('%s%d' % (sym[z].capitalize(), k), slice(s, o+ns))); o += ns
steps = [rd(f)[2] for f in files]
ncell = len(steps[0])
lo = int(sys.argv[0] and 46)
T = np.array([[c[0][3] for c in s] for s in steps])
print('T_e change over the steps, CSM cells: max |dlnT| %.3f' % np.abs(np.diff(np.log(T[:, lo-1:]), axis=0)).max())
rows = []
for l in range(lo-1, ncell):
    for g, sl in groups:
        f = np.array([s[l][1][sl].sum() for s in steps])
        if f.max() < 1e-4: continue
        d = np.diff(np.log(np.maximum(f, 1e-300)))
        rows.append((np.abs(d).max(), l+1, g, f, d))
rows.sort(key=lambda r: -r[0])
alt = sum(1 for r in rows[:100] if np.all(np.sign(r[4][1:]) == -np.sign(r[4][:-1])))
mono = sum(1 for r in rows[:100] if np.all(np.sign(r[4]) == np.sign(r[4][0])))
print('of the 100 largest ion-fraction swings (fraction >1e-4): alternating sign every step %d, one direction %d' % (alt, mono))
for r in rows[:25]:
    print('cell %2d %-5s ' % (r[1], r[2]) + ' '.join('%.2e' % v for v in r[3]))
