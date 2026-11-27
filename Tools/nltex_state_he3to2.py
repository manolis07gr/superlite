import sys
# move the He III population of cells c1..c2 into the He II ground state (output.nltex_state); keeps the line layout
atoms, sfile, c1, c2 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
sys.argv = ["x", atoms, sfile]
exec(open(__file__.replace('nltex_state_he3to2.py', 'nltex_state_ions.py')).read().split('steps = [rd(f)')[0])
g = dict(groups); i2 = g['He2'].start; i3 = g['He3'].start
L = open(sfile).read().splitlines(); nst = int(L[1].split()[4]); out = L[:4]; body = L[4:]; i = 0; cell = 0
while i < len(body):
    out.append(body[i]); i += 1; cell += 1
    lines = []; n = 0
    while n < nst: lines.append(body[i].split()); n += len(body[i].split()); i += 1
    if c1 <= cell <= c2:
        flat = [t for l in lines for t in l]; a, b = float(flat[i2]), float(flat[i3])
        flat[i2] = '%.7E' % (a + b); flat[i3] = '%.7E' % 0.0
        k = 0
        for l in lines:
            l[:] = flat[k:k + len(l)]; k += len(l)
    for l in lines: out.append('  ' + '  '.join(l))
open(sfile, "w").write('\n'.join(out) + '\n'); print('cells', c1, c2, ': He III -> He II ground, of', cell)
