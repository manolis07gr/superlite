import sys
# set T_e of cells c1..c2 (1-based) of an output.nltex_state to T (the warm start reads it with in_tesolve)
f, c1, c2, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
L = open(f).read().splitlines(); nst = int(L[1].split()[4])
out = L[:4]; body = L[4:]; i = 0; cell = 0
while i < len(body):
    h = body[i].split(); cell += 1
    if c1 <= cell <= c2: h[3] = '%.7E' % T
    out.append('  ' + '  '.join(h)); i += 1
    n = 0
    while n < nst: out.append(body[i]); n += len(body[i].split()); i += 1
open(f, 'w').write('\n'.join(out) + '\n'); print('cells', c1, c2, '-> T', T, 'of', cell)
