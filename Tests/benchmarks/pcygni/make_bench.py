import numpy as np
t = 20*86400.0; vph = 5.0e8; vmax = 1.5e9; n = 60
v = np.linspace(vph, vmax, n+1); x = v*t
rho = 1e-14
lines = ['# 1D spherical', '#%9d%11d%2d%3d%3d' % (n, 1, 1, 9, 1),
         '# x_left  x_right  vx_left  vx_right  mass  rad_temp  temp  lum  h']
L = 1e42
for i in range(n):
    m = rho*4/3*np.pi*(x[i+1]**3 - x[i]**3)
    lines.append('  '.join('%.6e' % q for q in (x[i], x[i+1], v[i], v[i+1], m, 1e4, 1e4, L, 1.0)))
open('input.str', 'w').write('\n'.join(lines) + '\n')
