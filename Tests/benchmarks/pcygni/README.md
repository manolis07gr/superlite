# Benchmark: pure-scattering Sobolev P-Cygni profile

Homologous envelope v = r/t (t = 20 d) from v_ph = 5000 km/s to 3 v_ph, all
continuum opacity off, H-alpha a pure-scattering Sobolev line of depth
tau = 3 (`in_nlteh_tau0`), absorbing photosphere (`in_innerbc = 'absb'`),
all energy from the inner boundary.  The emergent profile must match the
Sobolev formal integral with S = W(r) I_ph (Jeffery & Branch 1990).

```
python3 make_bench.py                         # input.str
python3 ../../../Tools/make_wlgrid.py -o .    # input.wlgrid, input.fluxwl
ln -s <superlite>/Data/* .
mpirun -np 1 superlite                        # 262k packets, ~10 s
python3 analytic.py                           # analytic_tau3.dat
```

Result (2026-10): chi^2/dof = 1.21 over 32 bins of 0.2 v_ph with Poisson
errors, largest deviation 2.9 sigma.  This checks the exact resonance
crossings, the Sobolev optical depth, the re-emission at line centre and the
comoving/lab frame transformations together.
