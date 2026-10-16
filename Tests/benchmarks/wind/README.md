# Benchmark: Sobolev line transfer in non-homologous winds

Pure-scattering H-alpha line (`in_nlteh_tau0`) in a spherical wind from an
absorbing photosphere (`in_innerbc = 'absb'`, radius 10^15 cm) to 10^16 cm,
all continuum opacity off, all energy from the photosphere.  The Sobolev depth
is tau(mu) = tau0 (v/r)/|Q(mu)| with Q(mu) = mu^2 dv/dr + (1-mu^2) v/r, so in a
constant-velocity wind photons cannot leave a resonance radially.  The emergent
profile is compared with the Sobolev formal solution for a two-level atom,
S = beta_c I_c / beta with direction-dependent escape (`sobolev.py`); the
velocity field is the same piecewise-linear one.

```
python3 make_bench.py const --tau0 30 --n2s 22   # or beta, homol
python3 ../../../Tools/make_wlgrid.py -o . --fine 6300 6800 20 --coarse 1 30000 1000
ln -s <superlite>/Data/* .
mpirun -np 1 superlite                            # 4M packets, ~40 s
python3 sobolev.py input.str 30 sobolev.dat 2500 2501
python3 compare.py . sobolev.dat --bin 100 --plot cmp.png
```

Results (2026-10, 50 bins of 100 km/s, errors from packet counts):

| Flow | tau0 | chi^2/dof | largest deviation |
|---|---|---|---|
| v = 2000 km/s | 3 | 1.17 | 3.0 sigma |
| v = 2000 km/s | 30 | 1.34 | 2.9 sigma |
| beta law v = v_inf (1 - 0.9 R_c/r) | 3 | 0.98 | 2.6 sigma |
| beta law | 30 | 0.99 | 2.3 sigma |
| homologous | 3 | 1.20 | 2.9 sigma |
| v = 2000 km/s, isotropic re-emission (stage 1) | 30 | 5.15 | 8.3 sigma |

The last row is the code before line photons were re-emitted with the
direction-dependent escape probability: it fills the absorption trough near
v_inf with photons that escape radially.
