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

## Non-monotonic flows (Rybicki & Hummer 1978)

Where the velocity along a ray is not monotonic, a photon meets the same line
at several points, and the source function at one point depends on the
others (Rybicki & Hummer 1978, ApJ 219, 654).  `sobolev_nl.py` solves the
generalized Sobolev problem of the same pure-scattering line: for every
radius and direction it finds the upstream resonances of the ray (points with
the same projected velocity), carries the core intensity through them
(I -> I e^-tau + S (1 - e^-tau)), and solves the linear system
S beta_bar = <beta(mu) I_inc(mu)> for S(r); then the emergent profile with all
the resonances of each ray.  For a monotonic flow it reduces to `sobolev.py`
(beta law: S identical to 6e-16).  `make_bench.py` has two non-monotonic laws:

* `shell`: homologous ejecta (v = 2000 km/s r/R_s, R_s = 3 R_c), a linear fall
  to 200 km/s over 0.1 R_s, then a constant-velocity circumstellar wind (the
  structure of a IIn: the CSM's narrow P-Cygni absorption against the broad
  emission of the ejecta);
* `accdec`: an acceleration to 2000 km/s at 3 R_c and a linear fall to
  800 km/s at 10 R_c.

```
python3 make_bench.py shell --tau0 3 --n2s 22    # or accdec
python3 ../../../Tools/make_wlgrid.py -o . --fine 6300 6800 20 --coarse 1 30000 1000
ln -s <superlite>/Data/* .
mpirun -np 1 superlite
python3 sobolev_nl.py input.str 3 nl 2500 251    # nl.S, nl.dat (~1.5 min)
python3 compare.py . nl.dat --bin 100
python3 compare_s.py . nl.S                      # the incident intensity per cell
```

Emergent profiles (2026-10, 4M packets, 50 bins of 100 km/s), against the
generalized solution and against the profile of the local Sobolev source
function (the core alone, no nonlocal coupling):

| Flow | tau0 | chi^2/dof | largest deviation | local S: chi^2/dof |
|---|---|---|---|---|
| shell | 3 | 1.14 | 2.9 sigma | 1.79 |
| shell, 40 km/s bins within 600 km/s | 3 | 0.93 | 1.8 sigma | 1.21 |
| accdec | 3 | 1.20 | 2.5 sigma | 27.8 (19.6 sigma) |

The local source function is off by up to 150 % (shell: the CSM lit by the
ejecta) and 1300 % (accdec: the decelerating outer wind lit by the faster
gas inside it); the Monte Carlo has the coupling.  The narrow absorption of
the shell's CSM at +180 km/s is resolved in the 40 km/s bins.

`compare_s.py` checks the resonance estimator of the line's incident
intensity in each cell (output.nlteh, Jres_Ha), which the rates take.  Its
terms e c/(nu |Q|) sum to the plain angle mean <I_inc>; Sobolev's rates need
<beta(mu) I_inc(mu)> = beta_bar S, and the two differ wherever the escape
probability and the incident intensity both vary over directions.  In a wind
of constant velocity the radial escape is small (|Q| -> 0) and the core is
seen along the radial directions: <I_inc>/S grows from 1.09 at the core to
127 at 8 R_c (tau0 = 30).  With `in_nltex_jbbeta` (default) each term has the
weight beta(mu)/beta_bar and the estimator is S (`--beta`).  Cells with at
least 100 events:

| Flow | estimator | against | cells | chi^2/dof | largest deviation |
|---|---|---|---|---|---|
| const, tau0 = 30 | plain | <I_inc> | 13 | 3.93 | 4.0 sigma |
| const, tau0 = 30 | beta-weighted | S | 13 | 0.81 | 1.8 sigma |
| shell | plain | <I_inc> | 60 | 3.83 | 8.6 sigma |
| shell | beta-weighted | S | 60 | 1.46 | 4.7 sigma |
| accdec | plain | <I_inc> | 54 | 1.27 | 2.8 sigma |
| accdec | beta-weighted | S | 54 | 1.10 | 2.8 sigma |

(errors from the number of events; the plain estimator's terms are large
for the rare crossings with small |Q|, so it scatters more than that.)

Where a cell has fewer than 30 events the rates take the J of the line's
group instead, weighted by its shares in eight bins of |mu| (the Monte Carlo
tallies them) with the bin means of the escape probability.  `compare_s.py
--rates` compares the incident intensity the rates take (Jb_Ha, estimator
or weighted group J) with S in every cell.  Constant velocity, tau0 = 30
(100 cells, 0-311 events each): 0.82-1.22 times S out to 9x10^15 cm with
the weighting; before, the plain estimator and the plain group J were 1.06
times S at the core, 2.4-3.6 at 1.2-2x10^15 cm, 6-13 at 2.3-3.5x10^15 and
22-93 beyond 4x10^15 cm.
