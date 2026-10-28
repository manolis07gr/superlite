# A4 and A7 at day 20 with the IIn CSM physics (#83)

Reruns of the two IIn models of `CSM_PHYSICS.md` §4.3 (A4, A7: STELLA
interaction models at day 20) with the current physics and defaults
(local detailed balance, Auger ionization, static Doppler lines in the
slow CSM, level dissolution, Lyman Rayleigh scattering, line overlap,
continuum destruction, adaptive damping), then their emergent spectra as
a formal integral with 12 electron-scattering passes.

| file | content |
|---|---|
| `A4/`, `A7/` | `input.str` (structure, from `Tools/stella2superlite_csm.py`), `input.wlgrid`, `input.fluxwl` (from `Tools/make_wlgrid.py`), `input.par` (template: `@ATOMS@`, `@CIE@`, `@INIT@` filled in by `run.sh`) |
| `A4/init/`, `A7/init/` | `output.nltex_state`: the warm-start populations (A4: `s15_A4`, two iterations of the current physics; A7: `s13_A7`) |
| `table_sn3.txt` | the ions and levels of the model atoms, for `Tools/make_nlte_atoms.sh` |
| `run.sh` | everything below in one script |

## Requirements

- gfortran, MPI (Open MPI), LAPACK/BLAS, CMake
- the CMFGEN atomic data (`atomic_data_21jun23` from the CMFGEN web page;
  not redistributable, so the model atoms are made from it)
- ChiantiPy with the CHIANTI database (`XUVTOP`) for the CIE spectra of
  the radiative shocks; or a `shock_cie.dat` made before (`CIE=...`)
- about 3 GB of memory per MPI rank and 2 GB of disk per model

## Run

```
git clone -b csm-physics https://github.com/manolis07gr/superlite.git
cd superlite
Runs/iin_day20/run.sh /path/to/atomic /path/to/work 4
```

This builds SuperLite in `work/build`, makes the model atoms in
`work/atoms` (`SL_TABLE=table_sn3.txt GFCUT=1e-3`: identical to the atoms
of the saved states), the CIE table in `work/shock_cie.dat`, and then per
model six iterations in `work/A4`, `work/A7` and the formal integral in
`work/A4_fi`, `work/A7_fi` (`output.fi_spectrum`). Optional environment
variables (`SL_BUILD`, `ATOMS`, `CIE`, `JMC`, `NITER`, `MPIRUN`) are
described at the top of `run.sh`.

`JMC` points to the Monte Carlo estimators of the saved states
(`A4/output.nltex_jmc`, `A7/output.nltex_jmc`, 30 MB compressed, given
separately because of their size): with them the first iteration starts
from the radiation field of the saved state; without them it builds its
macro-atom tables with J = B(T_R), and the later iterations recover.

Time: one iteration takes 45–60 min on 2 ranks (most of it the NLTE
solution of the metals), so six iterations are about 5–6 h per model on
2 ranks and less on 4; the formal integral takes minutes.

## Checks

- `grep "NLTE (it" work/A4/run.log`: the median change per iteration
  should fall below a few per cent by the last iterations; if not, run more
  (`NITER`), continuing from `work/<model>/output.nltex_state`.
- The spectra against the earlier runs and the observations:
  `Tools/superlite_spectrum.py` and §4.3 of `CSM_PHYSICS.md`.
