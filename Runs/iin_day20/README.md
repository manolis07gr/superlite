# A4 and A7 at day 20 with the IIn CSM physics (#83)

Reruns of the two IIn models of `CSM_PHYSICS.md` §4.3 (A4, A7: STELLA
interaction models at day 20) with the current physics and defaults
(local detailed balance, Auger ionization, static Doppler lines in the
slow CSM, level dissolution, Lyman Rayleigh scattering, line overlap,
continuum destruction, adaptive damping), then their emergent spectra as
a formal integral with 12 electron-scattering passes.

The model atoms are sn6 (`table_sn6.txt`, 133 ions): CMFGEN's to its last
ions (Fe XVI, the H-like C, N, O; `table_sn5.txt`) and, built from the
CHIANTI database (`Tools/chianti2superlite_atom.py`), the ions beyond, so
that every element reaches its bare nucleus. A4's shock-heated CSM
(6–8×10⁴ K, iron in Fe XIV–XIX) needs them: with sn3 (44 ions, to Fe VI)
every metal there sat in the top state of its atom; A7 does not
(CSM_PHYSICS.md §4.3).

The overlapping lines are those of the pairs with an H or He line
(`in_nltex_ovlz = 2`: Lyβ and O I 1025.76 Å, He II and O III): with all
903 546 pairs of the slow, dense CSM the metal blends made the NLTE
iteration swing by factors of 10–600 between iterations (A4 median
largest change 672, 21, 411, 18, 60); with the 9 384 pairs of H and He as
without the overlap (A7: 4.2 and 0.70 against 4.2 and 0.74; A4: 11.4,
14.9, 1.9 against 10.5, 14.3, 1.2).
The pairs of lines of one element are kept as well (`in_nltex_ovlsame`:
564 592 pairs in A7, the Fe II, Fe III, ... blends): the element's rate
matrix holds their exchange implicitly and the statistical equilibrium
converges (A7 median largest change 4.2, 0.77, 0.93, 0.52, 0.51, 0.33,
0.27; 1–6 of 68 cells at the iteration limit after the first iterations).
The pairs of two metals are left out: the statistical equilibrium, which
solves one element after the other, does not converge with them (A7
cells with the radiation field held fixed: at the iteration limit, and
solutions from LTE and from the saved state differ by factors of 2 to
10⁴; `Tools/cell_jtest.f90`, `SL_UNIQ`).

| file | content |
|---|---|
| `A4/`, `A7/` | `input.str` (structure, from `Tools/stella2superlite_csm.py`), `input.wlgrid`, `input.fluxwl` (from `Tools/make_wlgrid.py`), `input.par` (template: `@ATOMS@`, `@CIE@`, `@INIT@` filled in by `run.sh`) |
| `A4/init/`, `A7/init/` | `output.nltex_state`: the warm-start populations, of the sn3 atoms (A4: `s15_A4`; A7: `s13_A7`); `run.sh` extends them to the atoms of the run (`Tools/extend_nltex_state.py`: the added ions start in their ground states) |
| `table_sn6.txt` | the model atoms: the CMFGEN ions of `table_sn5.txt` (`Tools/make_nlte_atoms.sh`) and the CHIANTI ions |
| `table_sn3.txt`, `table_sn4.txt` | the smaller atom sets of the earlier runs |
| `run.sh` | everything below in one script |

## Requirements

- gfortran, MPI (Open MPI), LAPACK/BLAS, CMake
- the CMFGEN atomic data (`atomic_data_21jun23` from the CMFGEN web page;
  not redistributable, so the model atoms are made from it)
- the CHIANTI database (`XUVTOP`; v11 tested) for the ions beyond
  CMFGEN's data, and ChiantiPy for the CIE spectra of the radiative
  shocks (or a `shock_cie.dat` made before, `CIE=...`); python3 with numpy
  and scipy
- about 6 GB of memory per MPI rank (6011 states, 280 291 lines) and 3 GB
  of disk per model

## Run

```
git clone -b csm-physics https://github.com/manolis07gr/superlite.git
cd superlite
Runs/iin_day20/run.sh /path/to/atomic /path/to/work 4
```

This builds SuperLite in `work/build`, makes the model atoms in
`work/atoms` (the CMFGEN ions with `SL_TABLE=table_sn5.txt GFCUT=1e-3`,
the others from CHIANTI), the CIE table in `work/shock_cie.dat`, and then per
model six iterations in `work/A4`, `work/A7` and the formal integral in
`work/A4_fi`, `work/A7_fi` (`output.fi_spectrum`). Optional environment
variables (`SL_BUILD`, `ATOMS`, `CIE`, `JMC`, `NITER`, `MPIRUN`) are
described at the top of `run.sh`.

`JMC` points to the Monte Carlo estimators of the saved states
(`A4/output.nltex_jmc`, `A7/output.nltex_jmc`, 30 MB compressed, given
separately because of their size): with them the first iteration starts
from the radiation field of the saved state; without them it builds its
macro-atom tables with J = B(T_R), and the later iterations recover. The
estimators belong to the lines of their atoms: those of the sn3 states are
not used with sn6 (SuperLite says so and starts from J = B(T_R)).

Time (4 ranks on a recent desktop): with sn6 the NLTE solution takes most
of an iteration, 2–3 times that of sn3; about ten iterations from the
saved states (`NITER=10`; the added ions start far from their solution).
To continue a run, `INIT=<its work dir>`.

## Checks

- `grep "NLTE (it" work/A4/run.log`: the median change per iteration
  should fall to 0.1–0.2 (the floor of the Monte Carlo noise; the spectra
  then move by ~5 %); if not, run more (`NITER`, `INIT`).
- The spectra against observed IIn spectra near day 20 (rest frame;
  e.g. SN 2010jl at +19 d, SN 2006gy at +23 d), the models smoothed to the
  observed resolution:
  ```
  Tools/compare_obs_iin.py work/day20_vs_obs.png --res 600 \
    --model "A4=work/A4_fi/output.fi_spectrum" \
    --model "A7=work/A7_fi/output.fi_spectrum" \
    --obs "SN 2010jl +19 d=sn2010jl.p0019.dat" --obs "SN 2006gy +23 d=sn2006gy.p0023.dat"
  ```
  (prints the Hα equivalent widths; earlier runs: §4.3 of `CSM_PHYSICS.md`).
