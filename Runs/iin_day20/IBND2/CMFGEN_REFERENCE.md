# IBND2: a reference calculation with CMFGEN

The same model computed with CMFGEN tests SuperLite's NLTE solution and
spectrum against an established code, apart from whether this
parameterized model fits SN 2006jc.  Everything a second code needs is
below; `cmfgen_ref_table.txt` holds the structure cell by cell.

## Structure (`cmfgen_ref_table.txt`)

Written by `Tools/export_model_table.py Runs/iin_day20/IBND2 ... output.nltex`
(SuperLite's IBND4 state: 16 iterations, same-element overlaps off).

| Region | Cells | r [cm] | v [km/s] | Mass [Msun] | Composition |
|---|---|---|---|---|---|
| Ejecta (rho ~ r^-8) | 1-16 | 1.20e15 - 1.93e15 | 1200 - 1930 (homologous) | 0.20 | `he4mix` |
| Dense shell (CDS) | 17-20 | 1.93e15 - 2.00e15 | 2000 | 1.50 | `he4mix` |
| Smeared shock | 21-28 | 2.00e15 - 2.10e15 | 2000 -> 600 | | CSM (He 0.987) |
| Unshocked wind (Mdot 0.05 Msun/yr, 600 km/s) | 29-78 | 2.10e15 - 6.0e15 | 600 | | CSM (He 0.987) |

`he4mix` (mass fractions): He 0.50, C 0.04, N 0.005, O 0.216, Ne 0.06,
Na 0.002, Mg 0.04, Si 0.07, S 0.03, Ar 0.005, Ca 0.012, Fe 0.016,
Ni 0.004, after the he4p0 model of Dessart, Hillier & Kuncarayakti 2022
(A&A 658, A130).  No H.  Time since explosion: 24 d (only through the
velocity field: the flow is not homologous outside the ejecta).

The velocity field is not monotonic in r/v: constant across the shell,
falling across the shock, constant in the wind.  If CMFGEN needs a
monotonic or homologous v(r), the closest equivalent is to keep the
densities and give the shell and everything outside it the shell's
velocity gradient of the run (or 2000 km/s throughout with a small
gradient); say which, and SuperLite will be run on the same field.

## Energy input

- **6.0e42 erg/s deposited in the shell**, 1.5e42 erg/s in each of
  cells 17-20 (about 1.8e-3 erg/s/cm^3), all as fast electrons:
  SuperLite solves their Spencer-Fano degradation with all NLTE species
  (`in_nt_sf = t`, `in_nt_lumg = 1`), so the power goes to heating,
  ionization and excitation in the shares of that solution.  In CMFGEN
  this is non-thermal deposition with a prescribed local energy
  deposition rate (as radioactive decay, but with this profile instead
  of 56Ni), which is how Dessart et al. 2022 powered their CDS.
- **6.0e41 erg/s at the inner boundary** (1.2e15 cm): a Planck spectrum
  at 15 218 K, isotropic, the core reflecting (SuperLite `refl`).  A
  diffusion or hollow-core boundary with the same luminosity is the
  CMFGEN analogue.
- No shock emission: the smeared-shock cells add no power.
- Total L = 6.6e42 erg/s.  Steady state (no time-dependent terms).

## Physics settings of the SuperLite run to match

- NLTE for H, He, C, N, O, Ne, Na, Mg, Si, S, Ca, Fe with every ion
  up to the bare nucleus (CMFGEN ions from atomic_data_21jun23 with
  `GF_CUT = 1e-3`, the rest from CHIANTI); Ar and Ni in LTE.
- Lines in the Sobolev approximation with Doppler profiles of 10 km/s
  microturbulence (`in_nltex_vturb`), level dissolution on
  (`in_nltex_dis`), charge exchange on.
- Overlapping lines: pairs with an H or He line only (no same-element
  pairs of metals: they are off in this run).
- No clumping (IBNDC has the shell at f = 0.1).

## What to compare

Per cell: T_e, n_e, the ion fractions of He, O, Mg, Si, Ca, Fe (the
SuperLite values are the last columns of the table and in
`output.nltex`).  Emergent spectrum 3500-10000 A: He I 4471, 5016,
5876, 6678, 7065, the He II 4686 residual, Fe II below 5500 A, Ca II
H&K and 8500 A, Mg II 4481; the blue-to-red flux ratio
F(3700-5000)/F(6000-7300) (SuperLite: 2.5).  SuperLite's spectrum is a
formal integral with 12 electron-scattering passes, 10 km/s sampling.
