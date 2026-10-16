# SuperLite CSM-physics branch (`csm-physics`)

Goal: spectra of interacting supernovae — from no CSM to dense winds and shells of
any mass — that can be compared with observations, with physically correct NLTE line
formation, non-monotonic velocity fields and radiation transport through shocks.

This file records (1) what the released code gets wrong for that goal, (2) what
stages 1 and 2 of this branch change and how they are validated, (3) the roadmap
and (4) results.

## 1. Assessment of the released code (HEAD 176d0a9)

| # | Where | What | Consequence for interacting SNe |
|---|---|---|---|
| 1 | `SOURCE/sourceenergy.f90`, `SOURCE/interior_source.f90` | Every cell emits a blackbody *surface* at its radiation temperature: energy ∝ r² T_rad⁴, normalized to `in_L_bol`, Planck spectrum, directions outward only (μ = √r). | Optically thin CSM cells emit continuum they cannot emit; the energy goes where r is large, not where it is generated (shock, inner boundary). Line emission only appears after absorption. |
| 2 | `GAS/physical_opacity_nlte.f` | The "NLTE" H rate equations use J_ν = B_ν(T_gas) for every line, so bound-bound rates are in detailed balance at T_gas (excitation ≈ LTE). The photoionization fit gives ~10⁻²² s⁻¹ (negligible); the recombination term is evaluated at T_rad and lacks n_e. Ionization comes from LTE Saha (`GAS/eos_update.f`). | Balmer source functions near LTE at the local gas temperature → P-Cygni/absorption profiles, flat decrements; no photoionization–recombination physics. |
| 3 | `GAS/physical_opacity*.f` | Line opacity is the line strength smeared over the group width (no Sobolev saturation) and treated as true absorption with thermal re-emission over the whole spectrum (ε = 1). | An optically thick Balmer line is a thermalization sink: absorbed Hα energy is re-emitted mostly at other wavelengths. |
| 4 | `TRANSPORT1/transport11.f90` | At a Doppler (group-crossing) event the lab wavelength is reset with the Doppler factor *before* the move (`elabfactold`), so it jumps redward at every crossing; free-streaming photons drift redward by about the comoving velocity range they cross. The crossing distance uses ∇·v/3 (exact only for v ∝ r); blueshifts are decided from dv/dr alone. | Line profiles shifted and distorted (A4 day 20: Hα absorption trough at −1000 instead of −4000 km/s); non-monotonic flows (ejecta → shock → slow wind) treated approximately. |
| 5 | `transport11.f90` | Thomson scattering is coherent in the comoving frame. | No electron-scattering wings (the hallmark of SN IIn line profiles). |
| 6 | whole code | No shock physics: the gas temperature comes from the hydro code (single temperature, shocks smeared over zones, ~10⁴ K at day 20 in A4), no EUV/X-ray emission from shocks. | The unshocked CSM cannot be photoionized → no narrow recombination lines. |
| 7 | `mpimod_mpi.f` `reduce_gridtally` | The energy tally is reduced to rank 0 only, but every rank updates T_rad from it (`nlte_rad_est`). | With `in_niter > 1` and several MPI ranks, T_rad and the source distribution differ between ranks. |
| 8 | build/inputs | LAPACK never linked (CMake variable typo); gfortran ≥ 10 needs `-fallow-argument-mismatch`; `read_fluxgrid` never sets the number of bins (file flux grids never worked); Testsuite A4 10/20/30 d: `in_ndim` ≠ cells in `input.str`, unknown `ni56` column, day-5 `in_L_bol`. | — |
| 9 | atomic data | H lines at vacuum wavelengths (recomputed from level energies), metal lines at air wavelengths. | 83 km/s offset of the H lines relative to metals and to observed (air) spectra. |
| 10 | design | Time-independent snapshot; no non-thermal (γ-ray, X-ray) excitation/ionization; NLTE only for H; no light-travel time. | SNe II lose Hα too early without time-dependent ionization and non-thermal excitation; He lines LTE. |
| 11 | `GAS/physical_opacity*.f` | Free-free opacity counts every atom as a fully ionized ion of nuclear charge Z (`iz**2*gas_natomfr`, times n_atom n_e). | Neutral gas (e.g. CSM beyond the ionization front) gets free-free opacity and emission as if fully ionized; metals are weighted by Z² (Fe: 676). |
| 12 | `nltemod.f` | `nlte_dealloc` loops over NLTE data that is never allocated when the bb opacity is off. | Segmentation fault at exit of continuum-only runs (after the output is written). |
| 13 | `ionsmod.f` `ions_solve_eos` | Partition functions are computed only while they are all zero, i.e. once, at the temperature of the first cell of each rank, and reused for every cell and iteration. | LTE ionization balance and level populations of all elements use partition functions at the wrong temperature, different on each MPI rank. |
| 14 | `ionsmod.f` `ions_nebular` (`in_nebular`, experimental) | ζ and the ionization energy in δ are taken from the next ion up (data.zeta is indexed by the charge of the recombining ion: ζ = 0 for H); the Saha factor is evaluated at the gas temperature instead of T_R; W is geometric from one photospheric radius. | Not used by this branch; `in_tion = 'nebu'` is a corrected implementation with Monte Carlo W and T_R. |

## 2. Stage 1 (this branch)

All new physics is switchable; with the switches off the code reproduces the released
behaviour except for the bug fixes in rows 4 (exact Doppler, switch back with
`in_trn_doplegacy = t`), 7 and 8.

| Switch | Physics | Code |
|---|---|---|
| (default) `in_trn_doplegacy = f` | Exact distance to the next comoving group boundary (or H line) for the cell's linear v(r) of any slope; up to two stationary points of μv along a path, so blueshifts and repeated crossings in decelerating and non-monotonic flows are exact. Lab wavelength conserved in free flight. | `TRANSPORT1/dopplermod.f90` |
| `in_thm_thermal = t` | Thomson scattering off Maxwellian electrons: electron-rest-frame boost, isotropic or dipole, Compton recoil. Broadening √(2(1−cosθ)kT/m_ec²) per scattering, Kompaneets mean shift. | `TRANSPORT1/comptonmod.f90` |
| `in_nlteh = t` (+ `in_nlte = t`) | NLTE hydrogen, n-resolved atom (`in_nlteh_nlev` = 30): Johnson (1972) f-values and Gaunt factors, Milne-relation recombination, collision rates from `Atoms/data.nlte.h1`, Seaton ionization, three-body recombination, 2s–1s two-photon decay. Photoionization by the Monte Carlo J_ν; line rates in the Sobolev approximation with the **direction-dependent** gradient Q(μ) = μ²dv/dr + (1−μ²)v/r and the resonance estimator of the incident intensity. | `nltehmod.f90`, `GAS/nlteh_gas_update.f90` |
| (with `in_nlteh`) | H lines are Sobolev resonance events in the transport (found with the exact crossing, so a line can be met several times in a non-monotonic flow), interaction probability 1−e^{−τ(μ)}. Lucy's macro-atom with k-packets for the energy flow through H (recombination cascades, fluorescence, collisional excitation of hot gas), solved per cell as an absorbing Markov chain. H bound-free opacity from all levels (Balmer, Paschen continua). | `nltehtrmod.f90`, `transport11.f90` |
| `in_srctype = 'lumg'` | Energy injected where the hydro L(r) grows outward (`lum` column), plus L at the inner boundary; emitted isotropically with the local emissivity (k-packets). | `SOURCE/sourceenergy.f90`, `interior_source.f90` |
| `in_shock_dvmin`, `in_shock_tmax` | Shock cells (velocity drop > 50 km/s and 2 % per cell, ≥ `in_shock_dvmin` per shock) emit their injected energy as a radiatively cooling flow: T uniform between T_gas and min(T_postshock, `in_shock_tmax`), bremsstrahlung shape. Its EUV/X-rays photoionize the CSM and the dense shell. | `SOURCE/shock_temperatures.f90` |
| tools | `Tools/stella2superlite_csm.py`: STELLA profile → `input.str` with `lum`, whole CSM, shock-preserving zone merging, recommended `input.par`. `Tools/make_wlgrid.py`: 100 km/s groups and flux bins in the optical, coarse elsewhere. | `Tools/` |

### Validation

Unit tests (`Tests/unit/`, build with gfortran + LAPACK):

| Test | Checks | Result |
|---|---|---|
| `test_doppler` | 20 000 random cells (accelerating, decelerating, non-monotonic): crossing distance and direction against a brute-force scan | distance to 5×10⁻⁶ of the path (scan resolution), boundary hit to 10⁻¹¹ |
| `test_compton` | ⟨Δν/ν⟩ = (4kT − hν)/m_ec², ⟨(Δν/ν)²⟩ = 2kT/m_ec², exact recoil averages (1 Å), isotropy | 0.3 % (variance), within 1σ (mean) |
| `test_nlteh` | hydrogenic data; α_B against Hui & Gnedin (1997); J = B(T) ⇒ Saha–Boltzmann at 10⁶–10¹⁴ cm⁻³; Case B decrement | LTE to < 0.3 %; Case B Hα/Hβ = 2.74 (l-bundled atom at 10⁴ cm⁻³; Storey & Hummer 2.86 for l-resolved) |
| `test_macroatom` | macro-atom Hα/Hβ energy per recombination against the statistical-equilibrium emissivities | 2.720 ± 0.020 vs 2.737 |

Integration tests and the comparison with observed SNe IIn: see §4.

### Stage 2

| Switch | Physics | Code |
|---|---|---|
| (with `in_nlteh`) | A line photon leaves its Sobolev resonance in direction μ with probability ∝ β(μ) = (1−e^{−τ(μ)})/τ(μ), after resonance absorption and after macro-atom emission. Isotropic in homologous flows; in a constant-velocity wind photons cannot escape radially (τ → ∞). | `transport11.f90` (`line_mu`) |
| `in_bbsobolev = t` | Lines of the LTE elements: Sobolev expansion opacity with the local gradient, averaged over directions, (1/(λΔν))⟨\|Q\|(1−e^{−K/\|Q\|})⟩, instead of the line strength smeared over the group (which overestimates a thick line by its τ, up to 10⁶ in a slow CSM). | `bbsobmod.f90`, `GAS/physical_opacity*.f` |
| `in_bbeps = -1` | A line interaction thermalizes the photon with probability ε = X/(X+β), X = C_ul/A_ul = 20.6 λ³ n_e P(ΔE/kT)/√T (van Regemorter 1962), β the Sobolev escape probability; otherwise coherent scattering in the comoving frame (new IMC event; part of the extinction in the DDMC thresholds and leakage). `in_bbeps` ∈ [0,1] sets a constant ε; 1 is the released code. | `bbsobmod.f90`, `transport11.f90`, `diffusion11.f90`, `leakage_opacity11.f90` |
| `in_tesolve = t` | Electron temperature from the thermal balance with the radiation field of the previous iteration: LTE opacities 4π Σ κ(B(T)−J); H photoionization heating, recombination and collisional cooling with the H populations re-solved at every trial T; Compton heating/cooling with induced scattering; energy injected outside shocks. Damped (`in_te_damp`) between iterations. | `tebalmod.f90`, `GAS/te_update.f90` |
| `in_tion` | Ionization/excitation of the LTE elements: Saha–Boltzmann at the hydro temperature (`'hydr'`, default) or at T_e (`'elec'`), or the modified nebular approximation (Mazzali & Lucy 1993) with the radiation temperature (mean frequency) and dilution W of the Monte Carlo field of each cell and T_e (`'nebu'`; reduces to Saha for W = 1, T_e = T_R). Saha at T_e over-ionizes metals in a dilute field; at the hydro temperature they stay too neutral where the gas is photoionized and hot. | `GAS/eos_update.f`, `GAS/te_update.f90` (`nebular_wtr`), `ionsmod.f` (`ions_nebular_mc`) |
| `in_nt_emin = 100` | Photons above 100 eV make fast photoelectrons whose energy heats the gas, ionizes and excites H in the shares of Shull & van Steenberg (1985): in the transport (absorbed X-ray packets go to the k-pool, the H continuum or H n = 2), in the H rate equations (non-thermal ionization and 1→2 rates from the absorbed X-ray power) and in the thermal balance (only the heat share heats). | `nltehtrmod.f90`, `nltehmod.f90`, `tebalmod.f90` |
| (default) `in_fflegacy = f` | Free-free opacity from the ion charges of the EOS (NLTE H II fraction for H) instead of nuclear charges for all atoms (row 11). | `GAS/physical_opacity*.f`, `GAS/eos_update.f` |

Fixes: the H bound-free opacity of stage 1 stopped at the Lyman edge (the level loop
exited at n = 1 instead of skipping levels with thresholds above the frequency), so
the Balmer, Paschen, ... continua were missing from the transport (the rate equations
were right); OpenMP clauses; rows 11, 12 and 13.

Validation:

| Test | Checks | Result |
|---|---|---|
| `Tests/benchmarks/wind` | pure-scattering Sobolev line in constant-velocity, β-law and homologous winds against the formal solution with direction-dependent escape (`sobolev.py`) | χ²/dof = 1.17, 1.34 (v = const, τ₀ = 3, 30), 0.98, 0.99 (β law), 1.20 (homologous); isotropic re-emission (stage 1): 5.15 (8σ in the trough) |
| `test_bbsob` | angle-averaged expansion opacity against a brute-force angle integral for random accelerating, decelerating and compressing gradients; thin and thick limits; homologous (1−e^{−Kt})/t; van Regemorter table | 0.26 % (worst case), exact limits |
| `test_tebal` | Planck field J = B(T₀) with LTE opacities and NLTE H must return T₀; pure-H nebula around a 4×10⁴ K star | T₀ to 0.2 % at n_H = 10⁸ and 10¹³ cm⁻³ from both sides; 21 500 K |
| `test_nebular` | modified nebular approximation: W = 1, T_e = T_R reproduces the LTE Saha solution (all elements, 6000–20 000 K, 10⁹ and 10¹³ cm⁻³); ionization falls with dilution | fractions to 10⁻⁸; monotonic |
| all unit tests | `SL_BUILD=<build> Tests/unit/run_tests.sh` | pass |
| regression | stage-2 switches off and `in_fflegacy = t`: LTE runs bit-identical to stage 1 before the partition-function fix (row 13), which changes every LTE run | identical |

## 3. Roadmap

Done in stage 2: electron temperature from the thermal balance, non-thermal
ionization/excitation of H by shock X-rays, Sobolev opacity with thermalization for the
lines of the other elements, direction-dependent Sobolev escape (§2).

**Next — physics needed for reliable line strengths across regimes**
1. He I / He II macro-atoms with the same machinery (He I 5876/6678/7065, He II 4686: SN IIn/Ibn/Icn, flash features; the C-series). He II is hydrogenic: `nltehmod` with Z = 2 scaling is the quick route.
2. Ionization of the other elements beyond LTE: photoionization equilibrium with the Monte Carlo J (or the nebular approximation with Monte Carlo W and T_R) instead of Saha at one temperature (`in_tion`).
3. Time-dependent ionization (Utrobin & Chugai 2005; Dessart & Hillier 2008): needed for SN II Hα on the plateau (no-CSM models lose Hα too early).
4. Non-thermal deposition of ⁵⁶Co γ-rays (Kozma & Fransson 1992) with the same Shull & van Steenberg channels.
5. Shock structure: adiabatic vs radiative from the cooling time; cooling-function spectra (line cooling at 10⁵–10⁷ K) instead of bremsstrahlung shapes; T_e ≠ T_i behind fast shocks.
6. NLTE macro-atoms for Fe II, Ca II, O I, Mg II (the ε treatment of stage 2 is a two-level-atom approximation).
7. DDMC: macro-atom-aware absorption, non-thermal deposition and H lines (expansion opacity) in optically thick cells (at present DDMC re-emits absorbed energy thermally).

**Stage 3 — dynamics and geometry**
8. Time-dependent radiation transport (SuperNu's IMC with census) with light-travel time, for interaction phases where t_diff ≈ t_dyn.
9. 2D/3D CSM (disks, clumps, asymmetric shells) with the SuperNu 2D/3D transport.

**Stage 4 — validation and production**
10. The observed-line-metric harness (EW, decrement, absorption share, peak velocity on a common ±5000 km/s window; `balmer_validation`) as regression tests; grids for the A/B/C model series; comparison with CMFGEN/ARTIS spectra of the same hydro snapshots.

## 4. Results

### 4.1 Transport benchmark

`Tests/benchmarks/pcygni/`: pure-scattering line (τ = 3 at the photosphere,
τ ∝ r⁻⁷, homologous), photosphere as a Lambertian inner boundary
(`in_innerbc = 'thrm'`), against the analytic Sobolev profile (`analytic.py`):
χ²/dof = 1.21 over 32 velocity bins, largest deviation 2.9σ.

### 4.2 Stage-1 integration runs (STELLA A4 and A7 at day 20)

Line metrics on the common ±5000 km/s window (EW positive = net emission), measured
with the same code as the observations (`balmer_validation`):

| Run | L_out (erg/s) | EW5(Hα) Å | EW5(Hβ) Å | Hα absorption share | Hα peak (km/s) |
|---|---|---|---|---|---|
| legacy (A4) | 1.1e42 | +1.5 | +2.8 | 0.46 | 50 |
| exact Doppler only | 1.1e42 | −15.6 | −9.3 | 0.78 | 300 |
| + NLTE H, legacy surface source | 1.1e42 | −95 | −79 | 1.00 | −5000 |
| stage 1, all switches (A4, 0.2 M☉ CSM) | 2.2e42 | −79 | −49 | 1.00 | 300 |
| stage 1, all switches (A7, 2 M☉ CSM) | 1.0e43 | +17 | −55 | 0.31 | 3950 |
| observed SNe IIn, 10–30 d (median, n = 21/20) | — | ≈ +100 | ≈ +25 | 0 | 0 |

What the runs show:
* The exact Doppler treatment moves the A4 Hα absorption trough from −1000 to
  −3000…−5000 km/s, where the line-forming ejecta are.
* With NLTE H the CSM of A7 is photoionized by the shock radiation out to
  4.35×10¹⁵ cm (n_e ≈ 10⁹–10¹⁰ cm⁻³, τ_es ≈ 5); b₃/b₂ grows outward to ≈ 3, so
  the Hα source function exceeds B(T_gas). The cold dense shell stays neutral.
* Hα is still far too weak. A slow (300 km/s) optically thick CSM can emit an
  optically thick Sobolev line only at a rate ∝ S·|dv/dr|; with the hydro gas
  temperature (12 000 → 4700 K) this gives EW ≈ 10–20 Å. The observed ≈ 100 Å
  needs the CSM electron temperature from heating by the shock radiation instead of
  the hydro temperature (B_Hα rises 27× from 6000 to 25 000 K).
* The shock luminosity absorbed by metals is re-emitted as LTE metal lines (Ca II,
  Si II 6355, Fe II), broadened by electron scattering. Two causes: metal line
  opacity is smeared over the group without Sobolev saturation, which overestimates
  it by up to τ_line in the slow CSM, and lines are pure absorbers (ε = 1).

These define stage 2 below.

## 5. Running

```
cmake -DUSE_MPI=ON <src> && make
python3 Tools/stella2superlite_csm.py mesa.day020_post_Lbol_max.data -o run --tau-max 50
python3 Tools/make_wlgrid.py -o run
# input.par: in_grp_ng = 0, in_flx_ndim = -1,1,1 (written grids)
mpirun -np 4 superlite
```
Outputs: `output.flx_luminos` (spectrum per iteration), `output.nlteh` (populations,
departure coefficients, Sobolev optical depths, escape probabilities per cell),
`output.tebal_<rank>` with `in_tesolve` (heating and cooling terms per cell at the last
iteration).  Unit tests: `Tests/unit/run_tests.sh`.
