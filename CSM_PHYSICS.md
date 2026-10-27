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
| `in_srctype = 'lumg'` | Energy injected where the hydro L(r) grows outward (`lum` column), plus L at the inner boundary; emitted isotropically with the local emissivity (k-packets). The inner-boundary packets have a Planck spectrum at the input radiation temperature of the innermost cell (or `in_srcmax`), held fixed: following the Monte Carlo radiation temperature of that cell fed back on itself (harder injection → more trapping near the boundary → higher T_rad), +18 % in four iterations in the CMFGEN benchmark. | `SOURCE/sourceenergy.f90`, `interior_source.f90`, `boundary_source.f90` |
| `in_innerbc` (`refl`, `thrm`, `absb`) | The inner boundary of a grid with a hollow core (r_in > 0) is a surface of the innermost cell: a packet that reaches it is reflected with a Lambertian direction (`refl`, default), re-emitted with a Lambertian direction and the Planck spectrum at `in_srcmax` or the input T_rad of the innermost cell (`thrm`, comoving energy conserved), or absorbed (`absb`). Before (released code), the inward rays of the innermost cell crossed the core as part of the cell, with its opacity and their path in its J, and the reflection had a uniform μ: the J of the thin innermost cell of the CMFGEN benchmark (1.3×10¹¹ cm) was 2.1–2.8 times CMFGEN's; it now follows its neighbour, and no other cell or the spectrum changes by more than 2 %. Homologous flows hardly saw it: an inward photon met a group boundary inside the core. | `TRANSPORT1/transport11.f90` |
| `in_srcdiff` | The luminosity of the inner boundary (`lumg`) injected with the flux spectrum of the diffusion approximation, H_ν ∝ (∂B_ν/∂T)/χ_ν with the extinction of the innermost cell (Rosseland weighting; T = `in_srcmax` or the input T_rad of the innermost cell), as CMFGEN's inner boundary, instead of the Planck spectrum. For a grey opacity the Planck spectrum has (x/4)eˣ/(eˣ − 1) times the diffusion flux (x = hν/kT): at 17 000 K 2.1 times less at 1000 Å and 2.1 times more at 6000 Å. In the H/He benchmark it raises J at 925–1020 Å by at most 4 % in any cell (one iteration on the same populations, `bcG_d` against `bcG_ref`), a thermalizing boundary (`in_innerbc = 'thrm'`, `bcG_t`) by 2–10 %: the innermost cells reprocess the injected spectrum, and neither accounts for SuperLite's far-UV J being 0.7–0.9 times CMFGEN's (§4.4). | `SOURCE/boundary_source.f90` |
| `in_shock_dvmin`, `in_shock_tmax` | Shock cells (velocity drop > 50 km/s and 2 % per cell, ≥ `in_shock_dvmin` per shock) emit their injected energy as a radiatively cooling flow: T uniform between T_gas and min(T_postshock, `in_shock_tmax`) (isobaric cooling releases the same energy per unit T), bremsstrahlung shape. Its EUV/X-rays photoionize the CSM and the dense shell. | `SOURCE/shock_temperatures.f90` |
| `in_shock_cie` | At each T of the cooling flow the photons follow the spectrum of gas in collisional ionization equilibrium of the cell's composition (lines, free-free, free-bound, two-photon; 10^4.5–10^8 K) instead of bremsstrahlung: at 10⁵–10⁷ K line cooling exceeds bremsstrahlung 10–100 times and comes out in the EUV/soft X-rays and the UV/optical coronal lines. Per-element tables from CHIANTI 11 (Dere et al. 1997; Del Zanna et al. 2021) made with `Tools/chianti_cie_spectra.py` (ChiantiPy; the table is not distributed). The shock cells' emission per group (`shc_group_frac`: the cooling flow over T, CIE or bremsstrahlung, as the source draws it) is part of the emissivity of their source functions, so the photoionizations in their thick groups (local detailed balance) take the shock's EUV photons, as with the X-ray source (`in_xray_*`). | `shockciemod.f90`, `SOURCE/interior_source.f90`, `superlite.f90` (`shock_emnbf`) |
| tools | `Tools/stella2superlite_csm.py`: STELLA profile → `input.str` with `lum`, whole CSM, shock-preserving zone merging, recommended `input.par`. `Tools/make_wlgrid.py`: 100 km/s groups and flux bins in the optical, coarse elsewhere, with group boundaries at the main bound-free edges (H, He, abundant metals; `--from` adds them to an existing grid): in a group across an edge the cross-section is the group average, so the photons of the lines that converge on the edge from the red photoionize. In the H/He benchmark the n = 4 edge (14 588 Å) inside a 1000 km/s group lowered b₅ by up to 40 % below CMFGEN's in the inner wind (Pβ forms there); with the edge as a boundary b₅ agrees within 5 % (§4.4). Use edge-aligned grids for line work. | `Tools/` |

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
were right); OpenMP clauses; rows 11, 12 and 13.  NaN in the absorption or emission
opacity now stops the run (the released code printed a warning per cell and group and
went on).

Validation:

| Test | Checks | Result |
|---|---|---|
| `Tests/benchmarks/wind` | pure-scattering Sobolev line in constant-velocity, β-law and homologous winds against the formal solution with direction-dependent escape (`sobolev.py`); non-monotonic flows (a IIn-like shell: homologous ejecta to 2000 km/s, a fall to a 200 km/s CSM; an accelerating then decelerating wind) against the generalized Sobolev solution of Rybicki & Hummer (1978) with the coupling of all the resonances of a ray (`sobolev_nl.py`); the resonance estimator of the incident intensity per cell against ⟨I_inc⟩ and, β-weighted (`in_nltex_jbbeta`), against S (`compare_s.py`) | χ²/dof = 1.17, 1.34 (v = const, τ₀ = 3, 30), 0.98, 0.99 (β law), 1.20 (homologous); isotropic re-emission (stage 1): 5.15 (8σ in the trough). Shell 1.14 (40 km/s bins at the CSM's narrow absorption: 0.93), accelerating–decelerating 1.20, where the profile of the local Sobolev S (no coupling, off by up to 1300 %) has 27.8 (19.6σ); Rybicki & Hummer's decelerating wind v ∝ r^(−1/2) out to 100 R_c (τ₀ = 100) 1.65 against 15.2 for their disconnected approximation, line-centre flux 1.34 times the disconnected one (theirs: 1.35). The reference itself reproduces their Table 2 (S within 5 %, `rh78_table2.py`). Incident intensity: β-weighted against S 0.83, 1.47, 1.09, 1.00 (constant velocity, shell, accelerating–decelerating, Rybicki & Hummer; cells with ≥ 100 events), where ⟨I_inc⟩/S reaches 129, 30, 10 and 2.6 |
| `test_bbsob` | angle-averaged expansion opacity against a brute-force angle integral for random accelerating, decelerating and compressing gradients; thin and thick limits; homologous (1−e^{−Kt})/t; van Regemorter table | 0.26 % (worst case), exact limits |
| `test_tebal` | Planck field J = B(T₀) with LTE opacities and NLTE H must return T₀; pure-H nebula around a 4×10⁴ K star | T₀ to 0.2 % at n_H = 10⁸ and 10¹³ cm⁻³ from both sides; 21 500 K |
| `test_nebular` | modified nebular approximation: W = 1, T_e = T_R reproduces the LTE Saha solution (all elements, 6000–20 000 K, 10⁹ and 10¹³ cm⁻³); ionization falls with dilution | fractions to 10⁻⁸; monotonic |
| all unit tests | `SL_BUILD=<build> Tests/unit/run_tests.sh` | pass |
| regression | stage-2 switches off and `in_fflegacy = t`: LTE runs bit-identical to stage 1 before the partition-function fix (row 13), which changes every LTE run | identical |

### Stage 3: NLTE for all abundant elements

Stage 3 replaces the H-only machinery by one framework for any set of ions (the
built-in H atom of stage 1 is one of them; `nltehtrmod.f90` is gone, its transport and
macro-atom are `nltextrmod.f90`; `GAS/nlteh_gas_update.f90` is
`GAS/nltex_gas_update.f90`; `test_macroatom` is part of `test_nltex`).

| Switch | Physics | Code |
|---|---|---|
| `in_nltex_atoms`, `in_nltex_dir` | Model atoms converted from the CMFGEN atomic data (Hillier & Miller 1998; `Tools/cmfgen2superlite_atom.py`, `Tools/make_nlte_atoms.sh`; the data carry no license statement, so the converted files are not distributed): levels, f-values, collision strengths Ω(T) (missing pairs as in CMFGEN: van Regemorter with the Auer & Mihalas or ion Gaunt factors, Ω_set; Seaton ionization), photoionization cross-sections of every CMFGEN type; further routes to excited terms of the next ion are added to the route to the ground state (exact for a final term in Boltzmann equilibrium with the ground state). Default set: He I–II, C I–IV, N I–V, O I–VI (the high ions of flash-ionized CSM: N V 4604/4620, O V, O VI 3811/3834), Ne I–III, Na I–II, Mg I–III, Si I–IV, S I–IV, Ca I–III, Fe I–IV, Co II–IV, Ni II–IV (≈ 2000 states, 14 000 lines). | `nlteatommod.f90` |
| (with the atoms) | Statistical equilibrium of each element over all its ions and the next ion (top state): Sobolev bound-bound rates with the Q(μ) escape probabilities and resonance estimators of the incident J (a line with fewer than 30 crossings in a cell takes the J of its group without the photons of the group's lines that have kept their frequency history since their emission: in an expanding flow those only redshift away from their line, so the incident intensity is the continuum's, as in CMFGEN's Sobolev option; §4.4); photoionization and stimulated recombination by the Monte Carlo J of the groups; spontaneous recombination from the Milne relation, integrated accurately (the stage-1 H tables overestimated α of high n at high T by up to 70 %); collisional excitation, ionization and three-body recombination; two-photon decays; non-thermal ionization and excitation of H and He (Shull & van Steenberg 1985, He shares scaled with n_He/n_H). Populations, escape probabilities and n_e iterated together per cell. | `nltexmod.f90` |
| (with the atoms) | Every model-atom line is a Sobolev resonance in the transport; one macro-atom per element (Lucy 2002, 2003) with a shared k-packet pool, its internal transitions solved per cell as an absorbing Markov chain; bound-free absorption from every level with stimulated recombination. The absorbing level is chosen exactly with n_l σ_l(ν) and without a walk over all levels (13 000 in the ζ Pup set): per cell and each of 256 blocks of groups, the 256 levels of largest n_l max σ_l form a list, and a uniform number on [0, s + R) picks from the list below its weight s at the group. Above s, the other levels (bound R) are walked, and a number beyond their weight is drawn again. The walk happens with probability R/(s + R), 8×10⁻⁴ on average for the 44 ions at 25 000 K. | `nltextrmod.f90`, `TRANSPORT1/transport11.f90` |
| super-levels (the `--ftos` column of the model atoms) | CMFGEN's super-levels (F_TO_S assignments of the atomic data): the levels of a super-level share its departure coefficient, Boltzmann at T_e within it. The statistical equilibrium, the macro-atom (state energies: the Boltzmann mean of the levels; emitted energies those of the level pairs) and the thermal balance run over the super-levels; lines, bound-free opacity and emission channels keep the full levels. Collision rates between super-levels are summed without a loop over all level pairs (Ω_set pairs by running counts, own-Ω pairs as corrections). Lines within a super-level are thermal (absorbed into and emitted by the k-packet, cooling n_u A β hν). The atoms can then reach their ionization limits: Fe II 827 levels in 209 super-levels instead of 100 levels. | `nlteatommod.f90`, `nltexmod.f90`, `nltextrmod.f90`, `tebalmod.f90` |
| `in_nltex_chgexch` (default `chgexch.dat` in `in_nltex_dir`) | Charge exchange among the NLTE elements (CMFGEN's reactions: C, N, O, Fe, Co, Ni ions with H; C, N ions with He; Kingdon & Ferland 1996 fits; `Tools/cmfgen2superlite_chgexch.py`): forward rates into the product term with Boltzmann level shares, reverse rates by detailed balance with the model's own levels and ionization energies, the released energy heats the gas. In partially ionized H it sets O I/O II and Fe II/Fe III (Fe III + H⁰: ≈ 10⁻⁹ cm³ s⁻¹, orders of magnitude above radiative recombination). | `nlteatommod.f90`, `nltexmod.f90` |
| `DIE` of `Tools/make_nlte_atoms.sh` (`--die` of the converter) | Low-temperature dielectronic recombination as in CMFGEN (`[DIE_XzV]`, `RD_PHOT_DIE_V1`). Each transition of CMFGEN's `die_data` (C III, N IV, O V, Co II, Sc I–II; `DIE_FILE=htdr_data` for the larger C III and S V sets) runs from an autoionizing level, in LTE with the next ion, to a model level. It becomes a resonance of strength πe²f/mc in that level's cross-section, at the frequency of the stabilizing transition (Gaussian, 1000 km/s, `VSM_DIE`). Three things follow: the Milne recombination through it is the LTE population of the autoionizing level times A; photoionization through it counts as bound-free opacity; and its recombination photons fall at the stabilizing line. The default `wi` adds only the transitions from levels that cannot autoionize in LS coupling, since the others are already resonances of the Opacity-Project cross-sections. At 10⁴ K, C IV→C III recombination rises from 7.8×10⁻¹² to 1.6×10⁻¹¹ cm³ s⁻¹ and N V→N IV from 6.7×10⁻¹² to 1.1×10⁻¹¹. | `Tools/cmfgen2superlite_atom.py` |
| (with the atoms) | Lines of the line list of NLTE ions whose upper level lies above the model atom, and lines of the top ion, stay in the Sobolev expansion opacity, the lower level populated from the NLTE population of the highest model level at or below it with Boltzmann excitation at T_e: the truncated Fe-group atoms (35–50 % of the ionization energy) keep their UV blanketing. `in_nltex_bbext = f` leaves them out: only the lines of the model atoms, as a CMFGEN model with the same atoms (for benchmarks: with them SuperLite's J at 912–4000 Å in the full-metal wind was 1.2–1.7 times CMFGEN's on CMFGEN's own populations, without them 1.0–1.2, §4.4). | `GAS/physical_opacity_nlte.f` |
| `in_nltex_dis` | Level dissolution by the plasma microfield (Hummer & Mihalas 1988), as CMFGEN's `DO_DIS`, with the occupation-probability fit of Hubeny, Hummer & Lanz (1994) that CMFGEN's `OCCUPATION_PROB` uses: per cell b = (8.589×10¹⁴ N_ion^{1/3}/n_e)^{1.5}, a = 0.0009 n_e^{1/6}/√T, x = (1+a)^{3.15} (N_ion: all ions); per level of effective quantum number n (Z the charge seen by the outer electron) K = 16n/(3(n+1)²) (1 for n ≤ 3), β = (K Z³/n⁴)^{1.5} b, f = 1.091 (x + 4(Z−1)a) β²/(7.782 + xβ), w = f/(1+f), and w = 1 for n ≤ 2Z. The LTE population of a level carries w: the Saha–Boltzmann populations, the recombination and three-body rates into the level, its share in a super-level, and the electron-impact de-excitations (by the ratio of the states' mean w, for detailed balance). Every bound-free continuum goes on below its edge into the dissolved upper levels with the threshold cross-section times 1 − w(n*), n* = √(R Z²/(ν₀ − ν)), down to 0.8 ν₀ and where 1 − w > 10⁻⁴ (CMFGEN's `QUAD_MULTI` and `PHOT_DIS_PARAMETER`): in the photoionization and recombination rates, the opacity and emissivity, the macro-atom (the absorbing level of a bound-free absorption; recombination photons below the edge with their share of the level's recombinations) and the thermal balance (a dissolved photoionization takes h(ν₀ − ν) from the gas, a dissolved recombination gives it back). The lines are not changed, as in CMFGEN. At the densities of a IIn wind (n_e = 10⁹–10¹¹ cm⁻³) only levels above n ≈ 35–60 dissolve and the Balmer continuum moves 5–12 Å to the red; in a dense shell (n_e = 10¹³, A7's has 3–5×10¹²) the levels above n ≈ 20 dissolve and it reaches ≈ 3700 Å. In the H/He benchmark (four iterations from `cmf_hhe_e4`, formal integral with 12 scattering passes) it changes b₂–b₈ by < 2 % where the Balmer lines form, the flux at 3640–3660 Å by −1 to −2.4 % and the continuum elsewhere by < 1 %, and Hα, Hβ, Pβ by +1.6, +1.3, +0.7 % (108.0, 23.0, 54.2 Å against 106.3, 22.7, 53.8 Å without it; CMFGEN `ies` 100.9, 21.1, 55.2 Å); the transport of the iterations took 10–25 % longer. In A7 (four iterations from `s10_A7p` with 2¹⁶ packets, `s10_A7d` against `s10_A7nd`) the differences stay within the scatter of the unconverged CSM (±5–10 % in the bands of the formal integral, EW5(Hα) 248 against 257 Å), and an iteration took 20–25 % longer (the statistical equilibrium 25–40 %, the transport 10–15 %). Off by default (CMFGEN's IIn benchmark models have it on). | `nltedismod.f90`, `nltexmod.f90` (`nx_lte`, `nx_el_rates`), `nltextrmod.f90` (`nxt_tables`, `bfsig`, `nxt_bffreq`), `GAS/physical_opacity_nlte.f`, `tebalmod.f90` |
| `in_nltex_jbwin` (km/s, experimental; 0: off) | Incident intensity of the Sobolev lines from the flights in a window next to each line on the side its photons drift in from (to the blue of the line for a flight whose comoving wavelength grows), `in_nltex_jbwin` wide or up to the next line with a Sobolev depth ≥ 0.1: the path-length estimator of J_ν there, where it has at least 30 effective samples (the square of the sum of the contributions over the sum of their squares), before the resonance estimator. A drift across the window adds e c/(ν\|Q\|), as a resonance crossing; where electrons scatter a photon many times while it drifts through the window (a slow, electron-scattering CSM) every flight in the window adds to the estimate, and unlike the group J without the line's own photons it counts the photons of the line that drift back into it (scattered by electrons to the incident side, or emitted in a flow whose projected velocity gradient changes sign along their path). Aimed at the decelerating CSM of A7, where the resonance estimator of Hα is starved and the line falls back to the J of its group: on the held populations of `s10_A7p` with 2¹⁸ packets (one iteration, `s10_A7z300`, window 300 km/s) the resonance estimator has a median of 6 crossings per CSM cell (at most 14), the window 15–84 effective samples (median 42), at least 30 in 43 of the 53 CSM cells, where it gives J/B of Hα 0.15–0.60 (median 0.41) against the group J's 0.06–0.59 (median 0.31). (With 2¹⁶ packets, as the A7 runs from `s10_A7c` on had through `in_niter_early`, it had 2–19 samples and was never used.) | `nltextrmod.f90` (`nxt_window_tally`, `nxt_tally_finish`), `particle_advance.f90` |
| (with the atoms) | Core-shell X-ray absorption of the NLTE ions (Verner & Yakovlev 1995; the model atoms describe the valence shell), absorbed into fast photoelectrons. | `bfxsmod.f`, `GAS/physical_opacity_nlte.f` |
| `in_rayleigh` (with `in_bbsobolev`) | Rayleigh scattering by H I in the ground state, in the scattering opacity, as CMFGEN's `INC_RAY`. `'lyman'`: the Kramers–Heisenberg amplitude A(x) = Σ_n f_n x²/(x_n² − x²) + ∫ (df/dx′) x²/(x′² − x²) dx′ (x = ν/ν₁) over the Lyman resonances (hydrogenic f_n, n = 2–200) and the 1s continuum (the exact hydrogenic cross-section), σ = σ_T A² (Nussbaumer, Schmid & Vogel 1989; Lee & Kim 2004). Next to a resonance (n ≤ 4) it is capped at 10⁴ f_n² σ_T, as CMFGEN's `RAYLEIGH_SCAT`: 1732 σ_T within ≈ 1500 km/s of Lyα, where the lines carry their cores. Above the Lyman limit it is σ_T (1 + 1.66/√x) (Gavrila 1967), and far from the resonances Rayleigh's 1.265 σ_T x⁴. Group means over the sub-groups; coherent in the comoving frame. `'kurucz'` (default): Kurucz's (1970) fit, capped at Lyα, without the Lorentz wings of the Lyman lines (the earlier treatment). `'lymes'` (benchmark; needs `in_thm_thermal`): as `'lyman'`, but redistributed by the thermal electrons, as CMFGEN's incoherent scattering (`COH_ES = F`) treats it. `'none'`: off. In the IIn wind the Lyα wing gives the H I scattering a depth of 25–40 near Lyα above the photosphere (12 at 6×10¹⁵ cm), against 0.1–3.5 for the electrons. On CMFGEN's own populations it takes SuperLite's J at 1210–1220 Å from 7.0, 5.0, 2.4 and 1.5 times CMFGEN's to 1.6, 1.2, 1.05 and 0.86 times at 1.10, 1.13, 1.65 and 3.6×10¹⁵ cm (§4.4). `test_rayleigh`. | `rayleighmod.f90`, `GAS/physical_opacity_nlte.f`, `TRANSPORT1/transport11.f90` |
| `in_tesolve` | Thermal balance with all NLTE elements: photoionization heating, recombination and collisional cooling from the model atoms, charge exchange, the LTE opacities of the remaining elements, Compton, deposition; every element re-solved at the first trial T, then those with at least 1 % of the heating and cooling terms (a held population gives the wrong T dependence of radiatively pumped levels). | `tebalmod.f90`, `GAS/te_update.f90` |
| `in_niter_early`, `in_src_n2s_early` | Fewer source packets in the first iterations, while the gas state converges. | `superlite.f90` |
| `in_xray_table`, `in_xray_ff1/2`, `in_xray_t1/2`, `in_xray_v1/2`, `in_xray_smooth` | X-ray emission of shocks embedded in a wind as CMFGEN's `INC_XRAYS`, for winds and CMFGEN benchmarks. Up to two components, with filling factor f_k, shock temperature T_k and onset velocity V_k, emit η_ν = Σ_k e^{−V_k/v} f_k² Λ_ν(T_k) n_e,X N_atom. Λ_ν comes from CMFGEN's plasma table (`RS_XRAY_FLUXES` format; APEC spectra in `atomic/misc/gal_xray_hr.dat`, linear in log T). n_e,X counts H and He as fully ionized and 6 electrons per heavier atom; the densities are those of the smooth wind (`X_SM_WIND`). The power is mechanical, added to the injected energy and not to the thermal pool. The emission per group (`xs_group_emis`) is part of the emissivity of the source functions (`grd_emnbf`): the source function of the groups in local detailed balance and of the DDMC spectra, and the emissivity of the formal integral, as CMFGEN's source function holds its X-ray emissivity. In a group thick enough for local detailed balance the X-ray photons are absorbed near where they are emitted, so the photoionizations there take them only through S: without it the rates had J = S of the gas alone, and the X-ray-irradiated H/He wind kept the He ionization of the wind without X-rays (§4.4). | `xraysrcmod.f90`, `SOURCE/sourceenergy.f90`, `SOURCE/interior_source.f90`, `superlite.f90` (`xray_emnbf`) |
| `fcl` column of `input.str` | Microclumping as in CMFGEN (Hillier & Miller 1999): the gas sits in clumps filling a fraction f of each cell's volume, with empty space between them. The gas state (EOS, NLTE populations, n_e, thermal balance, opacities and emissivities per unit volume) is that of the clump density ρ/f. The transport takes the volume means, f times these. Sobolev depths are computed with the volume-mean opacity, f K/Q, both in the statistical equilibrium and in the expansion opacity of the LTE lines, so saturated lines are not diluted. Processes ∝ ρ² (recombination, free-free, collisional cooling) are therefore enhanced by 1/f at fixed mean density, and electron scattering is not. Clumped CMFGEN models convert with `Tools/cmfgen_model.py` (RVTJ holds the clump density); clumpy CSM. | `inputstrmod.f`, `GAS/`, `nltextrmod.f90`, `superlite.f90` (`clump_scale`) |
| `in_nltex_ali` | Accelerated lambda iteration of the NLTE bound-free transitions (Rybicki & Hummer 1991, local operator), as CMFGEN converges them. In an optically thick continuum the Monte Carlo J of a cell is mostly the cell's own emission, so iterating the statistical equilibrium against it moves an ionization front by about one optical depth per iteration. For each bound-free pair, the part of J that the cell's own recombination into the level makes, Λ* η/κ, is taken out of the photoionization. The recombinations that the level re-absorbs (share n_l σ/κ) are taken out of its recombination. Λ* = 1 − (1 − e^{−τ})/τ, from the volume-mean absorption depth across the cell. The corrections are equal at the populations with which J was made, so the solution does not change. The on-the-spot cycle drops out, and the populations follow the radiation from elsewhere at once. | `nltexmod.f90` (`nx_el_rates`), `nltextrmod.f90` (`nxt_ali_cell`) |
| `in_nltex_ldb` | Local detailed balance in the groups whose volume-mean absorption depth across a cell is at least this (e.g. 30; 0: off). The Monte Carlo J of such a group is the cell's own emission, sampled by the few packets emitted in it. Where the group carries a small share of the cell's emission (the Lyman continuum of a cool wind) the estimate is mostly zero, and the on-the-spot cycle (recombination to the ground state, re-absorbed in place) becomes net recombination: case A in place of case B. There the rates use J = Λ* S + J_ext. S = (κ_th B + η)/κ is the source function of the populations with which the packets ran. Λ* = 1 − (1 − 2E₃(τ))/(2τ) is the share of its own emission that a slab of depth τ keeps. J_ext is the radiation from elsewhere: the packets on IMC steps outside the cell of their last thermal or bound-free emission (DDMC packets in a lump are local). The same Λ* takes the cell's own cycle out of the rates, as in the ALI, so the iteration has the fixed point of the exact J and needs no Monte Carlo samples of the local emission. With DDMC (`in_puretran = f`, `in_trn_tauddmc` ≤ `in_nltex_ldb`) these groups are DDMC groups, whose radiation is the diffusion limit: there Λ* = 1 and J = S (to second order in 1/τ), as no external radiation is tallied in them. With DDMC, S holds the net emission of the NLTE lines inside the group: S = (κS_c + η_l)/(κ + κ_l), with η_l = Σ βn_uA hν/(4πΔν) and κ_l = Σ β(n_l B_lu − n_u B_ul) hν/(4πΔν) over the group's lines. A line packet emitted into a thick group of a DDMC cell joins the cell's lump, which J = S would otherwise leave out. Without DDMC, the line packets count in J_ext. These are the Fe III 3d⁶ ⁵D–3d⁵4p lines at 808–866 Å, C II 858 Å and O II 832 Å inside the Lyman continuum of the full-metal IIn wind, and He I 584 Å in an H/He wind. Their photons ionize H there, and before this term the rates never saw them (§4.4). The lines themselves see the continuum alone, Λ*S_c + J_ext (`nxt_jgc`): in their incident intensity where the resonance estimator has no packets, and in the continuum term of their escape (`in_nltex_cdest`). The photons the continuum takes from a line are that line's escape. With the group's S, which holds the line's own emission, C II 904 pumped itself: in the outer wind its upper level 2s2p² ²P rose to 150 times CMFGEN's, and its photons ionized H 20 times too much. DDMC packets in a cell's lump of thick groups are shared among them as the source spectrum (S Δν), not given to the group the packet carries. `output.nltex_ldb` lists per cell the parts Λ*S and J_ext and the Monte Carlo estimate at 228–504 and 504–912 Å. | `nltextrmod.f90` (`nxt_ldb`, `nxt_lsldb`, `nxt_ali_cell`), `particle_advance.f90` |
| `in_nltex_cdest` (default t) | Continuum destruction of the trapped photons of the NLTE lines (Hummer & Rybicki 1985). A photon emitted in a thick line scatters in its core many times before it escapes the resonance zone (β). If an absorption continuum is present, the continuum takes it first with probability P(r) per emission: P(r) = (1/√π)∫ r e^(−x²)/(r + e^(−x²)) dx = −r Li₁/₂(−1/r), where r = κ_c√π v_D/K is the continuum's absorption opacity over the line's opacity at its centre (tabulated as ln(P/(1−P)) to 7×10⁻⁴; small r: (2/√π) r √ln(1/r); large r: 1 − 1/(√2 r)). The rates and the macro-atom take the escape β + (1 − β)P. Their incident intensity takes P times the continuum's J in the line's group (without the lines' own emission, `in_nltex_ldb`), so in LTE the term cancels: P S_c = P S_line. κ_c is the absorption opacity of the line's group from the previous iteration, with the populations' volume, as K. v_D includes the microturbulence. The source function of the groups in local detailed balance holds the lines' net emission with β + P (`in_nltex_ldb`). Inside the Lyman continuum of the full-metal IIn wind the Fe III 3d⁶ ⁵D–3d⁵4p lines (808–866 Å) have r ≈ 5 and C II 858 has r ≈ 0.3, while their escape from the resonance zone is ~10⁻⁶. In CMFGEN's comoving-frame transfer, H I absorbs their photons and ionizes, about 1×10⁻⁶ s⁻¹ per H particle at 1.6×10¹⁵ cm (§4.4); with β alone SuperLite kept them trapped. Lyman lines are untouched (r ~ 10⁻¹⁵ even behind the O I ¹D continuum). He I 584 Å inside the H Lyman continuum has P ≫ β (r ~ 10⁻⁶–10⁻⁴; its β ~ 10⁻⁷), the classic He I 584 → H ionization of nebulae. | `nltexmod.f90` (`nx_pcont`, `nx_el_step`, `nx_solve_cell`), `nltextrmod.f90` (`nxt_tables`, `nxt_gas_solve`, `nxt_ldb`, `nxt_pc`) |
| `in_nltex_jbbeta` (default f; experimental) | Incident intensity of the Sobolev lines in the rates, weighted by the escape probability of each direction. Sobolev's mean intensity over the line is J̄ = (1 − β̄)S + ⟨β(μ) I_inc(μ)⟩, the mean over directions of the escape probability times the intensity that reaches the resonance. The resonance estimator Σ e c/(ν|Q|)/(4πV) is the plain mean ⟨I_inc⟩, and the rates took β̄⟨I_inc⟩: right where β or I_inc is isotropic (homologous flows, a diffuse field), wrong where the radiation of a core comes along the directions of small |Q(μ)|. In a wind of nearly constant velocity the radial escape is (dv/dr)/(v/r) times the tangential one, and the core is seen radially: in the CMFGEN benchmark wind (dv/dr = 0.02–0.10 v/r beyond 2×10¹⁵ cm) the core's pumping of a thick line is 0.49, 0.27, 0.10 and 0.05 times β̄⟨I_inc⟩ at 2.2, 3.1, 6.2 and 10×10¹⁵ cm, and in a constant-velocity wind (τ₀ = 30) ⟨I_inc⟩/S grows from 1.09 at the core to 127 at 8 R_c (`Tests/benchmarks/wind`). With the weight β(μ)/β̄ = (1 − e^(−τ))/(τβ̄) on each crossing (`nxt_jbw`; β̄ from the cell's escape table with the transport's K) the estimator is ⟨βI_inc⟩/β̄, and β̄J_b in the rates is Sobolev's term: in the pure-scattering benchmarks it is S within the noise (χ²/dof 0.81–1.46 for constant velocity, a IIn-like shell and an accelerating–decelerating wind), and its terms are bounded (≤ 1/β̄; the plain ones grow as 1/|Q|). Homologous flows: the weight is 1 (to 7×10⁻¹⁰). Most lines of a wind get too few crossings for the estimator (fewer than 30 per cell and iteration: in the IIn benchmarks the Lyman lines, Hα beyond 1.1×10¹⁵ cm and most metal lines), and the J of their group stands in for it, which has the same fault. The Monte Carlo therefore tallies the path length of every flight in eight bins of |μ| (edges 0, 0.5, 0.75, 0.88, 0.94, 0.97, 0.985, 0.995, 1; analytic along the straight flight, μ = z/√(p² + z²)), and the stand-in is weighted by Σ_b f_b β_b/Σ_b Δμ_b β_b, with f_b the group's share in bin b (the DDMC and local-detailed-balance part isotropic) and β_b the bin means of the escape probability (per-cell tables against K, `nx_betab_table`, `nx_bfac`): 1 for an isotropic field or a thin line; for a core seen within μ > 0.97 exact to 10⁻⁴ for any depth, within a bin to 7.5 % (`test_nltex`). In the constant-velocity benchmark (τ₀ = 30) the incident intensity the rates take, estimator or weighted group J, is 0.80–1.29 times S in every cell out to 10¹⁶ cm; before, it was 2.4–3.6 times S at 1.2–2×10¹⁵ cm and 22–112 times beyond 4×10¹⁵ cm. Off by default: in the CMFGEN H/He wind it takes b₁ further from CMFGEN's comoving-frame solutions, not closer. At 2.5×10¹⁵–1.5×10¹⁶ cm b₁ with it is 1.4–4.6 times `blk`'s; without it (`rayGq_hhe`) 0.86–1.5 times, and 0.9–1.3 times CMFGEN's Sobolev `sob` from 2.5×10¹⁵ cm out (§4.4). Probably because in a comoving-frame solution the core radiation that reaches a resonance has been processed by the electron scattering and the lines on the way, which the weighting of an isolated line leaves out. | `nltextrmod.f90` (`nxt_jbw`, `nxt_btab`, `nxt_bfac`, `nxt_jmsum`), `nltexmod.f90` (`nx_betab_table`, `nx_bfac`), `TRANSPORT1/transport11.f90`, `particle_advance.f90` |
| `in_nltex_dmpsec` (default f) | Adaptive damping of the NLTE populations between iterations. With the fixed weight `in_nlteh_damp` = 1/2 of the new solution (geometric mean of the old and new populations), an iteration whose solve answers a change of its input with −3 times that change (slope g = −3) oscillates with period 2 and does not decay. That is the full-metal IIn wind with continuum destruction (`in_nltex_cdest`): H and the metals that absorb in the Lyman continuum are strongly coupled, and the change between solves stayed at 1–4 (the largest relative change of a state above 10⁻⁴ per cell). With this option each cell takes the weight 1/(1 − g): g is the secant slope of the solve's output against its input along the last step, projected on the change of the log populations of the states above 10⁻⁴. That cancels the dominant mode of the iteration; a monotone approach (0 ≤ g < 1) keeps the weight 1. The weight is clipped to [0.1, 1] and averaged geometrically with the cell's last weight (Monte Carlo noise). Cells without a history, or with a step below 5 %, keep their last weight. The log prints the median and smallest weight. Each rank saves the history (weights and the log populations of the last solve's input and output) every iteration in `output.nltex_dmp_<rank>`; a run continued from its state with `in_nltex_init` reads those files from the same directory, matched by cell index (same cells and model atoms only), and keeps the weights instead of restarting at 1/2. Tested on toy01: a 2-iteration continuation of a 3-iteration run restored the history of all 100 cells; median weight 0.84 instead of 0.5. In the full-metal IIn wind with continuum destruction (`rayFs` → `rayFs2`) the change between solves fell from 2.1 (fixed weight 1/2: a period-2 cycle, x(H I) of cell 45 alternating between 10⁻³ and 1.5×10⁻²) to 0.58, 0.16 and 0.076 at a median weight of 0.2, and x(H I) of cell 45 settles at 3.3×10⁻³ (CMFGEN 9.3×10⁻⁴: the inner-wind offset of §4.4, not the iteration). Use it with `in_nltex_cdest`. Off by default because a monotone approach takes the full step (weight 1) where `in_nltex_damp` takes half. | `nltextrmod.f90` (`nxt_gas_solve`, `dmp_weight`, `nxt_write_dmp`, `nxt_read_dmp`) |
| `in_nt_sf`, `in_nt_sfemax` = 1000, `in_nt_sfde` = 1; `in_nt_lumg` | Fast electrons (from X-rays above `in_nt_emin`, or γ-rays and positrons) share their energy among heating, ionization and excitation of all NLTE species by the Spencer–Fano degradation (Kozma & Fransson 1992), as CMFGEN's non-thermal solver (Li, Dessart & Hillier 2012), in place of Shull & van Steenberg's shares for H I and He I. One electron is injected at `in_nt_sfemax` eV on a uniform grid of `in_nt_sfde` eV and degraded from the top down. A finished bin pushes its collisions to the bins the electrons land in (linear split), which conserves the number of electrons exactly. Secondaries follow Opal et al. (1971) with width 0.6 I. The Coulomb loss is n_e(2πe⁴/E) ln(2E/ħω_p) (Kozma & Fransson's form, as CMFGEN). Ionization acts on the levels of the ground term of each ion. It uses the subshell fits of Arnaud & Rothenflug (1985), read as `nt_ar.dat` from the model-atom directory (`Tools/cmfgen2superlite_nt.py` writes it from CMFGEN's `arnaud_rothenflug` data), or Lotz (1967) with the subshell thresholds of Verner et al. (1996) for ions without fits (Co). Excitation runs from every level above 10⁻⁴ of its element (in ions above 10⁻³) to every level above with a line or a tabulated collision strength. Allowed lines (A > 10³ s⁻¹, or no collision strength) use van Regemorter with the neutral or ion Gaunt factor; the others use the collision strength at the table's highest temperature, constant in energy. Ionization routes of ions below 10⁻³ of the most populated target are left out, keeping at least the ten most populated ions. Each cell is solved from its populations once per population solve and per thermal balance. The rate matrix gets the ionization rate of each state and the excitations between states (`nx_el_matrix`); the thermal balance gets the heat share. The Monte Carlo deposits a fast electron into heat, the continuum state of an ionized ion (energy I) or the upper state of an excitation (energy ε) by the cell's cumulative shares (`nxn_pick`, gathered across ranks). `in_nt_lumg` is the share of the `lumg` energy of each cell deposited by fast electrons: the decay deposition over the injected energy, e.g. 0.75 for the StaNdaRT toy01 at 19.5 d, where L_bol is 1.32 times the deposition. Those source packets enter through the same channels, and the cell's thermal heating keeps only their heat share. `output.nltex_nt_<rank>` lists per cell the deposited power and the shares of heat, ionization and excitation. | `nonthermmod.f90`, `nltentmod.f90`, `nltexmod.f90`, `nltextrmod.f90` (`nxt_nonthermal`), `GAS/te_update.f90`, `SOURCE/interior_source.f90` |
| `in_nltex_vturb`, `in_nltex_vturbmax` (km/s; the scattering with `in_thm_thermal`), `in_nltex_esrates`, `in_nltex_esdamp` | Electron scattering of line photons trapped in their Sobolev resonance zones. A photon of a thick line is absorbed and re-emitted about 1/β times before it escapes, while its comoving frequency drifts across the line core, |x| < x_c ≈ √ln(τ/√π) Doppler widths. That is a path of several Sobolev lengths v_D/\|Q\|, on which an electron can scatter it out of the line with a Doppler shift of hundreds of km/s: to the red it escapes, to the blue the flow brings it back to the line. Point resonances (the Sobolev approximation) have no such path. In a IIn wind σ_es v_D/\|Q\| ≈ 0.1–0.4 where Hα forms, so the trapped scattering is as likely as a scattering on the whole free path, and it builds the red shoulder of the Balmer lines (CMFGEN's comoving-frame transfer of the observer spectrum includes it). Per emission its probability is β_es = σ_es v_D ⟨y(τ_μ)⟩_μ/K × exp(−K Γλ/(4π²⟨\|Q\|⟩v_e)), with y = τG(τ) for complete redistribution over a Doppler profile (`nx_ytrap`; tabulated, and checked by a Monte Carlo of a resonance in `Tests/unit/test_estrap.f90`). The exponential is the escape of the scattered photon, shifted by v_e = √(2kT/m_e), through the Lorentz wings of the line (radiative damping Γ; off with `in_nltex_esdamp = f`, Doppler profiles as CMFGEN's iterations): Lyα and the low Lyman lines of a dense wind keep their photons. v_D² = 2kT/m + v_turb², with v_turb linear in v from `in_nltex_vturb` (v = 0) to `in_nltex_vturbmax` (outer edge), as CMFGEN's VTURB_MIN/MAX. By default the formal integral adds the emissivity n_u A β_es hν/4π of every line, redistributed with the angle-averaged thermal-electron kernel, to the resonances, with the populations of the iterations, which do not include it: as CMFGEN's observer spectrum adds the incoherent scattering to iterations with coherent scattering. CMFGEN's iterations with incoherent scattering change n₃ by < 10 % and Hα by 4 % in the H/He benchmark, so the populations need not see it. With `in_nltex_esrates` (experimental) the statistical equilibrium, the macro-atom and the transport have it as well. A line de-activation emits, with probability β_es/(β + β_es), the photon scattered by a thermal electron (comptonmod kernel about the line centre, new direction), and the rates and the macro-atom have the escape probability β + β_es. A photon that an electron scatters into the Doppler profile of a Sobolev line (in the transport, and the scattered line emission at its first step) is absorbed there with the depth of the profile over its further path, τ = (K/2\|q\|) e^(−x²) erfcx(a/2 − x), a = σ_es v_D/\|q\|, x the offset in Doppler widths, positive on the side it drifts in from (`nxt_esfar`): the Sobolev depth K/\|q\| ahead of a line and (K/\|q\|) erfc(\|x\|)/2 behind it without the scattering, and the profile against the electrons, K φ(x)/(σ_es v_D), where electrons scatter it many times per Doppler width of drift; the lines ahead have then been tested for its whole path, and it passes their resonances until it scatters again. The rates take this absorption from its tally per line and cell, as the incident intensity (β J_b + E c/(4πVνK))/(β + β_es) (the static lines of a slow cell: β_es J_es, J_es the group J averaged over the electron kernel). Where the tally has fewer than 30 events its expectation stands in for it, σ_es v_D ⟨W(K/\|Q(μ)\|)⟩_μ J_es/K (`nx_bew`), with W(t) = ∫₋₅⁵ (1 − e^(−t erfc(−x)/2)) dx the window in which `nxt_esfar` absorbs where the flow drifts the photon across the profile before it scatters again (tabulated, `nx_wwin`; `test_estrap`). W is wider than the core of β_es (W/y = 1.7 at t = 200, 1.1 at 10⁸) because the photons that land ahead of the core drift into it and pass the resonance. In the outer H/He wind Lyα gets a few hundredths of a packet per cell and iteration, Lyβ–Lyζ none, and even Hα has 30 events only below 1.075×10¹⁵ cm (where the tally is 0.8–1.36 times β_es J_es): before, those lines had the escape β_es without its inflow (a zero tally), or once in a while tens of times its mean, and x(H I) of the outer wind swung by 2.4 between iterations (`rayGe_hhe`); with the expectation by at most 15 % (`rayGf_hhe`, five solves from `rayG_hhe`). Its b₁ at 2.5, 3.6, 5.9 and 9.5×10¹⁵ cm: 8.9, 77, 807, 7110, against 6.0, 52, 583, 4150 without `in_nltex_esrates` (`rayG_hhe`) and CMFGEN's 4.2, 37, 312, 1620 (`ies`) and 3.9, 47, 496, 2650 (`blk`): the trapped-photon scattering in the rates makes the outer wind 1.4–1.7 times more neutral, away from `ies`. Redistributing the Rayleigh-scattered photons of the Lyα wing with the thermal electrons as well, as `ies` does (`in_rayleigh = 'lymes'`, `rayGel_hhe`, six iterations from `rayG_hhe`), raises b₁ there by another 1.3–1.8 (2.8–5.7 times `ies`' at 2.5–9.5×10¹⁵ cm), so `in_nltex_esrates` stays off by default. In the H/He benchmark (six iterations from the Sobolev state, Doppler profiles in the escape as in CMFGEN's iterations, `in_nltex_esdamp = f`; `cmf_hhe_esy`) the layer just above the photosphere ionizes as in `ies` but much less (x_HI ≤ 0.15 at 1.06–1.10×10¹⁵ cm, 0.20 without it, 0.006 in `ies`), the Balmer continuum comes 7–10 % closer to `ies`', and the formal integral with the Monte Carlo J gives Hα, Hβ, Pβ 103.4, 25.2, 59.1 Å (without it 91.3, 17.3, 55.4; `ies` 100.9, 21.1, 55.2); but b₁ of the outer wind rises to 3.7 times `ies`' (1150 against 311 at 5.9×10¹⁵ cm; 366 without it) and the converging Balmer series goes into emission (3650–4000 Å 17 % above `ies`, 3 % below without it). In the full-metal wind (four iterations from the hybrid, `cmf_m002re_dm16e`, formal integral with 12 scattering passes) Hα, Hβ, Pβ are 87.4, 15.8, 56.7 Å against CMFGEN's 87.2, 17.5, 50.5 (without it 98.5, 19.8, 52.7), the optical continuum 2–4 % brighter relative to CMFGEN's than without it, and b₁, b₃ of the outer wind 15 and 7 % lower (where CMFGEN's are already 6 and 1.2 times SuperLite's). Earlier versions, which took the incident intensity of the scattered photons from the group J and let only the Sobolev resonances absorb them, over-ionized the H/He wind 2–150 times (Hα 198 Å) and put the converging Balmer series of the full-metal wind into emission at 2.6 times CMFGEN's flux. Sharing a line's photons between the line and its wings in the transport alone (rates with β) keeps the populations but takes Pβ from 56 to 34 Å. The red shoulder depends on the microturbulence (CMFGEN H/He wind, EW of +300–2700 km/s on the air scale: 21, 25 and 38 Å for 1, 10–20 and 50 km/s; §4.4). The microturbulence also sets the Doppler widths of the lines in the formal integral (`in_fi_doppler`), and with `in_fi_nes` the trapped photons come from the iterated J of the formal integral instead of β_es. | `nltexmod.f90` (`nx_ytrap`, `nx_es_table`, `nx_bes`, `nx_wq`, `nx_jes`), `nltextrmod.f90` (`nxt_esfar`, `nxt_efsum`), `TRANSPORT1/nxemitmod.f90` (`nx_linees`, `nx_lineabs`), `TRANSPORT1/transport11.f90`, `particle_advance.f90`, `formalmod.f90` |
| `in_nltex_dopmc`, `in_nltex_dopz`, `in_nltex_dopes` (experimental) | Static Doppler profiles of the NLTE lines in the Monte Carlo transport of slow cells, for the overlapping line forest of a slow wind. Where the velocity changes across a cell by less than `in_nltex_dopmc` Doppler widths (velocity span and \|v/r\|Δr), a Sobolev resonance would lie beyond the cell. There a packet sees the active lines of the elements from Z = `in_nltex_dopz` (default 3) within 4 Doppler widths of its comoving frequency as opacity K/(√π v_D) e^(−x²) per unit length, except the line it was emitted in, whose trapping stays with the Sobolev escape of the macro-atom. Its flights end after a drift of half a Doppler width, so the opacity follows the comoving frequency and a drift through a profile has the Sobolev depth K/\|Q\|. An absorbed packet activates the macro-atom in the upper level of the line chosen by its share of the opacity, and line photons leave along the directions of their Sobolev escape. The incident intensity of these lines is the path-length estimator over a window 4–6 Doppler widths from the line, on the side the comoving frequency drifts in from (blue where the projected velocity gradient of the flight is positive); a drift through the window gives the Sobolev crossing e c/(ν\|Q\|). Inside the profile of a thick line the photons drifting in are absorbed in the wing, so the intensity there is about 1/τ of the incident one: the first version estimated it there, and the incident intensity of Hα in the H/He wind fell to 1 % of the Sobolev estimate, that of Hβ to 3 %. H and He lines keep their Sobolev resonances in slow cells (`in_nltex_dopz` = 3). In a dense wind an electron scatters a photon about as often as it drifts across a Doppler width. A static profile absorbs every photon scattered into the core of a thick line, while the macro-atom's Sobolev escape treats each absorption as a whole trapping episode, with no escape by electron scattering. Static H and He lines (`in_nltex_dopz` = 1) thus took the electron-scattering wings out of the Balmer lines of the H/He wind (EW of Hα 82 Å and Hβ 8 Å, against 102 and 22 Å with resonances and 101 and 21 Å in CMFGEN's comoving-frame model), cut the wind's Lyman-continuum J to a quarter and raised b₁ fourfold. With resonances, the scattered photons that point resonances miss and the absorptions that the Sobolev macro-atom over-converts appear to offset each other: that configuration agrees with CMFGEN. `in_nltex_dopes` (default t; needs `in_nltex_vturb` ≥ 0 and `in_thm_thermal`) gives the static lines of the slow cells the electron scattering of their trapped photons in the rates and the macro-atom, as `in_nltex_esrates` does for all lines: escape probability β + β_es, incident intensity (β J_b + β_es J_es)/(β + β_es), and the scattered emission (em_lines). The transport scatters photons into the static cores explicitly, and without the escape the macro-atom converted every one of them: in the full-metal benchmark the hybrid moved 2.6×10⁴¹ erg/s from the UV to the optical (1100–2000 Å at 0.52 of CMFGEN's luminosity, 5000–10 000 Å at 1.27–1.30); with it 1.03 and 1.07–1.09 after ten iterations, as with coherent scattering, and against CMFGEN's model with incoherent scattering in its iterations the spectral shape agrees within 4 % from 1100 Å to 2.5 µm but at 3000–3650 Å, Hα, Hβ and Pβ within 12, 14 and 4 % (§4.4). With it the static H and He lines of the H/He wind keep their wings (Hα 93, Hβ 18.5, Pβ 49 Å after four iterations, `cmf_hhe_dm5d`, against 82, 8, 38 Å without it), but b₁ of the outer wind stays 2.5–4 times CMFGEN's (84–92 against 37 at 3.6×10¹⁵ cm) and b₃, b₄ 10–30 % below, as without it: H and He keep their resonances by default. In the full-metal wind static H and He lines with this escape (`cmf_m002re_dm10s`) put the converging Balmer series into emission as well (1.36 times CMFGEN's flux at 3650–4000 Å). `test_dopmc`. | `nltextrmod.f90` (`nxt_static_line`, `nxt_static_tally`, `nxt_sobolev_neighbors`, `nxt_useb`, `nxt_es_cell`), `nltexmod.f90` (`nx_esm`), `TRANSPORT1/transport11.f90`, `TRANSPORT1/nxemitmod.f90`, `particle_advance.f90` |
| `in_nltex_overlap` (km/s; 0: off) | Overlapping lines in the Sobolev rates and the macro-atom. A photon trapped in the resonance zone of a thick line is absorbed also by the lines whose Doppler profiles reach its own: Lyβ's by O I 1025.76 Å, 11.7 km/s to its red (the Bowen fluorescence of O I 8446 and 11 287 Å); the converging members of a series; blends of multiplets. Point resonances keep each line's trapped photons to itself. The lines within this separation are a line's partners. Per cell the code solves the Sobolev zone of the line and those partners whose profiles reach it (`nx_ovl_zone`): complete redistribution over Doppler profiles (thermal and `in_nltex_vturb` widths), with the comoving frequency drifting to the red by \|Q(μ)\| per unit length. Each grid cell of a quarter Doppler width absorbs a photon in proportion to each line's depth there, averaged over the directions. Per emission this gives the escape from the zone, ρβ (ρ = 1 for a line alone; a blend of two equal lines escapes as one line of their summed depth), and the probability P_kj that partner j absorbs the photon. The photons line k gives to a partner j of the same element, n_u(k) A_k P_kj, are linear in k's upper population. They de-excite k and excite j in the element's rate matrix, so a blend converges like a single line. Partners of other elements enter through k's escape probability, ρβ + Σ_j P_kj, and their photons that k absorbs, Σ_j K_j S_j P_jk/K_k, with their populations of the moment; that coupling is one-sided (Lyβ and O I). The zone's incident intensity (that of its bluest line) is absorbed with the factor of detailed balance, ρβ + Σ_j P_kj − Σ_j K_j P_jk/K_k over all partners, so LTE stays exact. The macro-atom escapes the zone and hands the net photons of each pair to the partner's macro-atom (a channel that activates it). Zones are made at the first iteration of a cell's statistical equilibrium and again every tenth while a population has moved by more than 10 % since. Partners thinner than 10⁻³ for the mean projected gradient stay out of a zone. The thermal balance keeps the cell's zones for its trial temperatures. On cells of the full-metal benchmark that costs 0.14–0.19 s per cell (48–55 iterations, against 0.02–0.04 s and 12–37 without the overlap). With the partners of the element lagged as well, the blends of He I, C, N and Fe ran into the iteration limit. In the IIn wind O I absorbs a trapped Lyβ photon with probability χ(O I)/χ(Lyβ) ≈ 1.3×10⁻⁴ per emission, 300–3000 times its escape probability. So O I converts it long before it escapes, n = 3 empties into the ground state, and O I 3d ³D° is pumped (§4.4). `test_overlap`. | `nltexmod.f90` (`nx_ovl_init`, `nx_ovl_zone`, `nx_ovl_cell`, `nx_ovl_terms`), `nltextrmod.f90` (`ovl_tables`, `em_ovl`), `tebalmod.f90` |
| `in_nltex_init` | Warm start: the `output.nltex_state` (written after every iteration, so a run that stops can be continued) of a run on the same structure and model atoms (populations, and T_e with `in_tesolve`, matched by exterior mass as `in_nltex_prev`), and its Monte Carlo radiation field, `output.nltex_jmc` next to it (the J of the groups, averaged over this run's groups, and the incident J of the lines; written with the state). The first iteration keeps the populations instead of solving with J = B(T_R) and makes its macro-atom tables with the saved J, so it is an iteration of the saved state: a rerun with other output settings (formal integral, its line profiles) needs one iteration. Without the J file (runs before it existed) the first iteration's tables have J = B(T_R) and no line incident intensity: the line re-emission of its transport differs and its Monte Carlo J is not that of the state (the Hα electron-scattering wings of a formal integral with it were a third to a half of the converged ones); do not average it into the formal integral. `in_nltex_hold` (diagnostic) keeps the saved populations in every iteration, without statistical equilibrium: the transport and the formal integral on given populations, e.g. those of a CMFGEN model (`Tools/cmfgen2superlite_state.py`: the levels of its POP files summed into the super-levels of SuperLite atoms converted from the same data, averaged over each cell between CMFGEN's depths). That separates the radiative transfer from the rate equations (§4.4). | `nltextrmod.f90` (`nxt_read_init`, `nxt_read_jmc`, `nxt_write_jmc`, `nxt_gas_solve`, `nxt_hold`), `superlite.f90`, `Tools/cmfgen2superlite_state.py` |
| `in_fi`, `in_fi_wl1`, `in_fi_wl2`, `in_fi_dv`, `in_fi_npc`, `in_fi_navg`, `in_fi_doppler`, `in_fi_nes`, `in_fi_nesdv`, `in_fi_only` | Formal integral of the emergent spectrum after the last iteration, as CMFGEN's `cmf_flux`: noise-free line profiles at the spacing `in_fi_dv` [km/s] between `in_fi_wl1` and `in_fi_wl2` [cm] (`output.fi_spectrum`: λ [Å], L_λ, L_ν). The observer-frame invariant I/ν³ is integrated along rays of impact parameter p through the spherical grid, with the comoving source functions of the last iteration (first order in v/c). Continuum: absorption and thermal emission of the LTE opacities, the recombination emissivity of the NLTE levels, electron scattering of the Monte Carlo J redistributed by the thermal electrons with the angle-averaged kernel of Thomson scattering (isotropic, or dipole with `in_trn_thm_aniso`), the scattering part of the LTE expansion opacity. The Monte Carlo noise of J is that of the scattering emissivity, so of the continuum where electron scattering dominates (≈ 3 % per 100 km/s group with 2¹⁸ packets): `in_fi_navg` averages J over the last iterations of a converged model. NLTE lines: the source function of the level populations (complete redistribution) with Doppler profiles, the thermal width of the ion and the microturbulence of `in_nltex_vturb` (`in_fi_doppler`, default): opacity K/(√π v_D) e^(−x²) per unit length, which gives the Sobolev depth K/\|Q(μ)\| where the flow crosses the profile quickly and, in a slow flow, the depth along the path through the profile (the narrow P-Cygni absorption of a coasting wind, §4.4); a line whose profile a ray crosses within a tenth of the cell width is a point resonance there (fast flows: supernova ejecta, at the cost of the Sobolev integral); with `in_fi_doppler = f` Sobolev resonances wherever the comoving frequency crosses a line. With `in_nltex_vturb` also the electron scattering of the photons trapped in the lines (above), added where they are scattered first. `in_fi_nes` > 0 (with the Doppler profiles) iterates the electron scattering instead, as CMFGEN's `cmf_flux`: the comoving J of each cell on the grid of the spectrum (extended by the kernel and the flow) starts from the Monte Carlo J; each pass integrates the rays of all frequencies with the scattering emissivity of the previous J and tallies the new J along them (parallel rays sample every direction of a spherical shell: J = V⁻¹ Σ_p 2πp dp ∫ I dz, core rays also behind the inner boundary, with V the volume the rays sample, Σ_p 2πp dp times their chords in the cell, so that a uniform I gives J = I: the shell volume differs from it by 0.2–7 %, by 70 % in a thin innermost cell, an excess the scattering amplifies pass by pass). Inside a line J approaches its source function, so the trapped photons and their further scatterings through the flow (with its Doppler shifts) are in J, and the scattering emissivity has no Monte Carlo noise; the trapped-photon terms above are then left out. Each pass costs a formal integral; Ng's extrapolation every third pass; 12 passes converge the H/He IIn wind (τ_es ≈ 5) to 0.07 % of L and 0.05 % of the J of any cell per pass. `output.fi_jes`: per cell the frequency-integrated J of the last pass over the Monte Carlo J. `in_fi_nesdv` [km/s] runs the J and the rays of all passes but the last two on a coarser grid (the scattering emissivity is smooth over the electron kernel) and the last two on the grid of the spectrum: the cost of 10 passes falls by the ratio of the spacings. H/He benchmark with 60 km/s against 20 km/s: L and the continuum within 0.1 %, line cores unchanged, the electron-scattering wings of Hα and Hβ within 5 % (±0.25 Å per bin) — in the slow inner wind a ray samples one comoving frequency, so a coarse grid of rays samples the line cores unevenly, and the wing J of the coarse passes converges too slowly for the last two to correct it: for exploration, not for final spectra. `in_fi_only` computes only the formal integral of the saved state and Monte Carlo J of `in_nltex_init` (no transport; one iteration): a rerun with other spectrum settings, as CMFGEN's `cmf_flux` after a model. Rays through the inner boundary (`in_fi_npc`, uniform in p²) start with the Monte Carlo J of the innermost cell; two per cell outside. L_ν = 8π² ∫ I p dp; the frequencies are shared among the MPI ranks. Line wavelengths are in air above 2000 Å (row 9 of §1), so a CMFGEN observer spectrum (c/ν: vacuum) must be converted to air before profiles are compared (83 km/s at Hα). | `formalmod.f90` |
| `in_nltex_prev`, `in_nltex_dtprev` (with `lumg`) | Time-dependent statistical equilibrium of the NLTE elements over a sequence of models (Dessart & Hillier 2008, CMFGEN's time-dependent terms). Each model writes the populations of all states per cell (`output.nltex_state`). A later model reads them, matched by exterior mass, which is Lagrangian when both models end at the same mass shell, and solves dx/dt = M x per element as a backward-Euler step over dt: (M − 1/dt) x = −x_prev/dt, with particle conservation exact. Recombination therefore lags behind a falling radiation field, and ionization behind a rising one. The change of ionization and excitation energy, Σ n_el ε_s (x_prev − x)/dt, is part of the cell's injected energy, since the hydro L(r) contains what the gas releases. The macro-atom emits it, activated in the states that empty, so the stored ionization comes out as recombination lines and continua (the Hα of SNe II on the plateau). It is added where it exceeds the injected energy, and taken from that energy where the gas ionizes; the thermal pool heats only on the rest. | `nltexmod.f90`, `nltextrmod.f90`, `SOURCE/sourceenergy.f90`, `SOURCE/interior_source.f90`, `GAS/te_update.f90` |

Validation (`Tests/unit/run_tests.sh` with `SL_ATOMS` = the converted atoms):

| Test | Checks | Result |
|---|---|---|
| `test_nltex` | built-in H against the stage-1 H atom; Case B; the He I macro-atom and the H macro-atom against the statistical-equilibrium emissivities; every converted atom (14 elements, 44 ions, 2005 states; with N V, O V–VI and super-levels: 450–1100 states) in a Planck field at high density at 8000, 25 000 and 10⁵ K; charge exchange at LTE | rates to 2×10⁻³, LTE populations to 7×10⁻⁵; Hα/Hβ = 2.74 (Storey & Hummer 2.86 l-resolved), He I 5876/4471 = 2.53; macro-atom He I 2.34 ± 0.03 against 2.37, Hα/Hβ 3.15 ± 0.02 against 3.14; all elements Saha–Boltzmann to < 2 %; no net charge-exchange flow into any state to 10⁻¹⁴ |
| level dissolution (`test_nltex`) | occupation probabilities against CMFGEN's `OCCUPATION_PROB` formula (H at n = 8–30 and He II at n = 20, n_e = 10¹³ and 10¹⁵ cm⁻³, values computed independently); built-in H in a Planck field at 10⁴ K, n_H = 10¹⁵ cm⁻³ (n_e = 4.4×10¹⁴, 1 − w(n = 12) = 0.36); every converted atom at 2.5×10⁴ K, n_e = 10¹⁷; the bound-free level choice of the transport just to the red of the Lyman, Balmer and Paschen edges (915, 3700, 8400 Å) | w to 2×10⁻⁷; the statistical equilibrium returns the dissolved Saha–Boltzmann populations to 1.3×10⁻⁵ (H) and 2.8×10⁻⁴ (all elements); the level choice follows n_l σ_l with the dissolved continua (χ² 35, 3.4, 26 for 29, 5, 29 degrees of freedom) |
| bound-free level choice (`test_nltex`) | all atoms at LTE (25 000 K): the chosen levels at 3000, 900, 500 and 200 Å against n_l σ_l, with the default lists (2×10⁵ choices) and with lists of 4 levels (5×10⁴; a walk of the other levels in 21 % of the choices) | χ² = 39.8, 4.1, 45.1, 167.8 for 30, 4, 62, 182 degrees of freedom; with lists of 4: 25.3, 2.4, 33.3, 89.9 for 18, 2, 28, 105 |
| formal integral (`formalmod.f90`) | m002 (§4.4, IMC everywhere) with `in_fi` over 3500–13 500 Å at 20 km/s, against the Monte Carlo spectrum of the same iteration | L of 3500–5500, 5500–7500, 7500–13 500 Å within +3.1, +1.7 and +0.9 % of the Monte Carlo spectrum, which it follows without its noise |
| X-ray source (`xraysrcmod.f90`) | emitted X-ray luminosity of CMFGEN's ζ Pup model (f = 0.078 at 3 and 5 MK, onset 500 km/s, smooth-wind densities, microclumped wind) against CMFGEN's `OBSFLUX` | > 0.1 keV: 3.80×10³³ against 3.77×10³³ erg/s; > 1 keV: 3.76×10³² against 3.75×10³² erg/s (1 % and 0.3 %) |
| ALI (`test_nltex`) | built-in H at 10⁴ K, n_H = 10¹⁰: at the populations with which J was made the solution must not move; a lambda iteration of the Lyman continuum (J = J_ext + Λ* S, Λ* = 0.999) with and without ALI | fixed point kept to 3×10⁻⁷; converged in 2 iterations instead of 68, x(H II) the same to 7×10⁻⁶ |
| Continuum destruction (`test_nltex`) | P(r) against a quadrature of (2/√π)∫₀^∞ r/(1 + r e^(x²)) dx at r = 10⁻¹⁴–10¹² (on and between the table's points); H at 10⁴ K, 10¹⁷ cm⁻³ in LTE (J = B) with a continuum of 10⁻⁵ cm⁻¹ in every group; the Case B nebula with a continuum that takes every trapped photon where J = 0, against the same with a velocity gradient that makes β = 1 | P to 2×10⁻⁴; LTE kept to 7×10⁻⁹; P → 1 identical to the thin lines (x(H II) 0.94989) |
| Local detailed balance (`test_nltex`) | slab share Λ* against E₃ by quadrature (τ = 0.1–30); J of a one-cell grid (thick groups, DDMC lump share; with DDMC Lyβ inside a thick group: S = (κS_c + η_l)/(κ + κ_l) from the populations and β); the toy cell of the ALI test with a thick Lyman continuum (cell 10¹³ cm): J = J_ext + Λ*S, with the Monte Carlo estimate missing the cell's own emission in the groups above τ = 30 | Λ* to 1×10⁻⁷; J exact (line term to 2×10⁻¹⁶, moving J from S_c towards the line's source function); x(H II) 0.816 with the estimate alone (case A), 0.97447 for the exact J (125 lambda iterations) and with local detailed balance (4 iterations), the same to 10⁻¹⁴ |
| Resonance-estimator weights (`test_nltex`) | β(μ)/β̄ (`nxt_jbw`) in a homologous flow; its mean over directions for an isotropic field in flows of constant velocity, acceleration, deceleration and a gradient that changes sign over directions (K/max|Q| = 2×10⁻⁴–10⁷); its mean over the cone of a core (μ > 0.9) in a constant-velocity wind with a thick line, against (1 − μ_c)(2 + μ_c)/2; the weighting of the group J by its bins of |μ| (`nx_bfac`) for a core within μ > 0.97 (a bin edge) and μ > 0.9 (inside a bin), dv/dr = 0 and 0.05 v/r, K/qv = 10⁻³–10⁴, against the exact ⟨β⟩_core/β̄, and for an isotropic field | 1 to 7×10⁻¹⁰; unbiased to 6×10⁻⁴; 0.1450 (expected 0.1450); 8.6×10⁻⁵, 7.5 % and 0 |
| DDMC through the macro-atom (H/He wind, §4.4) | CMFGEN's H/He IIn-like wind (Sobolev lines, converged, radiative-equilibrium T_e), SuperLite at CMFGEN's T_e with DDMC (`in_trn_tauddmc` = 5) and local detailed balance, 6 iterations from the same state; Monte Carlo spectrum of the last two iterations against CMFGEN's observer-frame spectrum | EW of Hα 27 → 101 Å (CMFGEN 96), Hβ 0.7 → 21 Å (21), Pβ 8 → 29 Å (49; the Monte Carlo IR is noisy, 13–40 Å per iteration); the H I fraction stays within a factor 3 of CMFGEN's. With the continuum fallback of the lines' incident intensity: Hα 82, Hβ 16, Pβ 48 Å, the formal integral 81, 20, 43 Å. IMC everywhere: Hα 84, Hβ 16 Å, but H 3–200 times more ionized than in CMFGEN |
| time-dependent SE (`test_nltex`) | built-in H in a Case B nebula: the steady state with itself as the previous model, dt → 0, and one backward-Euler step of recombination in the dark (10⁴ K, n_H = 10⁴ cm⁻³, dt = 10⁸ s) against the root of x = x₀/(1 + α n_H x dt) | steady state kept to 1.4×10⁻⁷; frozen to 4×10⁻¹²; x(H II) = 0.7627 against 0.7626 |
| `test_tebal` | Planck field with LTE opacities, H, He I and He II must return T₀; pure-H nebula; fast-electron shares | T₀ to 0.2 %; 21 500 K; Shull & van Steenberg values |
| `test_superlevel` (`SL_ATOMS_SL`, `SL_ATOMS_FULL`) | Fe II–IV (400, 300, 300 levels) with CMFGEN's super-levels (136, 55, 51) against the same atoms with every level a state: photospheric, diluted nebular (W = 0.1), dense and cool cells | Fe ion fractions and Fe II level populations to ≤ 7 % (nebular), ≤ 2.5 % otherwise; the same atoms with super-levels and the CMFGEN-matched C, N, O, Ca, Fe set pass the LTE and charge-exchange checks of `test_nltex` |
| `test_dopmc` | static Doppler profiles of slow cells: the opacity of a line at its centre and one Doppler width off, of two overlapping lines, none beyond 4 widths, none of the packet's own line or of an H line (Z < `in_nltex_dopz`), which the resonance search of the cell finds instead; the absorbing line by its share of the opacity; the incident intensity of a drift through the window outside the profile, on the blue side for Q > 0 and the red side for Q < 0, against the Sobolev crossing e c/(ν\|Q\|); nothing from inside the profile, of the own line or of the H line | opacities to 2×10⁻¹²; share 0.66670 against 2/3 (3×10⁴ draws); the window against the crossing to 10⁻¹⁵ |
| `test_rayleigh` | H I Rayleigh cross-section against an independent Kramers–Heisenberg sum (Python/scipy, the continuum term by adaptive quadrature) at 1050–5000 Å; the static-polarizability limit 1.2654 σ_T x⁴; the cap at Lyα (1732.2 σ_T); Gavrila's fit; a group mean | to 5×10⁻⁶ |
| `test_overlap` (`SL_ATOMS_OVL`: atoms with `h1.atom`, `o1.atom`) | Sobolev zone of overlapping lines: one line against (1 − e^{−τ})/τ; Lyβ with O I 1025.76 against an independent evaluation (a five times finer grid); detailed balance of the zone; a blend; two lines three Doppler widths apart. With the H and O I atoms: the LTE limit with the overlap terms; the Bowen pumping of O I 3d ³D° in a dilute field; the macro-atom of that cell, whose H recombinations reach O I through the overlap channel | escape to 2×10⁻¹⁵; O I absorbs 1.2997×10⁻⁴ of the Lyβ photons (1.300×10⁻⁴), Lyβ 0.99979 of O I's (0.99980); detailed balance to 10⁻⁸; the blend escapes with 1/(τ₁+τ₂) exactly; LTE to 3.4×10⁻⁴; O I 3d/3P 2900 times higher with the overlap; 4.3 % of the H recombinations end as O I line packets with it, 0.1 % without (through the thermal pool). Without it every unit test gives the output of the code before it bit for bit |
| `test_xray` | X-ray source on a synthetic plasma table (CMFGEN's `RS_XRAY_FLUXES` format): the power of each cell against Σ_k f_k² Λ_k e^(−V_k/v) n_e,X N_atom V, with a component between two table temperatures; the emission per group: all groups together hold the power, a group over half a table bin half of that bin's share of each component, groups outside the table nothing | to 2×10⁻¹⁶ |
| `test_nontherm` | Spencer–Fano degradation of 3 keV electrons: pure H (Arnaud & Rothenflug ionization, van Regemorter Lyman series) against the heat, ionization and excitation shares of Shull & van Steenberg (1985) at x_e = 10⁻⁴–0.9; the work per ion pair of neutral H and He; grid convergence (1 and 0.5 eV bins); with the NLTE atoms (`nt_tests` of `test_nltex`, H at x_e = 0.01): the shares sum to 1, the ionization rate per atom gives the ionization share, and the Monte Carlo draws heat with the heat share | shares within 0.073; W(H I) 33.9 eV (36.3 eV, Dalgarno et al. 1999); W(He I) 45.8 eV (46 eV); heat share 0.2960 vs 0.2953; heat 0.308, sum 1 to 10⁻¹⁰, heat drawn 0.3085 |
| dielectronic resonances (`Tools/atom_recomb.py`) | recombination added by the resonances against the α columns of `die_data` for the transitions matched to model levels (×10⁻¹² cm³ s⁻¹ at 1, 2, 3×10⁴ K) | C III, 12 WI transitions: 7.95, 3.58, 2.30 against 7.94, 3.57, 2.30; its 40 AUTO transitions: 6.62, 3.50, 2.47 against 6.60, 3.50, 2.47; N IV, 11 WI: 4.12, 2.94, 2.15 against 4.07, 2.93, 2.15; `DIE=none` reproduces the atoms without them byte for byte; `test_nltex` passes with them |
| CIE table (`Tools/chianti_cie_spectra.py`, CHIANTI 11) | solar cooling function Λ_N = Σ_el (n_el/n_H) ε_el(T) from the per-element tables (Asplund et al. 2009 abundances) against Sutherland & Dopita (1993) and Gnat & Sternberg (2007) | log Λ_N = −21.41, −21.38, −21.89, −22.07, −22.54, −22.65 at log T = 5.0, 5.3, 6.0, 6.3, 7.0, 8.0 (literature −21.2 to −21.4 at the C/O peak, −21.9, −22.0, −22.6, −22.6): within 0.2 dex; C and O dominate at 10⁵ K, Fe at 10⁶–10⁷ K, free-free above (below 1 Å, ≈ 25 % of it at 10⁸ K, is outside the table) |
| cross-sections | converted valence-shell cross-sections against Verner et al. (1996) | O I at 23 eV 12.1 Mb (route A alone 3.5 Mb) against 11.9 Mb; C II, N III within 10 % to 1 keV; Fe II/III below Verner above 40 eV (CMFGEN's own data) |

## 3. Roadmap

Done in stage 2: electron temperature from the thermal balance, non-thermal
ionization/excitation of H by shock X-rays, Sobolev opacity with thermalization for the
lines of the other elements, direction-dependent Sobolev escape (§2).

Done in stage 3: items 1, 2 and 6 below (NLTE model atoms, macro-atoms and ionization
for H, He and the abundant metals), charge exchange, core-shell X-ray opacity.

**Next — physics needed for reliable line strengths across regimes**
1. ~~He I / He II macro-atoms~~ (stage 3).
2. ~~Ionization of the other elements beyond LTE~~ (stage 3).
2b. ~~Super-levels (CMFGEN's F_TO_S)~~ (stage 3); ~~low-temperature dielectronic recombination~~ (CMFGEN's `die_data` as cross-section resonances); high-temperature (Burgess) dielectronic recombination at T ≳ 5×10⁴ K for the ions without such data, and explicit autoionizing levels (CMFGEN's `auto_data`); ~~non-thermal ionization and excitation of the metals (Spencer–Fano)~~ (`in_nt_sf`).
3. ~~Time-dependent ionization~~ (`in_nltex_prev`: statistical equilibrium and release of the ionization energy over model sequences; the time-dependent thermal energy of the gas comes with the hydro L(r) of the `lumg` source).
4. ~~Non-thermal deposition of ⁵⁶Co γ-rays (Kozma & Fransson 1992)~~: `Tools/gamma_deposition.py` (Monte Carlo γ-ray transport with Compton scattering and photoabsorption, or gray absorption) gives the deposition of a homologous model for the `lumg` source, and `in_nt_lumg` routes it through the Spencer–Fano channels (`in_nt_sf`); §4.5.
5. Shock structure: adiabatic vs radiative from the cooling time; ~~cooling-function spectra (line cooling at 10⁵–10⁷ K) instead of bremsstrahlung shapes~~ (CIE spectra, `in_shock_cie`; non-equilibrium ionization behind the shock still missing); T_e ≠ T_i behind fast shocks.
6. ~~NLTE macro-atoms for Fe II, Ca II, O I, Mg II~~ (stage 3; the ε treatment of stage 2 remains for the LTE elements and the lines beyond the model atoms).
7. DDMC: the NLTE lines (Sobolev resonances exist only in IMC, so a DDMC group carries no NLTE line opacity) and thermal Compton broadening in optically thick cells. ~~Macro-atom re-emission and non-thermal deposition~~: with the NLTE macro-atoms every absorption in a DDMC group (the cell's lump of thick groups, by the rate κ S of each group) goes to the macro-atom or the thermal pool as on IMC steps, and X-ray photons make fast photoelectrons. A re-emission in a thick group of the cell keeps the packet in DDMC; any other leaves as an IMC packet from a random point of the cell, on the line it was emitted in. Before, DDMC re-emitted the absorbed energy with the continuum emissivity κ_LTE B + η of the cell, without the lines: the energy absorbed in the thick continua (the H Lyman and He continua of a wind are DDMC groups in every cell) never cascaded through the recombination lines, and the H/He benchmark lost three quarters of its Hα (§4.4). ~~NLTE source functions~~: with the NLTE macro-atoms DDMC weights its groups with the source function S = (κ_LTE B + η)/κ of the absorption opacity instead of the Planck function at T_e (η: the recombination emissivity of the NLTE levels from the Milne relation; `ddmc_spectra`). In LTE S = B. In a thick NLTE continuum the Planck weights forced J → B(T_e): the Lyman continuum of a wind, where S = B/b₁ with b₁ ≪ 1. For line work in dense CSM keep DDMC for the thick continua only: `in_trn_tauddmc` = `in_nltex_ldb` (e.g. 30), so that the optical and UV lines stay IMC (A4/A7 have only 4 cells with Δτ_es > 5, τ_es ≈ 50–65 in total) and the Lyman and He continua are DDMC groups in local detailed balance. IMC everywhere (`in_puretran = t`) with local detailed balance does not converge in the H/He benchmark and over-ionizes its wind (§4.4).

**Stage 3 — dynamics and geometry**
8. Time-dependent radiation transport (SuperNu's IMC with census) with light-travel time, for interaction phases where t_diff ≈ t_dyn.
9. 2D/3D CSM (disks, clumps, asymmetric shells) with the SuperNu 2D/3D transport.
9a. Continuum absorption of trapped line photons (Hummer & Rybicki 1985 line-photon destruction: the same trapped path as `in_nltex_vturb`, with the absorption opacity at the line, to the bound-free macro-atom or the thermal pool, and the photoionization it drives in the rates). In the H/He benchmark it matters only below 1.07×10¹⁵ cm, in the near-LTE layer above the photosphere (trapped Lyα photoionizing n = 2: β_abs/β > 1), and is below 0.2 % of β in the wind. Partial redistribution and pressure broadening of the Lorentz wings (only radiative damping now blocks the escape of scattered photons, as for Lyα).
9b. Formal integral as CMFGEN's `cmf_flux`: ~~Doppler line profiles~~ (`in_fi_doppler`), ~~electron scattering with the J of the formal solution~~ (`in_fi_nes`, Ng's extrapolation every third pass), ~~the dissolution of the high levels of H and hydrogenic ions~~ (`in_nltex_dis`, in the rates as well; Hummer & Mihalas 1988 occupation probabilities, bound-free opacity redward of the edges, as CMFGEN's `DO_DIS`). Still: an accelerated lambda operator for the scattering passes (the inner wind of a IIn, τ_es ≈ 5, converges slowly), the same scattering J in the statistical equilibrium (CMFGEN's incoherent scattering in its iterations; `in_nltex_esrates` is an experimental first step), the depth of the narrow P-Cygni trough of the H/He benchmark (0.7 of the continuum against CMFGEN's 0.25: CMFGEN's intensities on those rays fall below its own outer-wind source function, §4.4), and Stark profiles of the H and He lines (CMFGEN's `LIST_VGT`), which merge the high Balmer lines near 3646 Å. At the densities of a IIn their scale is small: the quasi-static (Holtsmark) width K_ij F₀ of H lines, K_ij = 5.5×10⁻⁵ (n_i n_j)⁴/(n_j² − n_i²) Å per esu (Underhill & Waddell 1959) and F₀ = 1.25×10⁻⁹ n_e^(2/3) esu, is 0.08, 0.9 and 8 km/s for Hα, H10 and H30 at n_e = 10¹² cm⁻³ and 0.4, 4 and 38 km/s at 10¹³ (Doppler widths 12–20 km/s), so they matter only for the series members above n ≈ 20 in the densest shells, which level dissolution already turns into continuum. On CMFGEN's own populations SuperLite's formal integral without them is within 3 % of `cmf_flux` at 3500–3800 Å relative to 4000–6500 Å (full-metal benchmark); with its own populations it is 5 % bright there in the H/He benchmark and 6–14 % in the full-metal one, which is the Balmer continuum of the populations (§4.4).
9c. Comoving-frame line transfer in slow flows: static Doppler profiles of the metal lines in the Monte Carlo (`in_nltex_dopmc`, experimental) with ~~the electron scattering of their trapped photons~~ (`in_nltex_dopes`: without it the photons scattered into the static cores were all converted, §4.4). Still: the noise of the static incident intensities (2–5 % per iteration against 0.6 % with resonances), static H and He lines (b₁ of the outer H/He wind 2.5–4 times CMFGEN's), and the outer wind of the full-metal benchmark, where CMFGEN's comoving-frame lines make H fourteen times more neutral than its Sobolev lines at 3.5×10¹⁵ cm and the hybrid four times: CMFGEN's J between Fe III 1214.56 Å and Lyα, which pumps n = 2, lies 4–13 times below the neighbouring continuum and below SuperLite's incident intensity of Lyα. CMFGEN's own incoherent scattering in its iterations makes H three times more neutral at 2×10¹⁵ cm (b₁ 13 against 4.0) and less neutral further out; the trapped-photon scattering of all lines in SuperLite's rates (`in_nltex_esrates`, four iterations from the hybrid) moves b₁ the other way (4.5 → 3.8). On CMFGEN's own populations SuperLite's J follows CMFGEN's in every band but the Lyα region (1200–1240 Å: 1.3–9 times CMFGEN's at 1.05–1.13×10¹⁵ cm), and leaving out the line-list lines beyond the model atoms (`in_nltex_bbext`) leaves the outer wind as it was (§4.4): the Lyα transport in the line forest (where CMFGEN's Fe III 1214.56 Å shadows Lyα) is what remains.

**Stage 4 — validation and production**
9d. ~~Thermal balance at sharp ionization fronts~~: the root search of the balance warm-started each trial temperature from another trial's statistical equilibrium, which at an H recombination front with the Lyman continuum in local detailed balance has a neutral and an ionized solution besides the continuing one; every trial now starts from the cell's populations (§4.4). Still: a linearization of T with the populations and the response of J, as CMFGEN's, would converge in fewer iterations.
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

### 4.3 Stage-2 runs (A7 and A4 at day 20)

Same metrics; each step adds to the previous one (10–12 iterations for the runs with
the electron-temperature solve; final runs: all stage-2 switches as written by
`Tools/stella2superlite_csm.py`):

| Run | EW5(Hα) Å | EW5(Hβ) Å | Hα absorption share | Hα peak (km/s) |
|---|---|---|---|---|
| A7 stage 1 | +17 | −55 | 0.31 | 3950 |
| A7 + Sobolev lines of the other elements, H bound-free fix | −57 | −48 | 1.00 | 5000 |
| A7 + T_e from the thermal balance, other elements at the hydro temperature | +970 | +295 | 0.00 | −100 |
| A7 + T_e, other elements in the nebular approximation (final) | +32 | −42 | 0.05 | 600 |
| A4 stage 1 | −79 | −49 | 1.00 | 300 |
| A4 all stage-2 physics, nebular approximation (final) | +241 | +9 | 0.00 | −50 |
| observed SNe IIn, 10–30 d (median) | ≈ +100 | ≈ +25 | 0 | 0 |

What the runs show:
* Sobolev lines of the other elements remove the narrow LTE metal emission (Ca II,
  Si II 6355, Fe II, the IR triplet) that the smeared opacity produced in the slow
  CSM; part of the stage-1 "Hα emission" of A7 was this contamination.  The
  transport is 4× faster.
* With the Balmer and higher continua in the opacity (stage-1 bug) and the dilute
  radiation field above the photosphere, n = 2 is overpopulated in the CSM
  (b₂ ≈ 6–10 ≈ 1/W, set by photoionization from n = 2), the Balmer continuum is
  optically thick (τ ≈ 1.5) with a source function ≈ B/b₂, and the UV below
  3646 Å is absorbed.  With the hydro temperatures Hα is in absorption.
* The electron temperature is decisive and is set by the other elements.  With
  their LTE populations at the hydro temperature, their photoabsorption of the shock
  EUV/X-rays dominates the energy exchange of the CSM and the shocked gas, T_e
  reaches 2–4×10⁴ K, and the Balmer lines are 5–10× stronger than observed with a
  flat decrement (≈ 1.7) but the observed shape: narrow core on electron-scattering
  wings.  In the nebular approximation, with their recombination emission computed
  from Saha equilibrium with the next ion at T_e (Kirchhoff's law with the actual
  populations misstates it), the outer CSM settles at 6000–12 000 K, close to the
  hydro temperatures, photoionization heating balancing recombination cooling.
* Final runs (nebular approximation): the CSM of both models is fully ionized with
  T_e ≈ 6000–9000 K in most of it (median 8900 K for A4, 8300 K for A7; T_rad of
  the dilute field 2400–11 000 K) and n = 2 strongly overpopulated (median b₂ ≈ 220
  in A4, 12 in A7).  A4 (0.2 M☉ CSM) now has a clean continuum and a strong narrow
  Hα on electron-scattering wings with no absorption — the SN IIn morphology — but
  the core is 7× higher than the observed median and Hβ is weak (Hα/Hβ ≈ 17
  against ≈ 3.6 observed).  A7 (2 M☉ CSM) has weak Hα (+32 Å) and Hβ in absorption.
  The observed median lies between the two models.
* The temperature of the shocked gas behind the cold dense shell is not converged in
  12 iterations (changes of ~20 % per iteration in a few cells, T_e up to
  4–5×10⁴ K): its cooling depends on metal ions that the nebular approximation
  describes only roughly.  The cold dense shell oscillated between iterations before
  the adaptive damping.

The line strengths of these models therefore hinge on the ionization and cooling of
the elements heavier than He in the photoionized CSM, which LTE-based treatments
bracket but cannot settle (EW(Hα) from ≈ 30 to ≈ 1000 Å for A7).  NLTE ionization
with consistent recombination and line cooling for the abundant metals (§3, items 2
and 6) is the next step; He I/He II (item 1) adds the main missing EUV opacity and
coolant of the CSM.  The steep decrement of A4 points to the Balmer line transfer in
gas with b₂ ≫ 1 (Hβ optically thicker and converted to Hα + Pα), to be checked
against an independent code (e.g. CMFGEN) on the same snapshot.

With the stage-3 physics (NLTE macro-atoms for all 43 ions of the model,
DDMC with local detailed balance for the thick continua, edge-aligned
groups, the CIE spectrum of the shock cells, T_e from the thermal balance,
the formal integral with Doppler profiles; `s9_A4`, 16 iterations and 2
more from its state) A4 has the observed decrement (R5 = 3.1–4.2 against
≈ 3.6) but EW5(Hα) = 340–390 Å, 3–4 times the observed median, and He I
5876 as strong as Hα (EW ≈ 130 Å within ±3000 km/s) with He II 4686 and
He I 7065: the shocked shell (T_e ≈ 1.8×10⁴ K, He fully He⁺) and the
CSM photoionized by the shock's CIE spectrum (He⁺ up to 6 % at 2×10¹⁵ cm,
H I 2×10⁻⁴–2×10⁻²). The solver matters here: IMC everywhere with uniform
515 km/s groups (the Lyman edge inside a group, `s8_A4`) left the same CSM
neutral (H I 0.8–0.97). There is no independent reference for this flow:
its velocity falls outward from the cold dense shell (4800 → 630 →
240 km/s), while CMFGEN's comoving-frame transfer assumes a monotonic flow.

Most of that He I emission was an artefact of the thick groups. DDMC
gives a packet absorbed in a cell's lump of thick groups to a group by
its rate κ S (`grd_spk`), and S held only the cell's own emission: at
T_e ≈ 2×10⁴ K the H and He I continua, not the He II continuum and the
EUV where the cooling flow of a shock cell emits. The shock's EUV
therefore cascaded through He I while the rate equations, whose J in
these groups is S, saw no extra ionizing radiation: He I at 10⁻³ of He
re-emitted 1.5×10⁴¹ erg/s in 5000–7000 Å (He I 5876, 6678) and
1.2×10⁴¹ erg/s at 1–2.5 µm, 17 % of the escaping luminosity in He
lines. With the emission of the shock cells in S (`s10_A4`, 6
iterations from `s9_A4_c2`) the shocked gas is 50–86 % He²⁺ between
1.52 and 1.65×10¹⁵ cm (T_e up to 3.7×10⁴ K), the edge of the cold
dense shell and the CSM out to 1.9×10¹⁵ cm are He⁺ (0.54–0.99; before
He⁰ with 0.09–0.17 He⁺ in the CSM), the He I lines vanish (8×10³⁹ erg/s
in 5000–7000 Å), He II 4686 grows, and the UV lines of C, N and O carry
4.6×10⁴¹ instead of 1.3×10⁴¹ erg/s, so that the optical luminosity
(3500–10 000 Å) falls from 1.05 to 0.88×10⁴² erg/s. EW5(Hα) is 190–200 Å
(Monte Carlo, formal integral), twice the observed median; part of the
340–390 Å before was the blueshifted He I 6678 within 5000 km/s of Hα.
He II 4686, blueshifted to −5000 km/s, fills the blue sideband of Hβ
(−16 000 to −12 000 km/s), so that the sideband continuum puts Hβ into
net absorption (EW5 −3 to −7 Å); against the red sideband alone Hβ is
67–80 Å in emission and R5 = 2.0–2.4 (`s9_A4_c2` on the same continuum:
2.7–2.8). The single-pass formal integral follows the Monte Carlo
spectrum within 2 % in every band of 3500–10 000 Å.

The formal integral of A4 against its own Monte Carlo spectrum. With
the Monte Carlo J for the electron scattering, the formal integral of
`s9_A4_c2` followed the Monte Carlo spectrum within 1–9 % in every band but
He I 5876 (0.78 of the Monte Carlo flux at 5800–5960 Å). With 12 passes
of the electron-scattering J (`in_fi_nes`) it lost flux in the bands of
the strong lines (He I 5876: 0.42, Hα: 0.64, 5500–5800 Å: 0.74,
3500–4500 Å: 0.83 of the Monte Carlo flux) while the continuum bands stayed
within 3 %; its J was 0.78–0.82 of the Monte Carlo J in every cell outside
the cold dense shell. The passes scattered photons into the cores of the
He I lines, whose source functions belonged to the inconsistent He I
populations above. With the shock in the source functions (`s10_A4`,
`fi_A4s10_c60`) the 12 passes take L 3 % below the single pass, their J is 0.95–1.05
of the Monte Carlo J in every cell (median 0.98), and the spectrum
follows the Monte Carlo one within 1–4 % in every band of 3500–10 000 Å
but 3500–4000 Å (0.92); the Balmer lines lose 10 % (Hα, F5) and 17–25 %
(Hβ, against the red sideband) of their flux to the scattering. Static
profiles of the metal lines in the slow CSM (58 of 82 cells,
`in_nltex_dopmc` 2 with `in_nltex_dopes`; `s11_A4h`, six iterations from
`s10_A4`) barely move the optical spectrum (EW5(Hα) 163–187 Å, 158 Å in the formal
integral with 12 scattering passes, `fi_A4h_q12`, within 3 % of the Monte
Carlo spectrum; with level dissolution, four more iterations after the
fixes of the MPI broadcast of the nebular ζ table, which the other ranks
had read beyond its end, and of the thermal balance (§4.4), `s12_A4`:
180 Å in the Monte Carlo, 160 Å in the 12 passes; Hβ
against the red sideband 48–54 Å, R5 2.7–2.8; He II 4686 a little
stronger and bluer), but the populations of the CSM still change by
22–29 % per iteration after six: for A4 the Sobolev lines stay the
default.

A7 with the stage-3 physics (`s10_A7`, 16 iterations from scratch, the
default reflecting inner boundary where the stage-2 A7 runs had the
thermalizing one): the CSM is fully ionized (x(H II) ≥ 0.995) at
T_e = 6000–9700 K with b₂ = 1.4–5 and b₃ = 1.7–12 (Hα thick: Sobolev depth
400–1600), and the photosphere lies in it (4.3×10¹⁵ cm). The Monte Carlo
spectrum has a narrow Hα core on weak wings, EW5(Hα) 97 Å, Hβ 47 Å,
R5 1.1 (observed median ≈ 100, 25 and 3.6): the decrement is far too
flat. It is not converged (the largest change of a state per cell,
median, 27–43 % per iteration in the last five), and its formal integral
does not follow the Monte Carlo: on the same populations (two more
iterations with 2¹⁶ packets, `s10_A7c`) the formal integral's Hα has 3–4.5 times the
Monte Carlo flux within ±500 km/s (EW5 248 against 97 Å; with 12
electron-scattering passes alike, `fi_A7c_c60`), where A4's agree within
7 %. With Sobolev resonances in the formal integral (`fi_A7c_sob`) the
core is 0.75 of the Doppler one's and the line (8.1×10⁴⁰ erg/s within
±5000 km/s) matches the Sobolev escape of the populations (net
7.5×10⁴⁰ erg/s, 5.3×10⁴⁰ of it from the CSM), while the Monte Carlo
emits 0.5–0.7 of that: the populations of the CSM are not those of the
Monte Carlo's energy flow through H. Their Hα incident intensity is the
continuum fallback (no resonance estimate in the CSM cells: the CSM
decelerates outward, so the photons of the faster gas inside reach it
blueshifted and drift further to the blue). A7's line strengths are not
reliable yet.

The cause is the electron scattering of the CSM. A7's CSM (1.9–5.3×10¹⁵ cm)
moves at 225–300 km/s and slows outward, so the velocity changes by
about one Doppler width across a cell, and electrons scatter a photon
many times while it drifts through a line: a photon that an electron
scatters into the core of Hα is absorbed there, which a Sobolev point
resonance, met only by a photon whose comoving frequency drifts through
it, does not do. Static Doppler profiles of all lines in the slow cells
absorb those photons, and the Monte Carlo and the formal integral agree
(EW5 of Hα in Å; single-pass formal integral with the Monte
Carlo J and its trapped-photon term, and with 12 scattering passes as
`cmf_flux`; the runs from `s10_A7p` on had 2¹⁶ packets per iteration,
`s10_A7y2` 2¹⁸):

| A7 day 20 | Monte Carlo | formal integral | 12 passes | R5 (MC, 12 passes) |
|---|---|---|---|---|
| Sobolev lines (`s10_A7p`) | 104 | 254 | | 0.95 |
| Sobolev lines, level dissolution (`s10_A7d`, four more iterations) | 121 | 248 | | 1.40 |
| `in_nltex_esrates` (`s10_A7y`, four more iterations) | 187 | 209 | | 1.47 |
| `in_nltex_esrates`, four more (`s10_A7y2`) | 176 | 216 | 253.6 | 1.62, 1.71 |
| static Doppler profiles of all lines in the slow cells (`in_nltex_dopmc` 2, `in_nltex_dopz` 1, with `in_nltex_dopes`; `s10_A7s`, six more iterations) | 209 | 300 | 222.5 | 1.45, 1.73 |
| the same with level dissolution and 2¹⁸ packets, eight more iterations (`s12_A7`, the recommended configuration of §5) | 236 | 290 | 232.7 | 1.47, 1.61 |
| the same after the fixes of the ζ broadcast and the thermal balance (§4.4), six more iterations (`s13_A7`) | 237 | | 226.8 | 1.45, 1.63 |

With static profiles the single-pass formal integral adds the
trapped-photon scattering of the Sobolev escape on top of the profiles'
own (300 Å); the 12 passes, which take the scattering from the J of the
formal solution, follow the Monte Carlo within 6 % (1.5 % after eight
more iterations with 2¹⁸ packets). With
`in_nltex_esrates` (Sobolev resonances that absorb the scattered photons
in their profiles) the single pass came close to the Monte Carlo, but the
12 passes give 1.44 times its Hα. A7's Hα is thus 220–235 Å with a flat
decrement (R5 1.45–1.75; the observed median is
≈ 100 Å and 3.6). The largest change of a state per cell is still
10–12 % per iteration (median over the cells; 15–17 % with 2¹⁶ packets),
which the noisy trace states dominate: EW5(Hα) of the Monte Carlo stays at
227–263 Å over the last seven iterations of `s12_A7` (236 ± 12 Å), as
A4's at 180–189 Å over those of `s12_A4`. For dense, slow CSM the static profiles
of all lines are needed; the default Sobolev lines are consistent where
the flow is fast (A4's formal integral follows its Monte Carlo within 3 %).

### 4.4 Benchmark against CMFGEN: IIn-like dense wind

CMFGEN models computed for this comparison (CMFGEN 2023 with the 21jun23
atomic data), all on one structure: R* = 10¹⁵ cm, L = 2.4×10⁴² erg/s,
Ṁ = 0.02 M☉/yr, v∞ = 250 km/s (β = 1), VTURB = 20 km/s. Electron
scattering is coherent in the comoving frame during CMFGEN's iterations
(`COH_ES`) and incoherent in its observer-frame spectrum (`cmf_flux`).
SuperLite runs on the same structure with T_e fixed to CMFGEN's, the same
model atoms (CMFGEN super-levels) and groups of 100 km/s at 3000–10 000 Å
(1000 km/s elsewhere). CMFGEN's observer spectra are in vacuum wavelengths
(c/ν) and SuperLite's lines at air wavelengths above 2000 Å: the CMFGEN
spectra are converted to air before profiles are compared (83 km/s at Hα;
the comparisons before this was done had CMFGEN's profiles 83 km/s to the
red, which moved part of its line core into its red wing).

Summary of the full-metal benchmark (details below). With the same atomic
data (`in_nltex_bbext = f`: no line-list lines beyond the model atoms) and
CMFGEN's own populations (`in_nltex_hold`), SuperLite's transport gives
CMFGEN's emergent spectrum within 5 % from 1500 Å to 2.5 µm and its J
within 4 % in the Lyman and He continua; its formal integral makes Hα, Hβ
and Pβ 6, 12 and 3 % stronger than `cmf_flux`. With its own populations,
like for like, Hα, Hβ and Pβ are 17, 30 and 5 % stronger, the Balmer
continuum 10–14 % brighter relative to the optical (b₂ 2–17 % low where
it forms), and H in the outer wind 2–4 times more ionized. Most of that
Balmer continuum is energy CMFGEN's model does not re-emit: at the fixed
T_e of its Sobolev model its gas absorbs 1.6×10⁴¹ erg/s (6.6 % of L) more
than it emits where the Balmer continuum forms, while SuperLite's Monte
Carlo re-emits what the gas absorbs; against CMFGEN's radiative-equilibrium
Sobolev model SuperLite's Balmer continuum agrees within 3 % (below).
Open: the Lyα transport in the slow line forest and the H ionization of
the outer wind (roadmap 9c).

Summary of the H/He benchmark. On CMFGEN's own populations SuperLite's
formal integral (12 scattering passes) gives `cmf_flux` within 1–4 % in
the continuum and Hα, Hβ, Pβ, He I 10830 within 2.5 %. With its own
populations (Sobolev lines in the iterations, as CMFGEN's Sobolev
option; Doppler profiles and 12 scattering passes in the formal
integral) Hα, Hβ and Pβ are 105.9, 22.6 and 55.0 Å against `ies`'
100.9, 21.1 and 55.2 (+5, +7, 0 %), with level dissolution 108.0, 23.0
and 54.2. Its T_e from the thermal balance stays within 0.5 % of
CMFGEN's radiative-equilibrium T_e in 64 of 69 cells and within 1.8 %
in all, the H recombination front included. Open: the incoherent
scattering of CMFGEN's iterations, which ionizes H just above the
photosphere (`in_nltex_esrates`, experimental, does part of it).

| CMFGEN model | lines | T_e | state |
|---|---|---|---|
| m002: H, He, C, N, O, Ca, Fe | comoving frame | held at its starting structure (≈ 5000 K beyond 2×10¹⁵ cm) | inner wind converged; beyond ≈ 2×10¹⁵ cm the H I ground state changes by factors of several per iteration |
| H/He | Sobolev | the same | outer wind oscillates by ≈ 10 %, ionized |
| H/He | Sobolev | radiative equilibrium (released after 25 iterations; ≈ 7400 K in the outer wind) | converged (< 10⁻⁴) |
| H/He | comoving frame | that T_e, fixed | converged |
| H/He | comoving frame, incoherent electron scattering in the iterations | that T_e, fixed | converged |
| m002 (`sn_iin_m002_re`, `_blk`, `_ies`) | Sobolev; comoving frame; comoving frame with incoherent electron scattering in the iterations | radiative equilibrium of the Sobolev model, then fixed | converged; the last: 4.8 % in the outer wind after 70 iterations |

* m002's 5000 K outer wind is CMFGEN's starting guess, never relaxed. At
  that temperature n = 2 is photoionized by the Balmer continuum ≈ 300
  times faster than it decays, so net recombination runs only through
  two-photon decay, Lyα escape and the escape of Lyman-continuum photons,
  with an effective recombination time (≈ 10 yr) longer than the flow
  time (≈ 8 yr): an ill-conditioned state that neither code holds well.
  The H/He models below are the reference.
* Sobolev against comoving-frame line transfer in CMFGEN's iterations,
  H/He at the same radiative-equilibrium temperature: populations within
  1–5 %, Hα EW 96.2 against 96.7 Å, Hβ 21.1 against 20.9 Å, Pβ 49
  against 51 Å. Incoherent electron scattering in CMFGEN's iterations:
  Hα 101, Hβ 21, Pβ 55 Å. The observer spectra of all three are
  comoving-frame transfer (`cmf_flux`) with the Doppler profiles of the
  lines, so they include the electron scattering of photons trapped in
  the resonance zones, which Sobolev point resonances lack (below).

The H/He wind at the radiative-equilibrium T_e (SuperLite: DDMC with
`in_trn_tauddmc` = 5 and local detailed balance from τ = 30, from the
same converged state; EW of the Monte Carlo spectrum of the last two
iterations, and of the formal integral):

| | H I fraction at 1.1 / 1.7 / 3.5 / 9.5×10¹⁵ cm | Hα | Hβ | Pβ (EW, Å) |
|---|---|---|---|---|
| CMFGEN (Sobolev) | 0.038 / 4.0×10⁻⁴ / 0.0047 / 0.095 | 96 | 21 | 49 |
| DDMC re-emitting with the continuum emissivity (before) | 0.12 / 5.6×10⁻⁴ / 0.0034 / 0.036 | 27 | 0.7 | 8 |
| DDMC through the macro-atom | 0.098 / 2.9×10⁻⁴ / 0.0024 / 0.031 | 101 (formal 137) | 21 (29) | 29 (38) |
| … and the continuum as the fallback incident intensity of the lines | 0.11 / 3.5×10⁻⁴ / 0.0029 / 0.038 | 82 (formal 81) | 16 (20) | 48 (43) |
| IMC everywhere (`in_puretran`), before the fallback | 1.6×10⁻⁴ / 5.4×10⁻⁵ / 2.2×10⁻⁴ / 0.0045 | 84 | 16 | 19 |

* Case A in the Monte Carlo Lyman continuum. Its estimate in a cool,
  thick wind is mostly zero, as the continuum carries a small share of the
  cells' emission, and the on-the-spot cycle turned into net
  recombination: m002's outer wind recombined (H I fraction 0.5–0.9 against
  CMFGEN's 0.01–0.2) in every transport mode. Local detailed balance
  (`in_nltex_ldb`) removes it.
* DDMC lost the recombination lines. The Lyman and He continua of the
  wind are DDMC groups in every cell, and DDMC re-emitted the energy
  absorbed there with the continuum emissivity of the cell, never through
  the macro-atom. In a thick Lyman continuum ≈ 60 % of the activations of
  the H continuum state recombine to the ground state, so most of the
  energy that photoionizes the wind circulates there; with the continuum
  emissivity it left as Balmer and Paschen continuum instead of cascading
  through the lines (case B). Every DDMC absorption now goes through the
  macro-atom or the thermal pool (roadmap 7): Hα rose from 27 to 101 Å,
  Hβ from 0.7 to 21 Å, with 4–5 times more DDMC steps (still a third of
  the IMC steps; the absorbing group in a lump is drawn by bisection).
* The lines' incident intensity. Hα has fewer than 30 resonance crossings
  per cell in the whole wind (≈ 9 photons per cell and iteration within
  25 km/s of it), so the rates and the macro-atom used the J of the line's
  1000 km/s group. In the outer wind that group carries 4–5 times more Hα
  photons than continuum, and those only redshift away from their line:
  the pumping by its own escaping photons raised n₃ (and S_Hα) 2–2.7
  times above CMFGEN's, and the formal integral with these populations
  gave Hα 137 Å against 101 Å of the Monte Carlo spectrum. The fallback
  is now the group J without the photons of its lines that have kept
  their frequency history (`nxt_tally_finish`), the continuum incident
  intensity of CMFGEN's Sobolev option. The Sobolev Hα emission of the
  populations dropped from 3.6 to 2.1×10⁴⁰ erg/s (CMFGEN 1.9×10⁴⁰), and
  the formal integral and the Monte Carlo spectrum now agree (Hα 81 and
  82 Å).
* Remaining in the populations: the outer wind 1.6–2.5 times more
  ionized; He I 5876/10830 stronger. IMC everywhere with local
  detailed balance does not converge in this wind (8–17 % changes per
  iteration) and drifts to an over-ionized state, also with the balance
  extended to τ ≥ 1: there J_ext adds 13–25 % to the Lyman-continuum J
  of the outer wind, where He is neutral; its source is not yet
  identified. Use DDMC for the thick continua (roadmap 7).
* m002 with all elements (H, He, C, N, O, Ca, Fe; CMFGEN rerun to
  convergence with Sobolev lines and T_e in radiative equilibrium,
  `sn_iin_m002_re`, 6544 lines against SuperLite's 6663 from the same
  atoms and cut; SuperLite with T_e fixed to it): the emergent spectrum
  agrees within 2–14 % in every band from 1250 Å to 1 μm (912–1250 Å:
  1.67 times CMFGEN's, 0.8 % of L), the continuum of 3500–13 500 Å within
  0–5 %, and b₂, b₃ of H within 10–20 % inside 3×10¹⁵ cm. The formal
  integral with 12 passes of the scattering J: Hα 102.8 against 92.4 Å,
  Hβ 22.8 against 19.8, Pβ 52.3 against 47.4, He I 10830 −3.4 against
  −3.5, He I 6678 −1.7 against −1.4 Å; the electron-scattering wings
  within 0.45 Å per bin of 150–1300 km/s (Hα; Hβ 0.35, Pβ 0.2 Å), the
  continuum 1.5–3.4 % above CMFGEN's at 5000–13 300 Å, and the excess in
  the line cores (Hα +5.5 Å at −150…+300 km/s), where SuperLite's
  outer-wind n₃ is 1.1–1.5 times CMFGEN's Sobolev model's, as in the H/He
  wind, whose outer populations agree with CMFGEN's comoving-frame `ies`
  model instead. (With the Monte Carlo J the wings and Hβ, Pβ lie in its
  noise: Hα 89, Hβ 16, Pβ 59 Å.) Inside the wind, SuperLite's J
  at 912–1240 Å and in the Lyman continuum is 1.6–4 times CMFGEN's
  (CMFGEN's J at depth from its `EDDFACTOR` records) at 1.1–2×10¹⁵ cm,
  while longward of 1240 Å it agrees within 10–18 %. Fe II has
  photoionization edges from its excited levels at 900–1600 Å, so Fe II
  (and H I and O I, tied to it by charge exchange and the Lyman
  continuum) is 2–3 times less abundant there (b₁ 2–5 times lower). (The
  narrow absorptions of Fe II 4924, 5018, 5169 Å, smeared over the
  100 km/s groups of the Monte Carlo spectrum (EW 0.6–1.2 Å), are in the
  formal integral with Doppler profiles: EW 1.6–2.0 Å against 1.8–2.2 Å in
  CMFGEN's Sobolev model and 1.6–2.1 Å in its comoving-frame model,
  troughs to 0.11 of the continuum against 0.05; O I 7774 1.4 Å against
  1.8 and 3.0 Å.) The H/He wind shows no such
  excess (912–1230 Å within 25 %), so it is the metal lines. Neither DDMC
  (`in_trn_tauddmc` 5 → 30: no NLTE lines in DDMC groups, roadmap 7) nor
  incoherent scattering (coherent as CMFGEN's `COH_ES`) changes it.
  CMFGEN's Sobolev option leaves the lines out of its continuum J (no
  blanketing, and no line photons in the J of the photoionizations),
  whereas SuperLite's J holds both, so the reference itself was in
  question there. CMFGEN's comoving-frame model at the same T_e
  (`sn_iin_m002_blk`: `GLOBAL_LINE = BLANK`, converged) settles it the
  other way: the blanketing of its line forest takes the far-UV J further
  down. Against it SuperLite's J (`in_trn_tauddmc` 30) is 6–16 times
  higher at 912–1100 Å and 2.3–3.5 times at 1100–1200 Å at
  1.1–6×10¹⁵ cm; in the Lyman continuum 1.2–3.8 times inside 2×10¹⁵ cm
  and 15–34 times at 2.5–6×10¹⁵ cm; 1.1–1.7 times at 1230–3000 Å; within
  13 % at 3700–5000 Å. With the comoving-frame lines CMFGEN's b₁ of H
  rises over its Sobolev model's (×1.9, 6.7, 17 at 2, 2.5, 3.6×10¹⁵ cm;
  b₂–b₄ within 15 %), and H, O and, by charge exchange, Fe recombine in
  the outer wind: H I 3.7×10⁻³, 7.2×10⁻², 0.17 of H at 2, 3.5,
  6×10¹⁵ cm against SuperLite's 8.5×10⁻⁴, 3.1×10⁻³, 1.0×10⁻²; Fe II
  0.19, 0.83, 0.93 of Fe against 0.04, 0.15, 0.42. The Sobolev Monte
  Carlo misses most of the blanketing of this slow wind (250 km/s,
  coasting outside 2×10¹⁵ cm): along its whole path through the coasting
  wind a photon's comoving frequency drifts by v(1 − μ) at most (μ where
  it starts; not at all for radial photons), so it meets only the
  resonances just to its red, while CMFGEN's Doppler profiles let
  overlapping lines absorb and fluoresce the line photons in place.
  Static Doppler profiles of the metal lines in the slow cells
  (`in_nltex_dopmc` 2, `in_nltex_dopz` 3; six iterations from the
  Sobolev state, microturbulence 10–20 km/s) take SuperLite's J to
  1.0–2.9 times CMFGEN's at 912–1100 Å, 0.5–0.7 times at 1100–1500 Å and
  0.8–1.3 times at 1500–3000 Å, the Lyman continuum beyond 2.5×10¹⁵ cm to
  7–16 times, b₁ there up by 2.2 and Fe II to 0.28 and 0.58 of Fe at 3.5
  and 6×10¹⁵ cm. With CMFGEN's iteration microturbulence (VTURB 20 km/s
  everywhere: wider profiles, more overlap; eight more iterations) H I
  is 3.0×10⁻³, 1.1×10⁻², 5.2×10⁻², 0.19 of H at 2, 3.5, 6, 9.5×10¹⁵ cm
  against CMFGEN's 3.7×10⁻³, 7.2×10⁻², 0.17, 0.19 (Sobolev:
  8.5×10⁻⁴…4.1×10⁻²), Fe II 0.15, 0.39, 0.77, 0.95 against 0.19, 0.83,
  0.93, 0.95, b₁ 4.8, 14, 72, 836 at 2, 2.5, 3.6, 5.9×10¹⁵ cm against
  4.0, 59, 666, 3033; the J at 912–1100 Å is 0.4–1.0 times CMFGEN's,
  but 1100–1500 Å falls to 0.2–0.6 times and the Lyα band to 0.1–0.9.
  Against CMFGEN's comoving-frame spectrum (`cmf_flux`, VTURB 10–20 km/s;
  EW within ±2700 km/s against the continuum at 4000–6000 km/s) Hα is
  82.9 Å (10–20 km/s) and 71.8 Å (20 km/s, formal integral at 10–20 km/s)
  against 80.4 Å (Sobolev: 98.1; CMFGEN's Sobolev model 91.5), Hβ 13.6
  and 17.3 against 17.4, Pβ 51.0 and 53.2 against 45.7; the narrow Fe II
  absorptions 1.2–1.4 Å, O I 7774 1.9 Å (CMFGEN 3.0 Å). The iterations
  with static profiles still change the populations by 4.5–6 % per
  iteration after 14 (largest change per cell, median; 0.6 % with
  resonances). Four more, counting the window crossings instead of the
  short flights of the slow cells for the incident-intensity threshold
  (`cmf_m002re_dm7`), change them by 4.8, 4.6 and 3.5 % and carry the
  outer wind further towards CMFGEN (b₁ 6.1, 18, 90, 883; H I 3.8×10⁻³,
  1.4×10⁻², 5.4×10⁻², 0.24; Fe II 0.19, 0.44, 0.78, 0.96), while Hα
  falls to 63.3 Å (core 45 Å against 53, wings 7 and 11 Å against 10
  and 18; Hβ 18.1, Pβ 52.8 Å): b₃ of the inner wind stays at 0.75–0.8 of
  CMFGEN's (1.23 against 1.62 at 1.1×10¹⁵ cm, 4.9 against 6.1 at
  1.65×10¹⁵ cm). The hybrid is not converged and stays experimental.
  Its far-UV deficit is not a transport inconsistent with its
  populations: a formal solution of 1100–2000 Å with 12 electron-
  scattering passes on the populations of `cmf_m002re_dm7`
  (`cmf_m002re_dm7_fiuv`) gives the emergent luminosity there within 5 %
  of the Monte Carlo one and a J 5–20 % above the Monte Carlo J in every
  cell of the wind (the trapped photons of the static cores). It comes
  with the electron scattering. The hybrid's emergent luminosity at
  1100–2000, 2000–3000 and 5000–10 000 Å is 0.52, 0.70 and 1.27–1.29 of
  CMFGEN's (Sobolev transport: 1.10, 1.13, 1.01–1.03): energy moves
  from the UV to the optical continuum, whose excess lowers the EW of Hα.
  With coherent scattering in the comoving frame (`in_thm_thermal = f`,
  as CMFGEN's iterations, `COH_ES`; `cmf_m002re_dm8c`, four iterations
  from `cmf_m002re_dm7`) the same populations give 0.87, 0.87 and
  1.15–1.17 in the first iteration and 0.92, 0.91 and 1.13–1.14 after three
  more, with populations changing by 6 and 1.5 % in the last two
  iterations instead of 3.5–4.8 %; b₃ of the inner wind rises to 1.32
  (1.1×10¹⁵ cm). Scattering off thermal electrons (±550 km/s at 10⁴ K)
  moves the photons between
  the lines into the cores of the static metal-line profiles, where the
  macro-atom converts them; coherent scattering leaves them between the
  lines. The macro-atom of a static line escaped by the Doppler escape β
  alone, not by the electron scattering of its trapped photons, which
  takes part of them out of the core unconverted. With that escape and
  the photons scattered into the core in the rates of the static lines
  (`in_nltex_dopes`; `cmf_m002re_dm9`, four iterations from
  `cmf_m002re_dm7` with incoherent scattering) the bands are 0.98, 0.91
  and 1.11–1.13 of CMFGEN's (whose emergent luminosity is 6 % below its
  inner-boundary luminosity at the fixed T_e, so that SuperLite's bands
  lie 6 % high on average), b₃ of the inner wind is 1.32, and the formal
  integral at 10–20 km/s gives Hα 76.9 Å against 80.4 Å (core 53.3
  against 52.6 Å, wings 9.2 and 14.4 against 10.2 and 17.6 Å), Hβ 12.4
  against 17.4 Å and Pβ 43.7 against 45.7 Å, with the populations still
  changing by 4.6 % per iteration. Six more iterations
  (`cmf_m002re_dm11`; 3.8, 4.3, 3.5, 2.7, 2.1 % per iteration) give the
  bands 1.03, 0.96 and 1.07–1.09 of CMFGEN's (0.97, 0.90 and 1.00–1.02 of
  its shape) and b₃ = 1.39 in the inner wind. Like for like, with 12
  electron-scattering passes in the formal integral (`in_fi_nes`, as
  CMFGEN's 10 incoherent passes of `cmf_flux`), the hybrid gives Hα
  96.3 Å (core 60.4, wings 13.2 and 22.7 Å), Hβ 19.4, Pβ 52.9 Å, and the
  Sobolev transport 104.8, 20.4, 52.9 Å, against CMFGEN's 80.4, 17.4,
  45.7 Å (its Sobolev model: 91.5, 19.5, 48.1 Å). CMFGEN's comoving-frame
  model with incoherent scattering in its iterations as well
  (`sn_iin_m002_ies`: `COH_ES = F` from the converged model, 70 more
  iterations, the last change 4.8 % in the outer wind; `cmf_flux` as
  before) loses 3, 11 and 7 % of its luminosity at 1100–2000, 2000–3000
  and 3000–3650 Å (the optical within 0.6 %), and gives Hα 85.8, Hβ 17.0,
  Pβ 51.1 Å; its H I is 3 times higher at 2×10¹⁵ cm and 0.4–0.6 times at
  3.5–9.5×10¹⁵ cm (b₁ 13, 85, 365, 1007 at 2, 2.5, 3.6, 5.9×10¹⁵ cm),
  Fe II 0.43 and 0.74 of Fe at 2 and 3.5×10¹⁵ cm. So the scattering
  shifts CMFGEN's UV by a tenth, not by half as the hybrid without the
  trapped-photon escape did. Against this model, the physically complete
  one, the hybrid's lines are 12, 14 and 4 % stronger (Sobolev
  transport: 22, 20, 4 %), and its spectral shape (normalized to its 9 %
  higher total: CMFGEN's emergent luminosity is below its inner-boundary
  one at the fixed T_e) agrees within 4 % from 1100 Å to 2.5 µm except
  in 3000–3650 Å (15 % high). Static
  H and He lines as well (`in_nltex_dopz` 1, six more iterations,
  `cmf_m002re_dm10s`) change the populations by 1.6–2.5 % per iteration
  and do not bring the outer wind closer (b₁ 14, 61, 685 at 2.5, 3.6,
  5.9×10¹⁵ cm). The outer wind does not follow the
  scattering: H I is
  2.9×10⁻³, 9.5×10⁻³, 3.9×10⁻² at 2, 3.5, 6×10¹⁵ cm with coherent
  scattering (b₁ 4.5, 13, 61, 585). There CMFGEN's comoving-frame J
  between
  Fe III 1214.56 Å and Lyα is 4–13 times below the continuum next to it
  from 1.1×10¹⁵ cm outward (1.1–1.7×10⁻⁸ erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹ at
  1.65–3.6×10¹⁵ cm), while SuperLite's incident intensity of Lyα is
  2–7 times that (Sobolev transport: 7–19 times): Lyα pumps n = 2, and
  the Balmer continuum ionizes it, less in CMFGEN (b₂/b₁ = 0.12 at
  3.6×10¹⁵ cm against 2.0 in its Sobolev model and 1.0 in the hybrid).
  With Sobolev lines both codes agree there (H I 3.1×10⁻³ and 1.0×10⁻²
  against CMFGEN's 5.1×10⁻³ and 1.4×10⁻² at 3.5 and 6×10¹⁵ cm); the
  comoving-frame lines raise CMFGEN's H I fourteenfold, the hybrid's
  fourfold, and O follows H by charge exchange in both. Without the
  trapped-photon escape, static profiles of all lines, H and He
  included, gave nearly the same far-UV J but took Hβ to 8 Å and Pβ to
  30 Å (the H/He benchmark, below). With DDMC only from τ = 30 the iterations are 2.5 times faster (fewer
  DDMC–IMC swaps) with the same populations.
* What sets the H ionization of the outer full-metal wind (roadmap 9c).
  Two pieces were missing. (1) Rayleigh scattering by H I in the Lorentz
  wing of Lyα (CMFGEN's `INC_RAY`): 1732 σ_T per H I atom within
  ≈ 1500 km/s of Lyα, a scattering depth of 25–40 above the photosphere
  (electrons: 0.1–3.5). On CMFGEN's own populations (`rayE_hold`,
  `in_rayleigh = 'lyman'`) SuperLite's J at 1210–1220 Å goes from 7.0,
  5.0, 2.4 and 1.5 times CMFGEN's to 1.6, 1.2, 1.05 and 0.86 times at
  1.10, 1.13, 1.65 and 3.6×10¹⁵ cm. The dip of CMFGEN's J between Fe III
  1214.56 Å and Lyα was this scattering, not the Fe III line. With its own
  populations (`rayF_nb`, six iterations) b₁ rises from 1.45 to 2.0 at
  1.65×10¹⁵ cm (CMFGEN 2.7), and the outer wind does not move. In the H/He
  wind (`rayG_hhe`) b₁ is within 1.05–2.6 times `ies`'. Redistributing the
  scattered photons with the electrons as CMFGEN's incoherent scattering
  does (`'lymes'`) brings that to 5–33 % (`rayG2_hhe_es`), and changes the
  full-metal wind by ≈ 10 % (`rayF2_nb_es`).
  (2) Overlapping lines. The rates of the outer wind (`SL_RATEDUMP`,
  against CMFGEN's `HIPRRR`) feed SuperLite's ground state of H with Lyα
  and drain it by charge exchange with O II (6.5×10⁻⁸ s⁻¹ per H atom at
  3.4×10¹⁵ cm), whose O I is photoionized from excited levels.
  CMFGEN's own models disagree with each other there: with Sobolev lines
  (`re`) b₁ = 38 at 3.6×10¹⁵ cm, half of SuperLite's, and with comoving-
  frame lines 663 (`blk`) and 364 (`ies`). In the H/He wind the three
  agree within 15 %, so the effect comes with the metals. In CMFGEN's
  comoving-frame models the Lyman lines feed the ground state 2.5–6 times
  more per H atom than in its Sobolev model, with 10–17 times more H I:
  25–100 times the Sobolev escape. And O I 3d ³D° is 10–100 times more
  populated, relative to the O I ground level, than in SuperLite (5–15
  times more than in CMFGEN's Sobolev model). O I 1025.76 Å lies 11.7 km/s
  to the red of Lyβ, inside its Doppler core, and absorbs the trapped Lyβ
  photons: the Bowen fluorescence of O I 8446 and 11 287 Å, which CMFGEN's
  spectrum shows in emission (EW 0.3 and 1.4 Å; SuperLite without
  the overlap −0.9 and −0.9 Å). With χ(O I)/χ(Lyβ) ≈ 1.3×10⁻⁴, O I takes a trapped Lyβ
  photon long before it escapes (β ≈ 1–4×10⁻⁷), and n = 3 empties into
  the ground state, with SuperLite's n₃ ≈ 10⁻⁶ s⁻¹ per H atom, ten times
  its Lyα flow. Point resonances keep each line's trapped photons to
  itself; `in_nltex_overlap` (§2) shares the zone.
* The far-UV radiation field of the IIn wind (925–1200 Å). It pumps the
  higher Lyman lines (Lyβ–Lyζ are in net absorption in the outer wind of
  both codes, from the continuum on their blue side, 16–1000 times their
  source functions), so it enters the ground-state balance of H.
  SuperLite's Monte Carlo J there is 0.68–0.90 times that of CMFGEN's
  model with incoherent electron scattering (`ies`) at 925–1020 Å and
  0.63–1.0 times at 1030–1190 Å (H/He wind; 0.55–0.86 times `blk`'s,
  whose coherent scattering keeps the far-UV photons out of the Lyman
  lines: with thermal electrons they random-walk in frequency and land in
  the lines). The innermost cell has J/B(T) = 0.7 at 950–1300 Å against
  CMFGEN's 1.02–1.07. Ruled out, each with one iteration on the same
  populations (`bcG_*`): the spectrum of the injected luminosity
  (`in_srcdiff`, the diffusion flux: ≤ 4 % in any cell), a thermalizing
  inner boundary (`in_innerbc = 'thrm'`: 2–10 %), and the lumping of the
  thick DDMC groups (`in_taulump` = 10³⁰: no change). The CMFGEN
  `EDDFACTOR` J of the full-metal `ies` is not consistent with its own
  populations to a few 10⁻³: every Lyman line and the Lyman continuum
  have the same J/S − 1 at a given radius (+0.003 to +0.008 at
  1.4–2.1×10¹⁵ cm), which a thick line cannot have. `blk` is consistent
  to ≤ 2×10⁻⁴: Lyα–Lyδ within 10⁻⁴ of S, the Lyman continuum +1×10⁻⁴
  at its edge. (An earlier figure of +0.012 to +0.017 came from an
  analysis script that applied the reduced-mass correction twice.)
  CMFGEN's own rates (`HIPRRR`) have P₁ = R₁ in the H/He wind (in both
  models), as SuperLite's local detailed balance. In the full-metal wind
  P₁ − R₁ = 1–2×10⁻⁶ s⁻¹ per H particle (both models): the metal lines
  inside the Lyman continuum (next item).
* The incident intensity of the Lyman lines and the ground state of H.
  With the β-weighted estimators (`in_nltex_jbbeta`, `rayGk_hhe`) b₁ of
  the H/He wind is 2.1 times `sob`'s at 1.6×10¹⁵ cm and 1.6–2.1 times
  beyond. Before the weighting (`rayGh_hhe`) it was 1.5 and 0.67–0.98
  times.
  The same binary without the weighting (`rayGq_hhe`, six iterations from
  `rayG_hhe`), b₁ at 2.5, 3.6, 5.9, 9.5 and 15×10¹⁵ cm:

  | | 2.5 | 3.6 | 5.9 | 9.5 | 15 |
  |---|---|---|---|---|---|
  | `rayGq_hhe` (no weighting) | 4.42 | 42.5 | 524 | 3540 | 4380 |
  | `rayGk_hhe` (weighting) | 9.05 | 68.2 | 905 | 10 600 | 13 600 |
  | CMFGEN `sob` | 3.40 | 48.5 | 560 | 5150 | 7090 |
  | CMFGEN `blk` | 3.53 | 49.2 | 489 | 2620 | 2930 |
  | CMFGEN `ies` | 3.78 | 38.8 | 308 | 1600 | 1770 |

  | `rayGkl_hhe` (no weighting, `in_rayleigh = 'lymes'`) | 3.84 | 38.0 | 403 | 3420 | 4400 |

  So the weighting is off by default (`in_nltex_jbbeta`). Redistributing
  the Rayleigh-scattered Lyα-wing photons with the thermal electrons as
  `ies` does (`'lymes'`, six iterations from `rayG_hhe`) brings b₁ to
  `ies`' within 2 % at 2.5–3.6×10¹⁵ cm and to 1.3, 2.1 and 2.5 times it at
  5.9, 9.5 and 15×10¹⁵ cm (without it 1.2–2.8 times). `'lymes'` copies
  CMFGEN's approximation (Rayleigh scattering is coherent in the atom's
  frame), so it stays an option for like-for-like comparisons with
  CMFGEN's incoherent-scattering models. The single-cell harness (`Tools/cell_jtest.f90`, J sampled 100
  km/s to the blue of each line) takes the rates out of the question: the
  difference is in the incident intensities of the H lines, and which line
  matters depends on CMFGEN's electron scattering. x(H I) of cell 47
  (1.9×10¹⁵ cm):

  | J in the harness | x(H I), cell 47 |
  |---|---|
  | SuperLite's own | 1.33×10⁻³ |
  | `sob` J, H lines / below 980 Å only | 3.9×10⁻⁴ / 4.5×10⁻⁴ |
  | `blk` J (coherent e.s.), H lines | 4.3×10⁻⁴ |
  | `ies` J (incoherent e.s.), H lines / Lyα and up / Lyβ and up / below 980 Å | 6.3×10⁻⁴ / 6.8×10⁻⁴ / 1.35×10⁻³ / 1.7×10⁻³ |
  | CMFGEN's own result: `sob`, `blk`, `ies` | 5.7, 5.9, 8.1 ×10⁻⁴ |

  With coherent scattering the continuum just blueward of the higher
  Lyman lines stays at its level (8×10⁻⁸ erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹), where
  SuperLite's estimator gets 3.8×10⁻⁸ at Lyγ (unweighted; 2.6×10⁻⁸
  weighted). With thermal electrons, which SuperLite has, the higher
  Lyman lines see less than SuperLite's estimators give. What is missing
  then is Lyα: `ies` has J ≈ 2×10⁻⁷ at 100–5000 km/s on either side of
  it, 0.73 of its source function, and no trough. `sob`, `blk` and
  SuperLite have a trough over that range (1–7×10⁻⁸). SuperLite's Lyα is
  fed at 0.15 S (0.23 with `in_nltex_esrates`, `rayGel_hhe`), so its net
  decay into the ground state is 3 times `ies`'s. The trough comes from
  the Rayleigh wing of Lyα (scattering depth 25–40 above the photosphere,
  item above). It forms between 1.08 and 1.11×10¹⁵ cm: inside that,
  SuperLite's J near Lyα is 0.7–1.1 times `ies`'s; at 1.1–1.3×10¹⁵ cm,
  0.1–0.25 times at 1000–3000 km/s from the line. `ies` has no trough
  because CMFGEN's incoherent scattering redistributes its whole scattering
  opacity with the thermal electrons, Rayleigh included (`comp_j_blank`:
  the source of the e.s. term is ESEC_EXT × RJ_ES with ESEC_EXT built from
  CHI_SCAT = ESEC + CHI_RAY). Scattering by an H atom is coherent to its
  thermal speed (10 km/s), and the coherent models (`blk`, `sob`) have the
  trough as SuperLite has. `'lymes'` reproduces CMFGEN's treatment for the
  like-for-like comparison: before the β weighting it brought b₁ to 5–33 %
  of `ies`'s (`rayG2_hhe_es`); with it, `rayGkl_hhe`. The trapped-photon
  scattering (Hummer & Rybicki's P(r), r = κ_es/κ₀ = 1.3×10⁻⁸ in cell 47)
  sends 40 % of the escaping Lyα photons out of the core, but in the
  Monte Carlo it acts only with `in_nltex_esrates` (`nx_linees`).
* The H ionization of the inner full-metal wind (1.2–2.5×10¹⁵ cm),
  where SuperLite's x(H I) was 3–16 times CMFGEN's (`blk`) after the
  overlap runs (`rayFo1b`, `rayFo1c`). `Tools/cell_jtest.f90` solves
  SuperLite's statistical equilibrium of those cells in CMFGEN's
  radiation field (`Tools/cmfgen_j_cells.py`: CMFGEN's comoving J
  averaged over the groups, the incident intensity of each Sobolev line
  CMFGEN's J five Doppler widths to its blue, every group with Lyman-
  continuum opacity in detailed balance with the current populations).
  It gives x(H I) 6–9 times CMFGEN's at 1.2–1.9×10¹⁵ cm (own J: 7–14
  times), with n = 2–30 within 4 % of CMFGEN's: only the ground state is
  off. (The first version of the tool fed CMFGEN's J to the group that
  straddles the Lyman edge, 911.44–911.84 Å. That J, CMFGEN's S_LyC,
  pinned n₁ to CMFGEN's and gave a false agreement within 1.2–2.0 times.)
  The ground state of thick Lyman lines and continuum is set by small net
  flows. In SuperLite (cell at 1.58×10¹⁵ cm) Lyβ photons absorbed by
  O I 1025.76 (Bowen) feed it, 5×10⁻⁷ s⁻¹ per H particle. Charge
  exchange O⁺ + H → O + H⁺ after O I photoionization drains it at the
  same rate. With both switched off x(H I) stays ~10 times too high.
  CMFGEN drains it by a net photoionization in the Lyman continuum:
  740 cm⁻³ s⁻¹ at 1.65×10¹⁵ cm, 0.75 % of the gross rate (`HIPRRR`; its
  O I and charge exchange are in balance, `OIPRRR`). In CMFGEN's J the
  net comes from J/S_LyC = 1.1–1.6 at the Fe III 3d⁶ ⁵D–3d⁵4p lines
  (808–866 Å), C II 858 and O II 832, while the continuum between them
  is within 10⁻⁴ of S_LyC. These are line photons emitted inside the
  Lyman continuum, which H I absorbs. SuperLite's lines there emit about
  as much (Fe III 325, C II 330, O II 51 cm⁻³ s⁻¹). With DDMC, however,
  the rates of the thick Lyman continuum took J = S of the continuum
  alone and never saw these photons (net 38 cm⁻³ s⁻¹). The source
  function of the groups in local detailed balance now holds the lines'
  net emission (§2, `in_nltex_ldb`). The neutral inner wind also makes
  O I 8–13 times CMFGEN's through charge exchange. Its 2p⁴ ¹D level
  (bound-free edge 1064 Å, with routes A and B of CMFGEN's data as
  CMFGEN takes them) then blankets 912–1064 Å: κ = 1–5×10⁻¹⁵ cm⁻¹ in the
  wind, 3–40 times the electron scattering. There SuperLite's Monte
  Carlo J is zero, starved of packets: in CMFGEN's model the band holds
  ~10⁻⁴ of the J at 1300–1500 Å. A one-iteration test from the same
  state without the static metal profiles raises the 1030–1190 Å J from
  0.63–0.67 to 0.98–1.17 times CMFGEN's; level dissolution adds little at
  912–930 Å.
* Against `ies` (thermal-electron scattering, as SuperLite's) rather than
  `blk`, the H ground state of the full-metal wind agrees within 20 %
  at 1.58–2.42×10¹⁵ cm (`rayFs2`, adaptive damping, continuum
  destruction); the discrepancy is the layer just above the photosphere,
  1.1–1.4×10¹⁵ cm (b₁ 1.35 and 1.04 at 1.17 and 1.30×10¹⁵ cm against
  0.41 and 0.30), and 2.4 times at 3.6×10¹⁵ cm. Three continuations of
  `rayFs2` (six iterations each) leave the layer where it was: without the
  O I 1025.76–Lyβ overlap (`in_nltex_overlap = 0`: 1.25, 0.93; the outer
  wind falls to 35 and 471 at 2.42 and 3.6×10¹⁵ cm against `ies`' 65 and
  373), local detailed balance only for groups with τ ≥ 1000
  (`in_nltex_ldb = 1000`: 1.20, 0.90) and the Rayleigh photons of the
  Lyman lines redistributed by the thermal electrons as in `ies`
  (`in_rayleigh = 'lymes'`: 1.32, 1.06). T_e and n_e there are CMFGEN's
  (the structure is fixed). In cell 41 (1.17×10¹⁵ cm) the Lyman-continuum
  J is 0.13 times `blk`'s, which the thick continuum ties to b₁ (J = S ∝
  1/b₁); the net photoionization of the ground state is 1700 cm⁻³ s⁻¹
  against CMFGEN's 2300 (`HIPRRR`, 1.19×10¹⁵ cm), 0.5 % of the gross rate
  in both. C, O and Fe ionization in the layer agree with CMFGEN's: the
  second number of each depth record of CMFGEN's `XzVOUT` files is the
  population of the ground *level* of the next ion (C II ²P°₁/₂ is 1/3 of
  C II, Fe III ⁵D₄ 1/2.6), not of the whole ion.
  The layer has two solutions of SuperLite's equations. Started from
  CMFGEN's `ies` populations and Monte Carlo J (`cmf_m002re_hold`), with
  the physics of `rayFs2` and the populations free (`rel_cmf`, `rel_cmf2`:
  eleven iterations), it stays at CMFGEN's ionization: b₁ 0.134, 0.076,
  0.49, 0.43 at 1.07, 1.10, 1.17, 1.30×10¹⁵ cm against `ies`' 0.126,
  0.072, 0.41, 0.30 (x(H I) 6.3–6.8×10⁻⁴ at 1.17×10¹⁵ cm in every
  iteration, `rayFs2` 1.8×10⁻³). The iteration from SuperLite's own start
  ends on the neutral branch. Further out the wind drifts away from
  CMFGEN's state in both runs (b₁ 25, 159, 2040 at 2.0, 2.42, 3.6×10¹⁵ cm
  after eleven iterations, still rising; `ies` 16, 65, 373): without the
  O I–Lyβ overlap `rayFs2` gives 9.7, 35, 471 there. From CMFGEN's state
  without the overlap (`rel_cmf_ovl0`, five iterations) every radius is
  within 7–30 % of `ies`: b₁ 0.135, 0.073, 0.50, 0.39, 1.56, 14, 65, 460
  at 1.07–3.6×10¹⁵ cm (`ies` 0.126, 0.072, 0.41, 0.30, 1.73, 16, 65, 373).
  The overlap's bookkeeping is consistent: at 2.42×10¹⁵ cm (cell 50 of
  `rel_cmf2`, `cell_jtest` with the flows into O I 3d ³D°) the net
  pumping of 3d ³D° from the ground (6.2×10⁻⁴ s⁻¹ per O particle) leaves
  by 11287 Å to 3p ³P (6.3×10⁻⁴), which is the net Lyβ decay of H
  (5.2×10⁻⁷ per H particle, O/H = 8.5×10⁻⁴). There H's ground state gains
  5.2×10⁻⁷ by Lyβ and loses 3.4×10⁻⁷ by charge exchange with O⁺, 0.9×10⁻⁷
  by Lyman-continuum photoionization and 0.9×10⁻⁷ through n = 2. CMFGEN
  pumps O I as hard (3d ³D°/ground 2–3 times SuperLite's, ¹D/³P within
  20 % of SuperLite's at 1.2–3.6×10¹⁵ cm), so its Bowen loss is not
  smaller: what keeps its H ionized is on the drain side. SuperLite's J
  there is 0.25 of `blk`'s in the Lyman continuum, 0 at 912–1000 Å (no
  packets), 0.3–0.5 at 1000–1240 Å and 0.8–1.1 longward of 1440 Å.
  Which part of J matters (`cell_jtest`, cell 50 of `rel_cmf2`, H solved
  with SuperLite's rates and other elements, J replaced by `blk`'s in one
  band): x(H I) 0.051 with SuperLite's J, 0.032 with all of `blk`'s, 0.044
  with its Lyman continuum only, 0.030 with its continuum below 1250 Å
  and 0.041 with its lines below 1250 Å. The groups at 912–1250 Å,
  which ionize O I 2p⁴ ¹D and C I and so feed the charge exchange that
  drains H's ground state, carry most of it. Their J is zero there
  because no packet crosses them: beyond 1.3×10¹⁵ cm 30–47 of the 47
  groups at 912–1064 Å are empty in an iteration. `in_nltex_jfill`
  (§5) gives such groups J from their neighbours. Interpolated in log J
  between the nearest nonzero groups (`rel_jfill`, six iterations from
  `rel_cmf2`) it overestimates the band 3–11 times at 1.6–2.6×10¹⁵ cm,
  since the nonzero groups hold one or two packets each, yet cell 50 then
  gives the same x(H I) with SuperLite's J as with `blk`'s (0.031 and
  0.030) and b₁ falls to 21 and 92 at 2.0 and 2.42×10¹⁵ cm (from 25 and
  159). Beyond 3×10¹⁵ cm the whole band is empty and b₁ stays at 1830
  (`ies` 373) at 3.6×10¹⁵ cm. Smoothing every group of a sparse window
  (`rel_jbox10`, `rel_jbox30`) broke the line incident intensities,
  which subtract each line's own share from the group J: the wind went
  neutral (b₁ 3.4 at 1.07×10¹⁵ cm) or over-ionized within six iterations.
  The fill now gives only the empty groups the window mean (zeros
  included). With 10 groups (`rel_jz10`) b₁ is 22, 104, 1340 at 2.0,
  2.42, 3.6×10¹⁵ cm; with 30 (`rel_jz30`) 18, 67, 360 and 760 at
  6×10¹⁵ cm (`ies` 16, 65, 373, 1020) with the inner wind unchanged. But
  the 30-group window carries the J of 1064–1200 Å, a thousand times the
  band's, into it: J at 1000–1064 Å is 15–38 times `ies`'s beyond
  3×10¹⁵ cm. The agreement depends on the window, so the fill is a
  diagnostic, not the remedy. `in_nltex_fbias` (§5) samples the band
  instead. Biased in every cell (`rel_fb100`) it lost the luminosity
  (2×10⁴¹ of 2.4×10⁴² erg s⁻¹ emerged after one iteration): a packet of
  the dense wind meets thousands of macro-atom events, and the product of
  its weights is carried by exponentially rare packets. Biased only in
  the cells where a quarter of the band had no flight (`rel_fb100b`, 100,
  q ≤ 0.1) the energy is kept to 10⁻³, the empty groups at 1.4–2.6×10¹⁵ cm
  fall from 18–45 to 2–28 of 47, but beyond 3×10¹⁵ cm the band stays
  empty and b₁ is that of `rel_cmf2` (29, 180, 2120 at 2.0, 2.42,
  3.6×10¹⁵ cm). The outer band is empty for a physical reason: its red
  edge, 1064 Å, is the edge of O I 2p⁴ ¹D (11.65 eV), so O I absorbs
  the whole band, and its J at 1000–1064 Å falls from 1.07 of `ies`'s at
  1.6×10¹⁵ cm to 0.08 at 3.4×10¹⁵ cm. O I follows H I through the
  charge exchange (SuperLite's O I fraction 0.06, 0.14, 0.20 at 2.4, 2.9,
  3.4×10¹⁵ cm against CMFGEN's 0.024, 0.038, 0.049), so a neutral wind
  shields the band that would ionize it: the same loop as the inner
  layer's two branches. The outer wind is the neutral branch of that
  loop, reached from CMFGEN's state through the Bowen drain; what keeps
  CMFGEN on the ionized one is still open. The transport is not it: on
  CMFGEN's populations held (`cmf_m002re_hold`) SuperLite's J at
  1000–1064 Å is 1.2–1.9 times `ies`'s out to 3.4×10¹⁵ cm (empty only
  beyond 4×10¹⁵ cm, where `ies`'s is 10⁻³ of the neighbouring bands'), so
  the band goes dark in the free runs because SuperLite's own O I is
  3–4 times CMFGEN's. The rates themselves agree: started from `ies`'s
  populations (`cell_jtest` on `ies_ovl`: `ies`'s state, the J of
  `cmf_m002re_hold`, the overlap on, `SL_MAXIT=3000` so that every solve
  converges) SuperLite's statistical equilibrium keeps x(H I) within 1 %
  of `ies`'s at 1.7, 2.4, 3.1, 3.7 and 4.4×10¹⁵ cm. There the ground
  state's budget is that of CMFGEN's HIPRRR: the Bowen inflow (Lyβ,
  5.5×10⁻⁷ s⁻¹ per H) leaves by net Lyman-continuum photoionization
  (3.4×10⁻⁷, a few 10⁻³ of the gross rates) and charge exchange (1.2×10⁻⁷);
  CMFGEN's net photoionization of the ground state is also about three
  times its net charge exchange. What holds the outer wind ionized is that
  small excess of J over S in a Lyman continuum thick across each cell
  (τ ~ 10³–10⁴). With the Lyman continuum in local detailed balance
  (J = S, Λ* = 1: any `cell_jtest` mode with the external J, even with no
  external J below 1250 Å) the same cells go to x(H I) 0.043–0.13, two to
  three times `ies`'s. In the free runs those groups are in local
  detailed balance (`in_nltex_ldb = 30`, DDMC), J = S with no radiation
  from elsewhere (J_ext = 0 in `output.nltex_ldb`), while on CMFGEN's
  populations the Monte Carlo J there is 1.07–1.42 S at 1.6–4.4×10¹⁵ cm.
  Restricting local detailed balance to τ ≥ 10⁶ (`rel_ldb6`, six
  iterations from `rel_cmf2`) does not keep it: the Monte Carlo J there
  falls to 0.14–0.23 of `ies`'s and b₁ is 150 and 2570 at 2.42 and
  3.6×10¹⁵ cm (the 1.3×10¹⁵ cm layer comes closest yet, 0.344 against
  0.303). On `ies`'s populations held with every option of the free runs
  (`hold_full`: overlap, Lyman Rayleigh scattering, continuum
  destruction) SuperLite's J is within 2–10 % of `ies`'s at 1200–1500 Å
  and 504–912 Å, and `cell_jtest` keeps `ies`'s state with it (x(H I)
  0.0226, 0.0392, 0.0487, 0.0567 against 0.0225, 0.0390, 0.0483, 0.0561
  at 2.4–4.4×10¹⁵ cm): CMFGEN's outer wind is a fixed point of
  SuperLite's own physics. Switching on only the Λ* of the groups in
  local detailed balance, with the same J (mode `fuv 912`), moves the
  cells to 0.030, 0.060, 0.085, 0.110. In a Lyman continuum thick
  across the cell the ground state's own cycle drops out and what is left
  is the absorption of the other species' recombination photons there
  (O⁺ → O I ground at 910 Å, inside H's continuum, a 10⁻³ share of the
  recombinations, comparable with the net term of a few 10⁻³), which the
  ALI holds at the J of the last transport while it lets H's own part
  follow the populations. The coupling of O and H through those photons
  and the charge exchange is then solved with one side frozen. Solving it
  with both sides free (`in_nltex_ldbse`: the groups with Λ* = 1 take
  J = (κ₀J₀ + δη)/(κ₀ + δκ) from the bound-free emission and opacity of the
  current populations of all elements in every iteration of the
  statistical equilibrium) does not keep the state either: x(H I) 0.0036,
  0.038, 0.070, 0.102, 0.139 at 1.7–4.4×10¹⁵ cm (`ies` 0.0052, 0.0225,
  0.039, 0.048, 0.056), after 2000–3000 iterations. The innermost cell
  moves to the ionized side, the others to the neutral: with J frozen the
  rates reproduce `ies`'s state to 0.5–0.8 % at every radius, and with
  the thick Lyman continuum's J following the populations that mismatch
  is amplified about a hundred times, toward whichever side it has. The
  outer wind of this model sits next to the fork between the two
  branches (a loop gain close to 1): its H ionization is set by a few
  10⁻³ of the gross rates and changes by factors of 2–5 for differences of
  a few 10⁻³ in them, below the agreement of two codes' atomic data and
  discretizations. It is not a quantity on which to benchmark the
  physics; the agreement of the rates at fixed J (0.5 %) is. The layer at
  1.1–1.4×10¹⁵ cm is the same kind of fork (its two branches above); the
  wind inside 1.1×10¹⁵ cm, where the continuum and most of the lines
  form, is not near one and agrees with `ies` to 1–10 % in b₁.
  `in_nltex_ldbse` stays off: it changes the answer only through the
  conditioning and multiplies the iterations of the solve.
* The Balmer jump of the full-metal wind. Against the optical continuum
  (4000–8300 Å), SuperLite's Balmer continuum (3000–3646 Å) is 11–16 %
  too bright with Sobolev lines (against CMFGEN's Sobolev model) and
  14–20 % with the hybrid (against `ies`), and 3650–4000 Å is 4–6 % too
  faint: L(3300–3646)/L(3650–4000) is 1.35 against 1.14 and 1.29 against
  1.03. In formal integrals of 2400–4600 Å (`fi_bj_coh`, `fi_bj_dm11`)
  SuperLite's Sobolev flux falls across the edge (5.8 → 4.6×10³⁸
  erg s⁻¹ Å⁻¹) where CMFGEN's rises (5.2 → 5.5). The excess is the H
  Balmer continuum itself: the ratio of the two codes' J steps at 3646 Å
  (hybrid, 1.07–1.6×10¹⁵ cm: 1.24–1.28 blueward, 0.97–1.00 redward). In
  the layers where it forms (1.03–1.13×10¹⁵ cm, τ_es 8 → 1) SuperLite's
  photoionization rate per n = 2 atom is 6–28 % above CMFGEN's and b₂
  5–24 % below (Sobolev: 2–14 % and 1–11 %), their product within 4 %
  inside 1.11×10¹⁵ cm: in both codes the Balmer continuum ionizes n = 2,
  and the brighter continuum and the smaller b₂ (S ∝ 1/b₂) hold each
  other up. In the H/He wind it is 8–13 % above `ies` (b₂ 3–5 % below
  where it forms, 1.03–1.07×10¹⁵ cm). On CMFGEN's own populations (below)
  SuperLite's Balmer continuum in the full-metal wind agrees with CMFGEN's
  within 4 %, and its J of the Balmer continuum in the wind is 15–21 %
  above CMFGEN's with the line-list lines of the NLTE ions beyond the
  model atoms (SuperLite's default, not in CMFGEN's model), 0–7 % without
  them. Most of the excess is the energy balance of CMFGEN's comparison
  models (next item). Ruled out: the inner boundary
  (SuperLite: a reflecting surface, the injected luminosity with a Planck
  spectrum; CMFGEN: the diffusion approximation): a thermalizing boundary
  at CMFGEN's core temperature (`in_innerbc = 'thrm'`, 17 027 K;
  `cmf_m002re_dm13t`, two iterations) changes J in the two innermost
  cells only and the spectrum not at all. T_e from SuperLite's thermal
  balance with the hybrid's J (`cmf_m002re_dm13e`) is within +0.6 % of
  CMFGEN's from the base to the photosphere (one cell +1.9 %) and 0.5–5 %
  below it beyond 2.5×10¹⁵ cm: SuperLite's radiation field balances the
  gas at CMFGEN's T_e, CMFGEN's own does not (next item). The deficit at
  3650–4000 Å is in the H/He wind as well (3 % in the Monte Carlo
  spectrum): level dissolution (`in_nltex_dis`, as CMFGEN's `DO_DIS`)
  makes up 0.4 % of it; CMFGEN also gives the H and He lines Stark
  profiles in its observer's spectrum (`LIST_VGT`), which merge the high
  Balmer lines near the edge (roadmap 9b). Deep inside (τ_es = 10–30)
  SuperLite's radiation field is not thermalized at CMFGEN's T_e: J is
  0.89–0.98 of B(T_e) (bolometric), 0.70 at 912–2000 Å, 1.07–1.14 at
  5000–8300 Å, and b₁ of H 1.3 (the Lyman continuum in local detailed
  balance holds S = B/b₁), where CMFGEN has J = B(T) within 1 % in every
  band and b₁ = 1.00: its diffusion boundary sets J = B(T) at the core,
  while SuperLite carries the luminosity with its own opacities (J/H at
  τ_es = 15–30 is 4–9 % below that of CMFGEN's Sobolev model). These
  layers lie below the photosphere of every band. On CMFGEN's LTE
  populations there (below) J at 1500 Å stays 0.7–0.8 of B, with dips to
  0.3–0.6 in the strongest line regions of 1700–2300 Å: a group's Monte
  Carlo J holds the photons between the lines and those drifting into
  the thick line cores, not the photons trapped in them, which the
  macro-atom's escape probability carries.
* The energy balance of the comparison models (CMFGEN's `GENCOOL`, net
  cooling integrated over the volume). CMFGEN's radiative-equilibrium
  Sobolev models are balanced (H/He `sob` 1.4×10⁴⁰ erg/s, m002 `re`
  2×10³⁹: < 0.6 % of L). Its comoving-frame models keep their T_e, and at
  that T_e their gas absorbs more than it emits inside 1.12×10¹⁵ cm, where
  the Balmer continuum forms, by up to 60 % of the local rates: H/He `blk`
  1.8×10⁴⁰ erg/s (0.7 % of L), `ies` 1.9×10⁴¹ (8 %), m002 `blk` 1.4×10⁴¹,
  `ies` 1.6×10⁴¹ (6.6 %); the emergent luminosity of `ies` is 3.3 % (H/He)
  and 4.3 % (m002) below that of the Sobolev model. SuperLite's Monte Carlo
  re-emits what the gas absorbs (k-packets and macro-atoms at the fixed
  T_e), so it emits this energy, much of it in the Balmer continuum: on
  `ies`' own populations (H/He, `cmf_hhe_holdies`, three iterations) its
  Balmer continuum is 5 % above `ies` relative to 4000–6500 Å, and its J
  at 3000–3646 Å 4–9 % above CMFGEN's at 1.04–1.10×10¹⁵ cm, while the
  formal integral of those populations (12 scattering passes,
  `fi_holdies_q12`), which emits with their emissivities as `cmf_flux`
  does, gives CMFGEN's spectrum: 3500–3700 Å within 1 %, 3700–13 000 Å
  0–4 % above (−2.5 to +1.4 % relative to 4000–6500 Å), Hα, Hβ, Pβ 103.3,
  21.6, 56.0 Å against 100.9, 21.1, 55.2, He I 10830 −3.7 Å against −3.7
  (H/He; the full-metal wind alike, below). CMFGEN's own
  comoving-frame models are fainter there than its Sobolev model by about
  their imbalance (3000–3646 Å: `blk` 0.98, `ies` 0.92 of `sob`; the
  optical within 1 %). Against the balanced models SuperLite's Balmer
  continuum agrees within 2–4 % relative to the optical: H/He
  (`cmf_hhe_sob_ldb`, Sobolev lines) 1.03 and 1.02 at 3000–3300 and
  3300–3646 Å of `sob`'s against 0.98–1.03 at 4000–8300 Å, and the default
  configuration (`cmf_hhe_e4c`) 1.04 and 1.04; m002 like for like
  (`cmf_m002re_nb1`) 1.05 against 1.06–1.08 in the optical of `re`. Of its
  excess against `ies` (H/He 8–13 %, m002 10–14 %), the larger part is
  thus the energy CMFGEN's fixed-T_e models do not re-emit. The rest is
  the incoherent scattering of Lyα in `ies`' iterations, which ionizes H
  just above the photosphere: at 1.06–1.10×10¹⁵ cm x_HI is at most 0.006
  in `ies`, against 0.10 in `blk` and 0.20 in SuperLite (0.17 with
  `in_nltex_esrates`), so SuperLite's n₂ is 7–43 % below `ies`' in the
  layer above the thermalization depth of the Balmer continuum (within 3 %
  of `blk`'s inside 1.07×10¹⁵ cm, 7–27 % below it outside).
* SuperLite's thermal balance against CMFGEN's radiative equilibrium
  (H/He, Sobolev lines in both; `cmf_hhe_teS`: T_e from `in_tesolve`, eight
  iterations from CMFGEN's T_e of `sob`, `cmf_hhe_teF`). T_e stays within
  0.5 % of CMFGEN's in 64 of the 69 cells and within 1.8 % in all (the
  innermost cell −1.8 %, the outer wind at 3.7–4.4×10¹⁵ cm −0.7 to −1.4 %),
  the H recombination front just above the photosphere (1.065–1.24×10¹⁵ cm,
  x_HI 0.02–0.2) included. Before (`cmf_hhe_teS`) two cells of the front
  ran 13 % hot and ionized fully: the root search of the balance
  warm-started the statistical equilibrium of each trial temperature from
  another trial's, and with the Lyman continuum in local detailed balance
  (J = S of the cell) the equations of a front cell have a nearly neutral
  (x_HII 0.01) and a fully ionized solution besides the one that continues
  its state (0.8): f(T) alternated between them from one trial to the
  next (`SL_TESCAN`, §5), and the roots were spurious. Every trial now
  starts from the cell's populations. The bands of the spectrum stay
  within 1–6 % of `sob`'s, as with CMFGEN's T_e, and the formal integral
  (12 passes, `fi_teF_q12`) gives Hα, Hβ, Pβ, He I 10830 108.3, 22.9, 54.5,
  −3.0 Å against `sob`'s 96.3, 21.2, 49.4, −3.2 (+12, +8, +10 %; `ies`
  100.9, 21.1, 55.2, −3.7). The
  two hot cells had changed the bands by < 2 % (relative to the total)
  and Hα, Hβ, Pβ of the formal integral with the Monte Carlo J by +13,
  +2, +8 % (91.3, 17.3, 55.4 → 103.6, 17.7, 59.7 Å, Monte Carlo noise of J
  included). In A4's and
  A7's CSM (ionized, x_HI ≤ 10⁻³) T_e varies by 4–8 % from cell to cell, as
  the noise of the heating rates. CMFGEN with T released from `ies`
  (`sn_iin_hhe_iesT`, 50 iterations; its linearization had to be rebuilt
  with T free, `FIX_BA`) does not converge either: T rises by 0–4 % over
  the Sobolev model's, most at the recombination front and in the outer
  wind, and its emergent luminosity from 97.6 to 99.5 % of the inner one,
  but its balance overshoots (net cooling 2.9×10⁴¹ erg/s inside
  1.12×10¹⁵ cm, locally up to 97 % of the rates) and the corrections stay
  at 6–20 % per iteration. SuperLite's T_e lies within 0.5–4 % of it.
* SuperLite's transport on CMFGEN's populations (`in_nltex_hold`,
  `Tools/cmfgen2superlite_state.py`: the populations of `ies` held in
  every iteration, the hybrid's transport, three or four iterations).
  With the line-list lines of the NLTE ions beyond the model atoms
  (default, `cmf_m002re_hold`) the Lyman and He continua (228–912 Å)
  have CMFGEN's J within 4 % in every cell (with SuperLite's own
  populations 0.6–13 times: there the rates differ) and the optical
  continuum CMFGEN's J within its 8 % lower luminosity, but J is 1.2–1.7
  times CMFGEN's at 912–2000 Å and 1.05–1.25 times at 2000–4000 Å from
  the photosphere outward, with Sobolev resonances for all lines as well
  (`cmf_m002re_holds`: 1.15–1.6), so not the static profiles, and the
  emergent 2000–2600 Å is 22 % faint (relative to the total). Without them
  (`cmf_m002re_holdnb`: `in_nobbopac`, all elements being NLTE; in general
  `in_nltex_bbext = f`), as in CMFGEN's model with the same atoms, the
  emergent spectrum agrees with CMFGEN's within 5 % from 1500 Å to
  2.5 µm (1000–1500 Å 7 % bright, 3650–4000 Å 5.5 % faint; normalized to
  CMFGEN's 9 % lower total), J within 0–4 % of the optical's ratio at
  2000–4000 Å and 3–12 % above it at 912–2000 Å in the photosphere and
  wind, the Balmer continuum's J within 7 %, and Hα has EW5 83–87 Å
  against CMFGEN's 86 Å (EW5 against the continuum at 12 000–16 000 km/s;
  Hβ 9.8 against 14.5 Å: the line itself as CMFGEN's, its sidebands 10 %
  brighter with SuperLite's stronger emission of Fe II 4924 and 5018 Å).
  So the transport reproduces CMFGEN's radiation field for the same
  populations and atomic data, and the extra lines account for the
  near-UV differences of all the full-metal comparisons above. What
  remains with the same populations: the J of the Lyα region (1200–1240 Å)
  1.3–9 times CMFGEN's at 1.05–1.13×10¹⁵ cm, and the deep layers (above).
* Like for like, SuperLite's own populations without those lines
  (`cmf_m002re_nb1`: the hybrid, six iterations from `cmf_m002re_dm11`,
  populations changing by 1.3–1.5 % per iteration in the last two): the
  near UV comes to CMFGEN's (1000–3000 Å 4–8 % bright, normalized to the
  total, where 2000–2600 Å was 14 % faint), the Balmer continuum's excess
  over the optical falls from 14–20 to 10–14 %, b₂ is 2–17 % below
  CMFGEN's where the Balmer continuum forms (was 3–24 %) and b₃ within 6 %
  (was 10 %). The outer wind does not change: H I 3.1×10⁻³, 1.1×10⁻²,
  3.5×10⁻² of H at 2, 3.5, 6×10¹⁵ cm against CMFGEN's 1.1×10⁻², 4.5×10⁻²,
  7.0×10⁻² (Fe III 0.84, 0.61, 0.30 of Fe against 0.57, 0.26, 0.17).
  Equivalent widths with 12 scattering passes in the formal integral, as
  `cmf_flux` (Å, within ±2700 km/s against the continuum at
  4000–6000 km/s, 200 km/s smoothing):

  | | CMFGEN `ies` | SuperLite on CMFGEN's populations (`cmf_m002re_holdnb_fi12`) | like for like (`cmf_m002re_nb1_fi12`) | hybrid with the extra lines (`cmf_m002re_dm11_fi12`) |
  |---|---|---|---|---|
  | Hα | 87.2 | 92.6 | 101.8 | 98.5 |
  | Hβ | 17.5 | 19.6 | 22.7 | 19.8 |
  | Pβ | 50.5 | 52.0 | 53.2 | 52.7 |
  | He I 10830 | −3.2 | −3.4 | −3.2 | −3.6 |
  | He I 6678 | −1.4 | −1.7 | −1.6 | −1.8 |
  | Ca II 8542 | −0.7 | −0.7 | −0.4 | −0.5 |

  On the same populations the formal integral makes Hα 6 % and Hβ 12 %
  stronger than `cmf_flux` (which has Stark profiles and dissolved high
  levels, roadmap 9b), and SuperLite's own populations add 10 and 16 %:
  b₂ (and b₃) where the lines form, with the transport now as CMFGEN's.

Line profiles of the H/He wind (Hα; EW in Å of −2700…−300, ±300 and
+300…+2700 km/s against the continuum fitted at 4000–6000 km/s on both
sides). SuperLite: formal integrals (`in_fi_only`) of one converged state
and Monte Carlo J (three iterations from the run above, with the J
restored at the warm start), CMFGEN: `cmf_flux` of its converged model.

| | blue | core | red | total |
|---|---|---|---|---|
| CMFGEN, VTURB 1 km/s | 12.1 | 58.2 | 21.3 | 92.1 |
| CMFGEN, VTURB 10–20 km/s | 14.2 | 57.5 | 24.6 | 96.8 |
| CMFGEN with comoving-frame lines and incoherent scattering in its iterations (`ies`), VTURB 10–20 km/s | 14.1 | 62.2 | 24.8 | 101.6 |
| SuperLite, Sobolev lines, trapped-photon scattering (10–20 km/s) | 14.0 | 69.5 | 20.8 | 104.4 |
| Doppler profiles (10–20 km/s), no trapped-photon scattering | 9.0 | 62.5 | 15.2 | 87.0 |
| Doppler profiles, trapped-photon scattering, 1 km/s | 11.7 | 65.6 | 18.1 | 95.6 |
| Doppler profiles, trapped-photon scattering, 10–20 km/s | 14.0 | 67.4 | 20.8 | 102.3 |
| Doppler profiles (10–20 km/s), 12 passes of the electron-scattering J (`in_fi_nes`) with Ng's acceleration (converged to 0.07 % per pass) | 15.0 | 64.4 | 24.3 | 104.5 |

EW (Å; the windows and 200 km/s smoothing of the earlier tables) of the other
lines, CMFGEN against SuperLite's formal integral with the Monte Carlo J and
with 12 passes of the electron-scattering J:

| | Hα | Hβ | He I 5876 | He I 6678 | He II 4686 | Pβ | He I 10830 |
|---|---|---|---|---|---|---|---|
| CMFGEN (VTURB 10–20 km/s) | 96.3 | 21.2 | 0.3 | −1.5 | 0.0 | 49.4 | −3.2 |
| CMFGEN `ies` (VTURB 10–20 km/s) | 100.9 | 21.1 | −0.0 | −1.6 | −0.0 | 55.2 | −3.7 |
| SuperLite, Monte Carlo J | 94.7 | 20.7 | −1.6 | −1.6 | 1.5 | 41.7 | −7.7 |
| SuperLite, 12 passes | 104.1 | 21.9 | 0.4 | −1.6 | −0.0 | 43.2 | −2.6 |
| SuperLite, edge-aligned groups, 12 passes | 105.9 | 22.6 | 0.4 | −1.6 | 0.0 | 55.0 | −3.1 |

* Doppler profiles: beyond 2×10¹⁵ cm the wind coasts (226–245 km/s;
  dv/dr ≈ v/r / 20), so along most directions a photon stays within the
  Doppler core of a line over a large part of the wind. With Sobolev
  point resonances only the material at the photon's projected velocity
  absorbs: a 20 km/s notch at −235 km/s instead of CMFGEN's P-Cygni
  trough at −300…−190 km/s. With the thermal and turbulent widths the
  trough appears (EW of −300…−150 km/s 7.3 → 1.4 Å; CMFGEN −0.7) and the
  emission peak moves 30–40 km/s to the red, the line photons escaping
  from the red wing of the profile (0…150 km/s 24.6 → 28.0 Å; CMFGEN
  26.0); the electron-scattering wings do not change. Hβ: trough at 0.10
  of the continuum (CMFGEN 0.08), peak 3.5 (3.3).
* Microturbulence 1 → 10–20 km/s changes Hα by +2.3 / +1.8 / +2.7 Å
  (blue / core / red) in SuperLite and +2.1 / −0.7 / +3.3 Å in CMFGEN.
* The trapped-photon scattering (5.6×10³⁹ erg/s against 1.7×10⁴⁰ erg/s
  of Hα escape) happens 82 % at τ_es > 0.3, so its photons scatter again
  in the expanding wind before they escape. Added where they are first
  scattered (with one more scattering, without Doppler shifts), a third
  of them land in the line core (+4.9 Å), whereas CMFGEN's core loses
  5 Å to the trapped scattering and its red wing gains twice its blue.
  Iterating the electron-scattering J in the formal integral
  (`in_fi_nes`), which contains the trapped photons and their further
  scatterings, halves the core excess. Converged (12 passes, Ng's
  extrapolation every third pass: the last pass changes L by 0.07 %)
  the electron-scattering wings agree with CMFGEN's `ies` to 0.35 Å in
  every bin of −2700…−300 and +300…+2700 km/s (Hβ 0.25 Å), and the weak
  features of the Monte Carlo J (He I 5876, He II 4686, He I 10830: its
  noise in the scattering emissivity, ±2–5 % of the continuum) are gone.
  The first pass (the scattering emissivity of the Monte Carlo J, no
  trapped-photon terms) gives the Monte Carlo L of 3500–13 500 Å (1.318
  against 1.317×10⁴² erg/s); the converged passes give 1.373×10⁴²
  (+4 %; the continuum is 3.7–4.8 % above CMFGEN's `ies` at 5000–13 300 Å
  and 9 % at 3700 Å): with the populations fixed, the photons that
  electrons scatter out of the line cores are not taken from the lines
  (the trapped-photon terms add 1.8 % to the formal integral with the
  Monte Carlo J), as in CMFGEN's `cmf_flux` of a model with coherent
  scattering. (Normalizing J by the shell volume instead of the volume
  the rays sample had raised L by 7 %, with the J of a cell still
  changing by 0.5 % per pass.)
* The Monte Carlo J of Hα varies between iterations (the line carries
  ≈ 0.5 % of the packets): with the J of two consecutive iterations the
  red wing of the formal integral is 20.8 and 25.7 Å. Average J over
  iterations (`in_fi_navg`) or iterate it (`in_fi_nes`).
* Against CMFGEN's most complete model of this wind (comoving-frame
  lines and incoherent electron scattering in its iterations, `ies`)
  the converged formal integral agrees within 2.2 Å in each part of Hα
  (blue, core, red 15.0, 64.4, 24.3 against 14.1, 62.2, 24.8 Å) and
  within 0.35 Å in every wing bin.
  The outer-wind populations that differ from CMFGEN's Sobolev model
  (b₃ 1.3–1.9 times, b₁ 2–3 times lower beyond 4×10¹⁵ cm) are those of
  `ies` (b₁ 725 against SuperLite's 685–727 at 7.4×10¹⁵ cm, b₃ 12.8
  against 14.1): there n = 3 is fed as much by Hα absorption from the
  overpopulated n = 2 and by the cascade from n = 4 as by recombination
  (a rate breakdown of the cell at 8.2×10¹⁵ cm), and the incident J of
  the Balmer lines includes the line photons scattered by electrons,
  which `ies` has and CMFGEN's Sobolev option (continuum J only) lacks.
  Just above the photosphere (1.05–1.10×10¹⁵ cm, v ≈ 2–30 km/s) `ies`
  empties the ground state of H (b₁ 0.01 against b₂ 0.37 at 1.08×10¹⁵ cm),
  where CMFGEN's Sobolev model and SuperLite have b₁ ≈ b₂ (0.22/0.34 and
  0.36/0.34): the Lyα photons that electrons scatter out of the line bring
  the energy of the hotter layers below (a net Lyα excitation of only
  5×10⁻⁷ of the gross rate balances the collisional de-excitation of n = 2
  there; CMFGEN's `HIPRRR`/`GENCOOL`). SuperLite's iterations do not
  include that scattering (`in_nltex_esrates`), so the layer is more
  neutral (x_HI up to 0.18 against 0.003, n_e 15 % lower), with n₂ and n₃
  within 25 % of `ies`.
* Pβ and the grid: with the n = 4 (Brackett) edge inside a group b₅ was
  10–30 % below CMFGEN's where Pβ forms (1.1–2.2×10¹⁵ cm; b₅ < b₆ there,
  unlike CMFGEN's monotonic b_n) and Pβ 22 % weak (43 Å). With the edges
  as group boundaries (`Tools/make_wlgrid.py`, §2) b₂…b₈ agree with
  `ies` within 4–20 % at 1.1–6×10¹⁵ cm (b₅ within 4 %), and the formal
  integral with 12 passes gives Pβ 55.0 Å against 55.2 Å, within 0.2 Å
  in every wing bin (Hα 105.9, Hβ 22.6 Å; the J of the passes is 4 %
  above the Monte Carlo J in every cell, 6 % at 3500–4650 Å).
* Remaining: the P-Cygni trough
  is shallower than both CMFGEN models' (0.7 of the continuum against
  0.25, ≈ 2 Å). At −220 km/s the line is thick on every ray out to
  p ≈ 6×10¹⁵ cm and ends at the outer edge, so a formal solution gives
  the source function of the outermost wind there: SuperLite's, and an
  independent one (Python, with either code's populations), give I ≈ S_L.
  CMFGEN's own intensities on those rays (its `IP_DATA`, which integrate
  to its spectrum within 0.5 %) are 10⁻⁹–2×10⁻⁷, up to 100 times below the
  Hα source function of its own outermost populations (5.7×10⁻⁷ erg/s/
  cm²/Hz/sr), although its observer-frame emissivity and opacity are the
  usual ones; the difference is in CMFGEN's treatment at the outer
  boundary and was not traced further (with CMFGEN's Sobolev model the
  outer-wind b₃ is also 1.3–1.9 times lower than SuperLite's, see above;
  b₂, n_e, T_e and the Hα escape probability agree, and Lyβ pumping is
  10⁻³ of the Hα escape there).
* X-ray-irradiated wind (the He I and He II lines of CSM ionized by
  shocks). The H/He wind at CMFGEN's Sobolev T_e with CMFGEN's shock
  X-rays (`INC_XRAYS`: two components, f = 0.128, 3×10⁶ and 5×10⁶ K,
  onset 100 km/s, smooth-wind densities; 5.0×10⁴⁰ erg/s, 2 % of L) in
  both codes. CMFGEN needs its X-rays switched on slowly (`XSLOW`:
  f from 5×10⁻⁴, doubled whenever a lambda iteration changes the
  populations by < 100 %); at f = 0.5 its full linearization diverges
  (corrections of 10⁹ % at one depth), at 0.128 it converges (10⁻⁴ % from
  the ramp's state at that f). Its `cmf_flux` of this model blows up near
  Lyα and runs from 3000 Å (`F_LAM_BEG`). The X-rays ionize the wind's
  He: He II from 2×10⁻⁵ to 0.99 of He beyond 1.1×10¹⁵ cm, He III 0.6–0.8 %,
  H I down 10 times; He I 10830 goes from −3.7 Å (absorption) to 18.5 Å.
  SuperLite first kept the He ionization of the wind without X-rays
  (He II 3×10⁻⁵): the He continua are thick (local detailed balance,
  J = S), and S lacked the X-ray emissivity, so the photons the source
  emits there never reached the photoionizations. With the X-ray
  emission in the source functions (`in_xray_*` row of §2): He II 0.99,
  0.99, 0.98, 0.96, 0.86 of He at 1.3, 2, 3.5, 6, 9.5×10¹⁵ cm against
  0.99, 0.99, 0.99, 0.98, 0.92; He III 0.64–0.26 % against 0.72–0.59 %;
  b(2³S) of He I within 1–5 % beyond 1.3×10¹⁵ cm; He I 10830 18.4 Å
  against 18.5 Å. Remaining: He I 2–6 times and H I 1.4–3 times more
  abundant than in CMFGEN; at 1.1×10¹⁵ cm (just outside the photosphere)
  He is still 0.76 neutral (CMFGEN 0.01); He I 5876 1.3 Å against 3.7 Å;
  the X-rays raise the H lines more than in CMFGEN (Hα 115.6 against
  103.6 Å, Hβ 25.2 against 22.7, Pβ 76.4 against 55.2; without X-rays
  102.4, 22.5, 65.7 against 100.9, 21.1, 55.2). Eight iterations from the
  state without X-rays (formal integral with the Monte Carlo J). The
  shock cells of the `lumg` source (`in_shock_cie`) emit into thick
  groups the same way, and their emission is now in the source functions
  as well (`shc_group_frac`).

### 4.5 StaNdaRT toy01: SN Ia-like ejecta at 19.5 d

The metal NLTE in a supernova rather than a wind: the StaNdaRT toy01
model (Blondin et al. 2022; `Tools/standart2superlite.py`, 100 zones,
the energy deposition of CMFGEN's model as the `lumg` source scaled to
its L_bol = 2.40×10⁴² erg/s), NLTE Si II–IV, S II–IV, Ca II–IV, Fe I–V,
Co II–V and Ni II–V, T_e from the thermal balance (`in_tesolve`), DDMC
for the thick continua with local detailed balance, 12 iterations.

Ion fractions at 19.48 d without (`toy01_ref`) and with the fast electrons
(`toy01_nt2`: `in_nt_sf`, `in_nt_lumg` = 0.754, CMFGEN's deposition over
L_bol, 12 iterations each), against CMFGEN's. Fe, Co and Ni exist only
below about 5000 km/s in toy01; outside that the ejecta are Si, S and Ca
(55, 35 and 10 % by mass).

| v [km/s] | 2000 | 5000 | 8000 | 11 000 | 14 000 | 18 000 | 22 000 | 26 000 |
|---|---|---|---|---|---|---|---|---|
| T_e SuperLite [K] | 25 450 | 10 220 | 6884 | 4884 | 3900 | 4022 | 4202 | 3828 |
| T_e CMFGEN [K] | 18 683 | 8600 | 6223 | 5097 | 4446 | 4665 | 5644 | 6010 |
| Ca IV, no fast electrons | 0.54 | 9e-13 | 7e-24 | 1e-36 | 4e-47 | 1e-53 | 5e-51 | 3e-65 |
| Ca IV, Spencer–Fano | 0.61 | 3.9e-4 | 1.1e-4 | 7.7e-5 | 1.8e-5 | 8.9e-5 | 1.8e-5 | 2.1e-5 |
| Ca IV CMFGEN | 0.091 | 3.5e-4 | 9.5e-5 | 6.2e-5 | 4.4e-5 | 8.3e-5 | 2.0e-4 | 4.7e-4 |
| Si IV, no fast electrons | 3.7e-3 | 4.7e-5 | 2e-13 | 2e-22 | 1e-27 | 1e-35 | 1e-37 | 1e-52 |
| Si IV, Spencer–Fano | 3.7e-3 | 2.5e-4 | 1.7e-5 | 4.3e-7 | 1.2e-7 | 5.0e-7 | 1.4e-6 | 6.2e-6 |
| Si IV CMFGEN | 0.80 | 2.2e-4 | 2.5e-6 | 7.7e-8 | 9.8e-8 | 4.3e-7 | 6.2e-7 | 1.3e-6 |
| S III, no fast electrons | 1.3e-4 | 0.79 | 5.6e-4 | 6e-8 | 2e-11 | 2e-16 | 2e-17 | 1e-26 |
| S III, Spencer–Fano | 1.2e-4 | 0.80 | 1.7e-3 | 6.2e-4 | 5.0e-4 | 1.3e-3 | 8.8e-4 | 1.2e-3 |
| S III CMFGEN | 5.5e-3 | 0.019 | 6.1e-4 | 4.5e-4 | 5.0e-4 | 1.1e-3 | 1.5e-3 | 2.2e-3 |
| S IV, Spencer–Fano | 0.96 | 9.8e-5 | 5.8e-8 | 2.1e-8 | 2.0e-8 | 1.1e-7 | 1.6e-7 | 4.4e-7 |
| S IV CMFGEN | 0.97 | 2.2e-6 | 2.3e-8 | 1.4e-8 | 2.0e-8 | 9.8e-8 | 1.7e-7 | 4.7e-7 |

* The fast electrons give the higher stages of the IME layers the floors
  of CMFGEN: Ca IV, S III and S IV within a factor 1.5 from 8000 to
  18 000 km/s, Si IV within 1–7. Without them those stages are 10⁻¹³–10⁻⁶⁵.
  The dominant stages and the core (Fe, Co, Ni IV at 2000 km/s 0.95–0.97
  against CMFGEN's 0.99; III at 5000 km/s) do not change.
* Ions below 10⁻³ of the most populated target (`nxn_omion`) are left out
  of the degradation but ionized by its spectrum (`nt_solve`'s rate per
  particle; checked in `test_nontherm`). Without that, trace elements (Fe in
  the IME layers of a model with mixing) would get no floor.
* Not from the fast electrons: SuperLite's IME layers are 500–2200 K
  cooler than CMFGEN's from 14 000 km/s out (the same with and without
  them), so S I and Ca II dominate there where CMFGEN has S II and Ca III
  (S II 0.27–0.39 against 1.00 at 22 000–26 000 km/s); the core is 6800 K
  hotter (and Si IV 3.7×10⁻³ against 0.80 at 2000 km/s, while S IV
  agrees). The population change per iteration was still 20 % in the
  last iterations (`toy01_nt3` continues it).
* T_e against the other codes of StaNdaRT at 19.5 d. In the IME layers
  (14 000–26 000 km/s) SuperLite's 3800–4300 K is within their spread
  (Sedona, SuperNu, URILIGHT: 3300–5400 K); only CMFGEN warms outward, to
  6000 K. These layers are converged (the root of the thermal balance
  within 1–4 % of T_e; secant steps, `in_te_secant`, end at the same
  temperatures), and there the line radiation the gas absorbs and
  thermalizes sets T_e (collisional terms 100–1000 times the deposition).
  The core is too hot: 24 000–25 500 K at 2000 km/s against 18 700–20 800
  K. The CMFGEN-like super-level atoms (`t01`: Fe IV 1000 levels in 100
  super-levels against 50 levels) take 900–1600 K off (`toy01_t01`),
  DDMC only above τ = 30 (`in_trn_tauddmc`) another 1000–2300 K
  (`tdd30`). The source was CMFGEN's deposition scaled to its L_bol
  (×1.33), so the radiation stored at earlier times was injected as heat
  where the γ-rays deposit; with the deposition alone the core is up to
  2300 K cooler (`toy01_dep`, but L_bol too low), and with that share
  emitted as radiation (`in_lumg_rad` = 0.2457, a Planck spectrum at the
  cell's T_R, not in the heating) 750 K cooler (`toy01_rad`: the core is
  thick, so most of it is absorbed there anyway). All three together
  (`toy01_comb`, seven iterations from `toy01_t01`):

  | v [km/s] | 2000 | 3000 | 4000 | 5000 | 8000 |
  |---|---|---|---|---|---|
  | other codes [K] | 18 700–20 800 | 15 100–16 900 | 10 700–12 700 | 8600–9800 | 6100–7200 |
  | v2 atoms, τ_DDMC = 5 [K] | 25 452 | 20 934 | 16 213 | 10 221 | 6884 |
  | t01 atoms, τ_DDMC = 30, `in_lumg_rad` [K] | 20 835 | 18 471 | 14 513 | 9716 | 6296 |

  The rest of the excess was the far UV of the thick core, as in DDC10
  (§4.6): its groups had no local detailed balance, so the Monte Carlo
  J of the continua thick across a cell was the few packets emitted in
  them. With `in_nltex_ldb = 30` (eight iterations from `toy01_comb`,
  `toy01_ldb`, and four more, `toy01_ldb2`):

  | v [km/s] | 2000 | 3000 | 4000 | 5000 | 6000 | 8000 | 12 000 |
  |---|---|---|---|---|---|---|---|
  | CMFGEN [K] | 18 683 | 15 081 | 10 734 | 8600 | 7466 | 6223 | 4832 |
  | Sedona, SuperNu, URILIGHT [K] | 18 700–20 800 | 15 400–16 900 | 11 000–12 700 | 8600–9800 | 7600–8500 | 6100–7200 | 4800–5700 |
  | `toy01_comb` [K] | 20 835 | 18 470 | 14 513 | 9716 | 7804 | 6296 | 4712 |
  | `toy01_ldb2` [K] | 17 341 | 15 558 | 12 714 | 9022 | 7631 | 6305 | 4714 |

  From 3000 km/s out T_e is within the codes' spread (12 714 K at
  4000 km/s against Sedona's 12 674); at 2000 km/s it is 1300 K below
  CMFGEN's. The optical rises by 15–28 %: L at 3500–4500, 4500–5500 and
  5500–6500 Å is 3.47, 2.97, 2.01×10⁴¹ erg s⁻¹ (before 2.71, 2.41, 1.58;
  CMFGEN 4.74, 3.94, 2.51; Sedona 4.67, 4.15, 3.72). For supernova ejecta
  use super-level atoms, `in_trn_tauddmc = 30`, `in_nltex_ldb = 30`,
  `in_lumg_rad` = 1 − deposition/L_bol and `in_nltex_taumin = 1d-2` (§5).
* An earlier version of this section compared Fe II beyond 5000 km/s;
  SuperLite has no Fe there, so its fractions in those cells are only the
  initial values.

**γ-ray deposition** (`Tools/gamma_deposition.py --mc --pe`, 200 000 packets,
100 zones, against the deposition of CMFGEN in StaNdaRT). The share of the
decay power deposited outside 1000 km/s (the extent of CMFGEN's grid) and the
deposition here over CMFGEN's at a few velocities:

| model, t [d] | deposited, here / CMFGEN | 2000–14 000 km/s | 17 500 km/s | 21 700 km/s | 26 400 km/s |
|---|---|---|---|---|---|
| toy06, 14.6 | 0.961 / 1.004 | 0.95–0.98 | 1.03 | 1.26 | 1.66 |
| toy06, 19.5 | 0.908 / 0.953 | 0.95–1.05 | 0.96 | 1.10 | 1.66 |
| toy06, 31.4 | 0.714 / 0.745 | 0.94–1.00 | 1.07 | 0.86 | 1.19 |
| toy06, 41.0 | 0.562 / 0.579 | 0.95–0.97 | 0.94 | 0.98 | 1.34 |
| toy06, 197 | 0.068 / 0.071 | 0.95–0.98 | 0.96 | 0.97 | 0.99 |
| toy01, 19.5 | 0.949 / 0.974 | 0.96–1.10 | 1.03 | 1.82 | 2.49 |

The Monte Carlo traps 2–5 % less than CMFGEN overall and reproduces its
profile to 5–10 % where 99.9 % of the energy goes. Before about 20 d the
outermost, thinnest layers get up to 2.5 times CMFGEN's deposition. Gray
absorption (κ = 0.03 cm²/g) does as well from 40 d on (total within 4 %),
but before that it carries too much energy outward: 2–7 times CMFGEN's in
the outer layers of toy01 at 19.5 d. A Compton-scattered γ-ray is degraded
near where it scatters. Without photoabsorption (below ~150 keV in Fe-group
matter) the deposited share at 41 d is 0.545.

### 4.6 StaNdaRT DDC10: O-rich outer ejecta at 18.7 d

CMFGEN's delayed-detonation model DDC10 of StaNdaRT (`Tools/standart2superlite.py`,
100 zones from 874 to 30 000 km/s, the ⁵⁶Ni chain decayed to 18.72 d; CMFGEN's
deposition scaled to its L_bol = 1.35×10⁴³ erg/s, 0.8645 of it through the fast
electrons), the CMFGEN-like super-level atoms of C I–IV, O I–IV, Ne I–III, Mg I–III,
Si I–IV, S I–IV, Ca I–IV, Fe I–VI, Co and Ni II–V, 12 iterations (`ddc10_t1`).
O dominates from 18 000 km/s out (0.57–0.67 by mass), C beyond 24 000 km/s.

| v [km/s] | 3000 | 9000 | 15 000 | 18 000 | 22 000 | 26 000 |
|---|---|---|---|---|---|---|
| T_e SuperLite [K] | 20 861 | 15 275 | 9129 | 7662 | 7696 | 7906 |
| T_e CMFGEN [K] | 19 857 | 12 695 | 9328 | 7801 | 8873 | 8338 |
| T_e Sedona, SuperNu, URILIGHT [K] | 20 000–22 700 | 12 100–13 600 | 7300–7700 | 5800–6500 | 4900–5700 | 5000–5200 |
| O III SuperLite | 0.99 | 0.39 | 0.022 | 0.0066 | 0.043 | 0.048 |
| O III CMFGEN | 0.99 | 0.16 | 0.065 | 0.11 | 0.27 | 0.20 |
| C III SuperLite | 0.91 | 1.00 | 0.048 | 0.020 | 0.072 | 0.089 |
| C III CMFGEN | 0.50 | 0.99 | 0.15 | 0.57 | 0.84 | 0.58 |
| Ca IV SuperLite | 0.094 | 0.0081 | 8.6e-4 | 1.5e-4 | 1.5e-3 | 1.4e-3 |
| Ca IV CMFGEN | 0.39 | 0.039 | 0.018 | 0.021 | 0.047 | 0.061 |

n_e within 15 % of CMFGEN's everywhere. In the O-rich layers only SuperLite
and CMFGEN have the higher stages: there the codes without non-thermal
ionization at this epoch (ARTIS, Sedona, SuperNu, URILIGHT) have O III
10⁻²⁷–10⁻⁵, Ca IV 10⁻³⁵–10⁻⁶ and C III 10⁻¹⁷–10⁻⁶ beyond 15 000 km/s, so
these stages come from the fast electrons. SuperLite's are below CMFGEN's by
2–16 (O III), 3–30 (C III) and 15–140 (Ca IV), with T_e and n_e close to
CMFGEN's; the fast electrons put 9 % of their energy into ionization there
and 85 % into heat. Inside, the core is too ionized out to 9000 km/s (Si IV
0.92 against 0.17, Fe IV 0.96 against 0.30), as toy01's core is too hot.
The deficit is the Monte Carlo far UV, not the fast electrons. SuperLite's
deposition per volume is CMFGEN's to 2–7 % from 3000 to 26 000 km/s, and with
the recombination coefficients of the atoms (O III → O II 2.8×10⁻¹² cm³ s⁻¹ at
7700 K, C III → C II 7.4×10⁻¹²) the fast electrons alone give C III/C II ≈ 0.1–0.2
at 22 000 km/s, against CMFGEN's 5. In the O-rich cells beyond 15 000 km/s every
group from 150 to 912 Å had J = 0 (no packet crossed it) while the continua of
O II (353 Å), Ne II, C II (508 Å) and Ca III (243 Å) are thick across a cell
(τ ~ 10⁴): the recombinations into their ground states went unreabsorbed (case A
in place of case B) and the far-UV emission of the gas was lost. Local detailed
balance (`in_nltex_ldb = 30`, with DDMC from τ = 5; six iterations from
`ddc10_t1`, `ddc10_ldb`) gives those groups J = S:

| v [km/s] | 3000 | 9000 | 15 000 | 18 000 | 22 000 | 26 000 |
|---|---|---|---|---|---|---|
| T_e SuperLite [K] | 18 711 | 14 184 | 9118 | 7661 | 7709 | 7612 |
| O III SuperLite (CMFGEN) | 0.996 (0.99) | 0.98 (0.16) | 0.092 (0.065) | 0.098 (0.11) | 0.21 (0.27) | 0.24 (0.20) |
| C III | 0.74 (0.50) | 0.97 (0.99) | 0.066 (0.15) | 0.19 (0.57) | 0.31 (0.84) | 0.24 (0.58) |
| Ca IV | 0.035 (0.39) | 0.23 (0.039) | 0.023 (0.018) | 0.017 (0.021) | 0.043 (0.047) | 0.058 (0.061) |
| Ne III | 0.93 (0.98) | 0.85 (0.18) | 0.16 (0.093) | 0.12 (0.12) | 0.28 (0.23) | 0.35 (0.22) |
| Si IV | 0.47 (0.41) | 0.95 (0.17) | 0.015 (0.010) | 0.008 (0.008) | 0.018 (0.027) | 0.018 (0.028) |
| Fe IV | 0.99 (0.97) | 0.98 (0.30) | 0.016 (0.017) | 0.039 (0.087) | 0.075 (0.22) | 0.052 (0.089) |

From 15 000 km/s out O III, Ne III, Ca IV and Si IV are now within 1.6 of
CMFGEN's (2–140 below before), Fe IV within 3 and C III within 1.4–3. Use
`in_nltex_ldb` (30, with DDMC) for supernova ejecta as for the IIn wind. The
intermediate layer at 9000 km/s, already too ionized and too hot, is more so
with it (O III 0.98 against 0.16, Ca IV 0.23 against 0.039): J = S there is
the source function of gas 1500 K hotter than CMFGEN's, the same excess as
toy01's core (§4.5, #98).

### 4.7 CMFGEN ζ Pup-like O-star wind (T fixed)

CMFGEN's sample ζ Pup-like model (`sample_grid/zpup_like_model`: clumped
wind, f = 0.05, its X-ray source at 3 and 5 MK, f_X = 0.078), converted with
its atom set, with T fixed to CMFGEN's. `zpup_sl` (DDMC from τ = 5) ran 10
iterations; with DDMC from τ = 50 the first transport did not end (the thick
He II and Lyman continua stay IMC and packets random-walk). Ion fractions
against CMFGEN's at 1.5×10¹² cm (inner wind) and 1.2×10¹³–10¹⁴ cm (outer):

| | inner, CMFGEN | inner, SuperLite | outer, CMFGEN | outer, `zpup_sl` | outer, `zpup_aug` |
|---|---|---|---|---|---|
| He II | 0.0026 | 0.0006 | 0.50–0.99 | 0.26–0.92 | 0.34–1.00 |
| N V | 0.29 | 0.69 | 3.4–4.4×10⁻⁴ | 1.4×10⁻⁵–2.5×10⁻⁴ | 2.0–3.3×10⁻⁴ |
| N VI | 9.7×10⁻⁴ | 3.2×10⁻³ | 7–9×10⁻⁴ | 10⁻¹⁰–10⁻⁸ | 1.1–1.3×10⁻⁴ |
| O V | 0.083 | 0.38 | 2.5–3.6×10⁻⁴ | 1.9×10⁻⁵–3.9×10⁻⁴ | 2.2–4.1×10⁻⁴ |
| O VI | 3.4×10⁻⁴ | 4.7×10⁻⁴ | 0.7–3.4×10⁻⁴ | 10⁻⁹–10⁻⁸ | 3–8×10⁻⁵ |
| Ca top (Ca VIII) | 4.4×10⁻⁵ | 4.8×10⁻⁵ | 0.6–1.4×10⁻⁷ | ~10⁻¹² | 0.5–1.4×10⁻⁷ |

* Local detailed balance (`in_nltex_ldb = 30`, `zpup_ldb5`, five iterations
  from `zpup_sl`) brings the outer He to CMFGEN's (He II 0.88, 0.99 at 3×10¹³
  and 10¹⁴ cm against 0.93, 0.99; He III 0.0061 against 0.0067) and changes
  the inner wind little.
* The outer wind's highest stages are X-ray ionization with the Auger
  electron. SuperLite absorbed the X-rays in the core shells but gave the
  energy only to the fast electrons; with `in_nltex_auger` (a core-shell
  ionization moves an ion with ≥ 4 electrons up two stages, a Li-like one
  up one; four iterations from `zpup_ldb5`, `zpup_aug`) N V, O V and the top
  stages of Ca and Fe come within 1.5–8 of CMFGEN's, from 10³–10⁶ below, and
  N VI and O VI within 7 and 2–5.
* Still open: the inner wind (1.5–6×10¹² cm) is too ionized (He II 4–5, C IV
  3–5 too low; N V, O V 2–5 too high), and from 1.2×10¹³ cm out O IV stays
  0.9 where CMFGEN has 0.16–0.74 (its O III 0.25–0.84). What it is not:
  the clumping (SuperLite's filling factor, T and n_e follow CMFGEN's to
  20 %), and the He II Lyα–O III 303.80 Å overlap (`in_nltex_overlap` 150
  with `in_nltex_ovllmax` 5, 86 603 pairs; `zpup_ovl`: no change). CMFGEN's
  He II ionizes through n = 2 (photoionization 3.0 against recombination
  1.2×10¹⁰ cm⁻³ s⁻¹ at 1.5×10¹² cm), its ground state is on the spot
  (2.53 against 2.41). SuperLite's transport on CMFGEN's populations held
  (`zpup_hold`, `Tools/cmfgen2superlite_state.py`) and `cell_jtest` with
  that J (`SL_IONS`, `SL_AUGER`) give at 1.5×10¹²–3.3×10¹³ cm He II
  0.44–1.25, C IV 0.24–1.0, N V 0.12–3.2, O V 0.09–2.6 and O III 0.28–0.6 of
  CMFGEN's (with `cell_jtest` before it used the clump densities: He II
  0.04–0.15, O III 0.01–0.02, which these notes first read as a J that is
  far off). The J of CMFGEN's populations is therefore close; what remains
  is O III, 1.7–3.6 times low throughout. Its He II ground photoionization per ion is
  1.2–2.9 times CMFGEN's (PRRR over POPHE) at 1.5–12×10¹² cm. In this
  clumped wind no cell is thick in the He II continuum (no group below
  504 Å reaches local detailed balance even from τ = 5, `zpup_ldbs`, which
  only over-ionized the outer He), so that J is the Monte Carlo's: 7–34 of
  the 39 groups at 200–228 Å empty and the others up to 14 times the mean.
  Biasing the emission into 150–228 Å (`in_nltex_fbias` 30, `zpup_fb`)
  fills those groups (1–24 empty) and leaves the inner wind as it was: the
  noise is not the cause. On CMFGEN's populations SuperLite's He II
  photoionization per ion is CMFGEN's to 1 % in the photosphere (1.23–
  1.25×10¹² cm, 52 000–39 000 K; 3.3 times at the inner boundary, whose
  123 000 K Planck spectrum the photosphere thermalizes) and 1.7–2.9 times
  from 1.27×10¹² cm out, where the flow passes ~100 km/s: there the
  Sobolev point resonances let through the far UV that the Doppler
  profiles of the line forest (Fe V–VI, O IV–V, …) block in CMFGEN's
  comoving-frame transfer. With static Doppler profiles of all lines in
  the slow cells (`in_nltex_dopmc` 2, `in_nltex_dopz` 1; 44 of 83 cells;
  three solves from `zpup_aug`, `zpup_dop`) He II is 0.030, 0.19, 0.94 at
  3, 6×10¹² and 3×10¹³ cm (CMFGEN 0.025, 0.14, 0.93; before 0.007, 0.083,
  0.69), C IV 0.32 and 0.73 at 3 and 6×10¹² cm (0.50, 0.89; before 0.11,
  0.47), N V 0.41 and O V 0.22 at 1.5×10¹² cm (0.29, 0.083; before 0.69,
  0.38); O III of the outer wind doubles (0.12 against 0.69 at 3×10¹³ cm).
  For O-star winds use static profiles of all lines in the slow cells, as
  for the slow CSM of the IIn models.

## 5. Running

Build type: CMake now builds Release (-O3) unless `CMAKE_BUILD_TYPE` is given;
before, a build without it had no optimization (-O0). On toy01 (t01 atoms, 2
ranks, 2¹⁸ packets) -O3 takes 1234 s for two iterations against 1820 s at -O0
(transport 1.43 times, the NLTE update 2 times faster), with identical results
for the same random numbers. Memory: the per-line arrays of the transport are
(lines × cells) on every rank; with the super-level atoms of a supernova (3.8×10⁵
lines, 100 cells) each takes 300 MB. The tallies of `in_nltex_jbwin` and
`in_nltex_esrates` are allocated only with those options and the MPI sums are
in place: 3.4–3.6 GB per rank instead of 4.5 GB. Two such runs of two ranks each
exceeded the 16 GB container (one was killed for memory): plan about 4 GB per rank.
Active lines: the transport steps at the Sobolev resonance of every line whose
depth exceeds `in_nltex_taumin` (default 10⁻⁴) in some cell, in every cell; the
log prints how many. With the toy01 atoms that is 92 028 of 379 700 lines and
most of the transport's steps (gprof: 5×10⁸ steps for 2×2¹⁵ packets). With
`in_nltex_taumin = 1d-2` 47 379 lines remain and the transport takes 260 s per
iteration instead of 500 s (2¹⁸ packets, four iterations from `toy01_comb`,
`tau2` against `tau4`); T_e within 0.7 %, ion fractions within 1–4 % (the noise),
the emergent luminosity within 0.1 % and every band from 1500 Å to 2 μm within
2 %. Use 10⁻² for supernova ejecta; the IIn benchmarks were made with 10⁻⁴.
Noise: `in_nltex_javg` averages the Monte Carlo group J over iterations (weight of
the previous ones; off by default). In toy01 (six iterations from `tau2`, weight
0.5 against none) the median of the largest population change per cell fell from
0.19, 0.16, 0.12, 0.12 to 0.13, 0.12, 0.09, 0.11 in iterations 3–6, while the T_e
scatter between iterations stayed at 0.3–0.4 % (0.8–1.1 % in the outer half) and
outer trace stages (Ca III at 14 000–22 000 km/s) still move by ±30 %: the rest of
the noise is in the line estimators and the coupling to T_e. More packets (2¹⁹
instead of 2¹⁸) are the remedy for final models; the transport is now 2.8 times
faster than at -O0 with τ_min = 10⁻⁴.
Starved groups: `in_nltex_jfill = n` gives every group within n groups of one
that no packet crossed in a cell (J = 0 after the LDB) the mean J over the n
groups on each side, weighted by the group widths (off by default). It is meant
for the far UV of thin winds, where a zero J switches off the photoionization of
levels with edges there (§4.4: O I ¹D and C I at 912–1064 Å in the IIn wind).
Band biasing: `in_nltex_fbias = B` (band `in_nltex_fwl1`–`in_nltex_fwl2` Å,
default 912–1064) chooses the emission channels of a macro-atom state and the
k-packet outcomes that emit in the band with probability q = max(P, min(B P, 0.1))
instead of their share P, the packet weighted by P/q, the others by (1 − P)/(1 − q),
in the cells where a quarter of the band's groups had no flight in the last
iteration. Unbiased (`test_nltex`: weighted frequencies within 0.5σ); off by
default, and with it off every pick draws the random numbers it drew before.
Biasing every cell loses the energy in the weights (§4.4).

```
cmake -DUSE_MPI=ON <src> && make
python3 Tools/stella2superlite_csm.py mesa.day020_post_Lbol_max.data -o run --tau-max 50
python3 Tools/make_wlgrid.py -o run
# NLTE model atoms and charge-exchange reactions from the CMFGEN atomic data
# (atomic_data_21jun23 of https://sites.pitt.edu/~hillier/web/CMFGEN.htm);
# DIE=none leaves out the dielectronic resonances (as CMFGEN with [DIE_XzV] = F,F)
Tools/make_nlte_atoms.sh <CMFGEN atomic dir> atoms
# input.par: in_grp_ng = 0, in_flx_ndim = -1,1,1 (written grids);
#   in_nlte = t, in_nlteh = t, in_grp_ngs = 1,
#   in_nltex_dir = 'atoms', in_nltex_atoms = 'he1.atom he2.atom c1.atom ... fe4.atom'
mpirun -np 4 superlite
```
IIn models with the physics of the benchmarks (§4.4) and the A4/A7 runs
(§4.3); the switches are those of §2:
```
 in_srctype = 'lumg'                  ! L(r) of the hydro model, shocks where v drops outward
 in_shock_cie = 'shock_cie.dat'       ! shock emission as CIE spectra (Tools/chianti_cie_spectra.py)
 in_grp_ng = 0, in_flx_ndim = -1,1,1  ! groups of Tools/make_wlgrid.py, edges as boundaries
 in_src_n2s = 18                      ! 2^18 packets per iteration
 in_niter_early = 8, in_src_n2s_early = 16   ! from scratch only: 0 in warm starts
 in_thm_thermal = t                   ! thermal electrons: the scattering wings
 in_nlte = t, in_nltex_dir = 'atoms', in_nltex_atoms = 'h1.atom he1.atom ... fe6.atom'
 in_tesolve = t, in_tion = 'nebu'     ! T_e from the thermal balance
 in_trn_tauddmc = 30d0, in_nltex_ldb = 30d0
 in_nltex_vturb = 10d0, in_nltex_vturbmax = 20d0
 in_nltex_dopmc = 2d0                 ! static metal lines in slow cells
 in_nltex_dis = t                     ! level dissolution, as CMFGEN's DO_DIS
 in_fi = t
```
In a dense, slow CSM (a velocity change of about a Doppler width per cell
where electron scattering dominates, as A7's) also `in_nltex_dopz = 1`
(static H and He lines there): with Sobolev H lines the Monte Carlo and
the formal integral differ by a factor of two in Hα, with
`in_nltex_esrates` by 1.44 (§4.3). The emergent spectrum is the formal integral of the
last state with 12 scattering passes, as CMFGEN's `cmf_flux`:
`in_fi_only = t, in_fi_nes = 12, in_fi_nesdv = 60d0, in_nltex_init =
'<run>/output.nltex_state'` (it reads `output.nltex_jmc` next to it).
Warm starts (`in_nltex_init`) continue a converging model; the log prints
the source particles of every iteration.

Outputs: `output.flx_luminos` (spectrum per iteration), `output.nltex` (ion fractions
of every NLTE element and n_e per cell), `output.nlteh` (H populations,
departure coefficients, Sobolev optical depths, escape probabilities per cell),
`output.tebal_<rank>` with `in_tesolve` (heating and cooling terms per cell at the last
iteration; with the environment variable `SL_TESCAN=<cell>,<cell>,...` also
`output.tescan_<cell>`, the balance function f(T) around the old T_e in every
iteration, with n_e and the top-state fraction of the first element); with `SL_RATEDUMP=<cell>,...` `output.ratedump_<cell>`, the rates of the statistical equilibrium of those cells at every iteration (the net flows into the ground state of H from each state and by charge exchange; per Lyman line (the first six) and Hα the escape probability without and with β_es, the incident intensities of the resonance estimator, of the trapped-photon e.s. term and in the rates, the source function, the net radiative rate and, with `in_nltex_esrates`, the e.s. term of the expectation `nx_bew`; and per element and ion stage the flows up and down and the net charge exchange); with `SL_FIJCELL=1` and `in_fi_nes` also `output.fi_jcell`, the J of the formal integral per cell on its comoving grid.  Unit tests: `Tests/unit/run_tests.sh`.
