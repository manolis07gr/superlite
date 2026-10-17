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
| `in_shock_dvmin`, `in_shock_tmax` | Shock cells (velocity drop > 50 km/s and 2 % per cell, ≥ `in_shock_dvmin` per shock) emit their injected energy as a radiatively cooling flow: T uniform between T_gas and min(T_postshock, `in_shock_tmax`) (isobaric cooling releases the same energy per unit T), bremsstrahlung shape. Its EUV/X-rays photoionize the CSM and the dense shell. | `SOURCE/shock_temperatures.f90` |
| `in_shock_cie` | At each T of the cooling flow the photons follow the spectrum of gas in collisional ionization equilibrium of the cell's composition (lines, free-free, free-bound, two-photon; 10^4.5–10^8 K) instead of bremsstrahlung: at 10⁵–10⁷ K line cooling exceeds bremsstrahlung 10–100 times and comes out in the EUV/soft X-rays and the UV/optical coronal lines. Per-element tables from CHIANTI 11 (Dere et al. 1997; Del Zanna et al. 2021) made with `Tools/chianti_cie_spectra.py` (ChiantiPy; the table is not distributed). | `shockciemod.f90`, `SOURCE/interior_source.f90` |
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
| `Tests/benchmarks/wind` | pure-scattering Sobolev line in constant-velocity, β-law and homologous winds against the formal solution with direction-dependent escape (`sobolev.py`) | χ²/dof = 1.17, 1.34 (v = const, τ₀ = 3, 30), 0.98, 0.99 (β law), 1.20 (homologous); isotropic re-emission (stage 1): 5.15 (8σ in the trough) |
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
| (with the atoms) | Lines of the line list of NLTE ions whose upper level lies above the model atom, and lines of the top ion, stay in the Sobolev expansion opacity, the lower level populated from the NLTE population of the highest model level at or below it with Boltzmann excitation at T_e: the truncated Fe-group atoms (35–50 % of the ionization energy) keep their UV blanketing. | `GAS/physical_opacity_nlte.f` |
| (with the atoms) | Core-shell X-ray absorption of the NLTE ions (Verner & Yakovlev 1995; the model atoms describe the valence shell), absorbed into fast photoelectrons. | `bfxsmod.f`, `GAS/physical_opacity_nlte.f` |
| (with `in_bbsobolev`) | Rayleigh scattering by H I in the ground state (Kurucz 1970 fit, capped at Lyα), in the scattering opacity; CMFGEN's `INC_RAY`. | `GAS/physical_opacity_nlte.f` |
| `in_tesolve` | Thermal balance with all NLTE elements: photoionization heating, recombination and collisional cooling from the model atoms, charge exchange, the LTE opacities of the remaining elements, Compton, deposition; every element re-solved at the first trial T, then those with at least 1 % of the heating and cooling terms (a held population gives the wrong T dependence of radiatively pumped levels). | `tebalmod.f90`, `GAS/te_update.f90` |
| `in_niter_early`, `in_src_n2s_early` | Fewer source packets in the first iterations, while the gas state converges. | `superlite.f90` |
| `in_xray_table`, `in_xray_ff1/2`, `in_xray_t1/2`, `in_xray_v1/2`, `in_xray_smooth` | X-ray emission of shocks embedded in a wind as CMFGEN's `INC_XRAYS`, for winds and CMFGEN benchmarks. Up to two components, with filling factor f_k, shock temperature T_k and onset velocity V_k, emit η_ν = Σ_k e^{−V_k/v} f_k² Λ_ν(T_k) n_e,X N_atom. Λ_ν comes from CMFGEN's plasma table (`RS_XRAY_FLUXES` format; APEC spectra in `atomic/misc/gal_xray_hr.dat`, linear in log T). n_e,X counts H and He as fully ionized and 6 electrons per heavier atom; the densities are those of the smooth wind (`X_SM_WIND`). The power is mechanical, added to the injected energy and not to the thermal pool. | `xraysrcmod.f90`, `SOURCE/sourceenergy.f90`, `SOURCE/interior_source.f90` |
| `fcl` column of `input.str` | Microclumping as in CMFGEN (Hillier & Miller 1999): the gas sits in clumps filling a fraction f of each cell's volume, with empty space between them. The gas state (EOS, NLTE populations, n_e, thermal balance, opacities and emissivities per unit volume) is that of the clump density ρ/f. The transport takes the volume means, f times these. Sobolev depths are computed with the volume-mean opacity, f K/Q, both in the statistical equilibrium and in the expansion opacity of the LTE lines, so saturated lines are not diluted. Processes ∝ ρ² (recombination, free-free, collisional cooling) are therefore enhanced by 1/f at fixed mean density, and electron scattering is not. Clumped CMFGEN models convert with `Tools/cmfgen_model.py` (RVTJ holds the clump density); clumpy CSM. | `inputstrmod.f`, `GAS/`, `nltextrmod.f90`, `superlite.f90` (`clump_scale`) |
| `in_nltex_ali` | Accelerated lambda iteration of the NLTE bound-free transitions (Rybicki & Hummer 1991, local operator), as CMFGEN converges them. In an optically thick continuum the Monte Carlo J of a cell is mostly the cell's own emission, so iterating the statistical equilibrium against it moves an ionization front by about one optical depth per iteration. For each bound-free pair, the part of J that the cell's own recombination into the level makes, Λ* η/κ, is taken out of the photoionization. The recombinations that the level re-absorbs (share n_l σ/κ) are taken out of its recombination. Λ* = 1 − (1 − e^{−τ})/τ, from the volume-mean absorption depth across the cell. The corrections are equal at the populations with which J was made, so the solution does not change. The on-the-spot cycle drops out, and the populations follow the radiation from elsewhere at once. | `nltexmod.f90` (`nx_el_rates`), `nltextrmod.f90` (`nxt_ali_cell`) |
| `in_nltex_ldb` | Local detailed balance in the groups whose volume-mean absorption depth across a cell is at least this (e.g. 30; 0: off). The Monte Carlo J of such a group is the cell's own emission, sampled by the few packets emitted in it. Where the group carries a small share of the cell's emission (the Lyman continuum of a cool wind) the estimate is mostly zero, and the on-the-spot cycle (recombination to the ground state, re-absorbed in place) becomes net recombination: case A in place of case B. There the rates use J = Λ* S + J_ext. S = (κ_th B + η)/κ is the source function of the populations with which the packets ran. Λ* = 1 − (1 − 2E₃(τ))/(2τ) is the share of its own emission that a slab of depth τ keeps. J_ext is the radiation from elsewhere: the packets on IMC steps outside the cell of their last thermal or bound-free emission (DDMC packets in a lump are local). The same Λ* takes the cell's own cycle out of the rates, as in the ALI, so the iteration has the fixed point of the exact J and needs no Monte Carlo samples of the local emission. With DDMC (`in_puretran = f`, `in_trn_tauddmc` ≤ `in_nltex_ldb`) these groups are DDMC groups, whose radiation is the diffusion limit: there Λ* = 1 and J = S (to second order in 1/τ), as no external radiation is tallied in them. DDMC packets in a cell's lump of thick groups are shared among them as the source spectrum (S Δν), not given to the group the packet carries. `output.nltex_ldb` lists per cell the parts Λ*S and J_ext and the Monte Carlo estimate at 228–504 and 504–912 Å. | `nltextrmod.f90` (`nxt_ldb`, `nxt_lsldb`, `nxt_ali_cell`), `particle_advance.f90` |
| `in_nltex_vturb`, `in_nltex_vturbmax` (km/s; the scattering with `in_thm_thermal`), `in_nltex_esrates`, `in_nltex_esdamp` | Electron scattering of line photons trapped in their Sobolev resonance zones. A photon of a thick line is absorbed and re-emitted about 1/β times before it escapes, while its comoving frequency drifts across the line core, |x| < x_c ≈ √ln(τ/√π) Doppler widths. That is a path of several Sobolev lengths v_D/\|Q\|, on which an electron can scatter it out of the line with a Doppler shift of hundreds of km/s: to the red it escapes, to the blue the flow brings it back to the line. Point resonances (the Sobolev approximation) have no such path. In a IIn wind σ_es v_D/\|Q\| ≈ 0.1–0.4 where Hα forms, so the trapped scattering is as likely as a scattering on the whole free path, and it builds the red shoulder of the Balmer lines (CMFGEN's comoving-frame transfer of the observer spectrum includes it). Per emission its probability is β_es = σ_es v_D ⟨y(τ_μ)⟩_μ/K × exp(−K Γλ/(4π²⟨\|Q\|⟩v_e)), with y = τG(τ) for complete redistribution over a Doppler profile (`nx_ytrap`; tabulated, and checked by a Monte Carlo of a resonance in `Tests/unit/test_estrap.f90`). The exponential is the escape of the scattered photon, shifted by v_e = √(2kT/m_e), through the Lorentz wings of the line (radiative damping Γ; off with `in_nltex_esdamp = f`, Doppler profiles as CMFGEN's iterations): Lyα and the low Lyman lines of a dense wind keep their photons. v_D² = 2kT/m + v_turb², with v_turb linear in v from `in_nltex_vturb` (v = 0) to `in_nltex_vturbmax` (outer edge), as CMFGEN's VTURB_MIN/MAX. By default the formal integral adds the emissivity n_u A β_es hν/4π of every line, redistributed with the angle-averaged thermal-electron kernel, to the resonances, with the populations of the iterations, which do not include it: as CMFGEN's observer spectrum adds the incoherent scattering to iterations with coherent scattering. CMFGEN's iterations with incoherent scattering change n₃ by < 10 % and Hα by 4 % in the H/He benchmark, so the populations need not see it. With `in_nltex_esrates` (experimental) the statistical equilibrium, the macro-atom and the transport have it as well: escape probability β + β_es with the incident intensity (β J_b + β_es J_es)/(β + β_es), J_es the group J averaged over the electron kernel about the line (photons scattered into the line core), and a line photon scattered by a thermal electron from the line centre (comptonmod kernel, new direction) with probability β_es/(β + β_es). In the H/He benchmark this does not converge and over-ionizes the wind: 2–8 times more than CMFGEN with incoherent scattering in its iterations after four iterations, 10–150 times after six with CMFGEN's Doppler profiles (`in_nltex_esdamp = f`); Hα 198 Å against 102 Å. Sharing a line's photons between the line and its wings in the transport alone (rates with β) keeps the populations but takes Pβ from 56 to 34 Å. The red shoulder depends on the microturbulence (CMFGEN H/He wind, EW of +300–2700 km/s on the air scale: 21, 25 and 38 Å for 1, 10–20 and 50 km/s; §4.4). The microturbulence also sets the Doppler widths of the lines in the formal integral (`in_fi_doppler`), and with `in_fi_nes` the trapped photons come from the iterated J of the formal integral instead of β_es. | `nltexmod.f90` (`nx_ytrap`, `nx_es_table`, `nx_bes`, `nx_wq`, `nx_jes`), `nltextrmod.f90`, `TRANSPORT1/nxemitmod.f90` (`nx_linees`), `formalmod.f90` |
| `in_nltex_init` | Warm start: the `output.nltex_state` of a run on the same structure and model atoms (populations, and T_e with `in_tesolve`, matched by exterior mass as `in_nltex_prev`), and its Monte Carlo radiation field, `output.nltex_jmc` next to it (the J of the groups, averaged over this run's groups, and the incident J of the lines; written with the state). The first iteration keeps the populations instead of solving with J = B(T_R) and makes its macro-atom tables with the saved J, so it is an iteration of the saved state: a rerun with other output settings (formal integral, its line profiles) needs one iteration. Without the J file (runs before it existed) the first iteration's tables have J = B(T_R) and no line incident intensity: the line re-emission of its transport differs and its Monte Carlo J is not that of the state (the Hα electron-scattering wings of a formal integral with it were a third to a half of the converged ones); do not average it into the formal integral. | `nltextrmod.f90` (`nxt_read_init`, `nxt_read_jmc`, `nxt_write_jmc`, `nxt_gas_solve`), `superlite.f90` |
| `in_fi`, `in_fi_wl1`, `in_fi_wl2`, `in_fi_dv`, `in_fi_npc`, `in_fi_navg`, `in_fi_doppler`, `in_fi_nes`, `in_fi_only` | Formal integral of the emergent spectrum after the last iteration, as CMFGEN's `cmf_flux`: noise-free line profiles at the spacing `in_fi_dv` [km/s] between `in_fi_wl1` and `in_fi_wl2` [cm] (`output.fi_spectrum`: λ [Å], L_λ, L_ν). The observer-frame invariant I/ν³ is integrated along rays of impact parameter p through the spherical grid, with the comoving source functions of the last iteration (first order in v/c). Continuum: absorption and thermal emission of the LTE opacities, the recombination emissivity of the NLTE levels, electron scattering of the Monte Carlo J redistributed by the thermal electrons with the angle-averaged kernel of Thomson scattering (isotropic, or dipole with `in_trn_thm_aniso`), the scattering part of the LTE expansion opacity. The Monte Carlo noise of J is that of the scattering emissivity, so of the continuum where electron scattering dominates (≈ 3 % per 100 km/s group with 2¹⁸ packets): `in_fi_navg` averages J over the last iterations of a converged model. NLTE lines: the source function of the level populations (complete redistribution) with Doppler profiles, the thermal width of the ion and the microturbulence of `in_nltex_vturb` (`in_fi_doppler`, default): opacity K/(√π v_D) e^(−x²) per unit length, which gives the Sobolev depth K/\|Q(μ)\| where the flow crosses the profile quickly and, in a slow flow, the depth along the path through the profile (the narrow P-Cygni absorption of a coasting wind, §4.4); a line whose profile a ray crosses within a tenth of the cell width is a point resonance there (fast flows: supernova ejecta, at the cost of the Sobolev integral); with `in_fi_doppler = f` Sobolev resonances wherever the comoving frequency crosses a line. With `in_nltex_vturb` also the electron scattering of the photons trapped in the lines (above), added where they are scattered first. `in_fi_nes` > 0 (with the Doppler profiles) iterates the electron scattering instead, as CMFGEN's `cmf_flux`: the comoving J of each cell on the grid of the spectrum (extended by the kernel and the flow) starts from the Monte Carlo J; each pass integrates the rays of all frequencies with the scattering emissivity of the previous J and tallies the new J along them (parallel rays sample every direction of a spherical shell: J = V⁻¹ Σ_p 2πp dp ∫ I dz, core rays also behind the inner boundary, with V the volume the rays sample, Σ_p 2πp dp times their chords in the cell, so that a uniform I gives J = I: the shell volume differs from it by 0.2–7 %, by 70 % in a thin innermost cell, an excess the scattering amplifies pass by pass). Inside a line J approaches its source function, so the trapped photons and their further scatterings through the flow (with its Doppler shifts) are in J, and the scattering emissivity has no Monte Carlo noise; the trapped-photon terms above are then left out. Each pass costs a formal integral; Ng's extrapolation every third pass; 12 passes converge the H/He IIn wind (τ_es ≈ 5) to 0.07 % of L and 0.05 % of the J of any cell per pass. `output.fi_jes`: per cell the frequency-integrated J of the last pass over the Monte Carlo J. `in_fi_only` computes only the formal integral of the saved state and Monte Carlo J of `in_nltex_init` (no transport; one iteration): a rerun with other spectrum settings, as CMFGEN's `cmf_flux` after a model. Rays through the inner boundary (`in_fi_npc`, uniform in p²) start with the Monte Carlo J of the innermost cell; two per cell outside. L_ν = 8π² ∫ I p dp; the frequencies are shared among the MPI ranks. Line wavelengths are in air above 2000 Å (row 9 of §1), so a CMFGEN observer spectrum (c/ν: vacuum) must be converted to air before profiles are compared (83 km/s at Hα). | `formalmod.f90` |
| `in_nltex_prev`, `in_nltex_dtprev` (with `lumg`) | Time-dependent statistical equilibrium of the NLTE elements over a sequence of models (Dessart & Hillier 2008, CMFGEN's time-dependent terms). Each model writes the populations of all states per cell (`output.nltex_state`). A later model reads them, matched by exterior mass, which is Lagrangian when both models end at the same mass shell, and solves dx/dt = M x per element as a backward-Euler step over dt: (M − 1/dt) x = −x_prev/dt, with particle conservation exact. Recombination therefore lags behind a falling radiation field, and ionization behind a rising one. The change of ionization and excitation energy, Σ n_el ε_s (x_prev − x)/dt, is part of the cell's injected energy, since the hydro L(r) contains what the gas releases. The macro-atom emits it, activated in the states that empty, so the stored ionization comes out as recombination lines and continua (the Hα of SNe II on the plateau). It is added where it exceeds the injected energy, and taken from that energy where the gas ionizes; the thermal pool heats only on the rest. | `nltexmod.f90`, `nltextrmod.f90`, `SOURCE/sourceenergy.f90`, `SOURCE/interior_source.f90`, `GAS/te_update.f90` |

Validation (`Tests/unit/run_tests.sh` with `SL_ATOMS` = the converted atoms):

| Test | Checks | Result |
|---|---|---|
| `test_nltex` | built-in H against the stage-1 H atom; Case B; the He I macro-atom and the H macro-atom against the statistical-equilibrium emissivities; every converted atom (14 elements, 44 ions, 2005 states; with N V, O V–VI and super-levels: 450–1100 states) in a Planck field at high density at 8000, 25 000 and 10⁵ K; charge exchange at LTE | rates to 2×10⁻³, LTE populations to 7×10⁻⁵; Hα/Hβ = 2.74 (Storey & Hummer 2.86 l-resolved), He I 5876/4471 = 2.53; macro-atom He I 2.34 ± 0.03 against 2.37, Hα/Hβ 3.15 ± 0.02 against 3.14; all elements Saha–Boltzmann to < 2 %; no net charge-exchange flow into any state to 10⁻¹⁴ |
| bound-free level choice (`test_nltex`) | all atoms at LTE (25 000 K): the chosen levels at 3000, 900, 500 and 200 Å against n_l σ_l, with the default lists (2×10⁵ choices) and with lists of 4 levels (5×10⁴; a walk of the other levels in 21 % of the choices) | χ² = 39.8, 4.1, 45.1, 167.8 for 30, 4, 62, 182 degrees of freedom; with lists of 4: 25.3, 2.4, 33.3, 89.9 for 18, 2, 28, 105 |
| formal integral (`formalmod.f90`) | m002 (§4.4, IMC everywhere) with `in_fi` over 3500–13 500 Å at 20 km/s, against the Monte Carlo spectrum of the same iteration | L of 3500–5500, 5500–7500, 7500–13 500 Å within +3.1, +1.7 and +0.9 % of the Monte Carlo spectrum, which it follows without its noise |
| X-ray source (`xraysrcmod.f90`) | emitted X-ray luminosity of CMFGEN's ζ Pup model (f = 0.078 at 3 and 5 MK, onset 500 km/s, smooth-wind densities, microclumped wind) against CMFGEN's `OBSFLUX` | > 0.1 keV: 3.80×10³³ against 3.77×10³³ erg/s; > 1 keV: 3.76×10³² against 3.75×10³² erg/s (1 % and 0.3 %) |
| ALI (`test_nltex`) | built-in H at 10⁴ K, n_H = 10¹⁰: at the populations with which J was made the solution must not move; a lambda iteration of the Lyman continuum (J = J_ext + Λ* S, Λ* = 0.999) with and without ALI | fixed point kept to 3×10⁻⁷; converged in 2 iterations instead of 68, x(H II) the same to 7×10⁻⁶ |
| Local detailed balance (`test_nltex`) | slab share Λ* against E₃ by quadrature (τ = 0.1–30); J of a one-cell grid (thick groups, DDMC lump share); the toy cell of the ALI test with a thick Lyman continuum (cell 10¹³ cm): J = J_ext + Λ*S, with the Monte Carlo estimate missing the cell's own emission in the groups above τ = 30 | Λ* to 1×10⁻⁷; J exact; x(H II) 0.816 with the estimate alone (case A), 0.97447 for the exact J (125 lambda iterations) and with local detailed balance (4 iterations), the same to 10⁻¹⁴ |
| DDMC through the macro-atom (H/He wind, §4.4) | CMFGEN's H/He IIn-like wind (Sobolev lines, converged, radiative-equilibrium T_e), SuperLite at CMFGEN's T_e with DDMC (`in_trn_tauddmc` = 5) and local detailed balance, 6 iterations from the same state; Monte Carlo spectrum of the last two iterations against CMFGEN's observer-frame spectrum | EW of Hα 27 → 101 Å (CMFGEN 96), Hβ 0.7 → 21 Å (21), Pβ 8 → 29 Å (49; the Monte Carlo IR is noisy, 13–40 Å per iteration); the H I fraction stays within a factor 3 of CMFGEN's. With the continuum fallback of the lines' incident intensity: Hα 82, Hβ 16, Pβ 48 Å, the formal integral 81, 20, 43 Å. IMC everywhere: Hα 84, Hβ 16 Å, but H 3–200 times more ionized than in CMFGEN |
| time-dependent SE (`test_nltex`) | built-in H in a Case B nebula: the steady state with itself as the previous model, dt → 0, and one backward-Euler step of recombination in the dark (10⁴ K, n_H = 10⁴ cm⁻³, dt = 10⁸ s) against the root of x = x₀/(1 + α n_H x dt) | steady state kept to 1.4×10⁻⁷; frozen to 4×10⁻¹²; x(H II) = 0.7627 against 0.7626 |
| `test_tebal` | Planck field with LTE opacities, H, He I and He II must return T₀; pure-H nebula; fast-electron shares | T₀ to 0.2 %; 21 500 K; Shull & van Steenberg values |
| `test_superlevel` (`SL_ATOMS_SL`, `SL_ATOMS_FULL`) | Fe II–IV (400, 300, 300 levels) with CMFGEN's super-levels (136, 55, 51) against the same atoms with every level a state: photospheric, diluted nebular (W = 0.1), dense and cool cells | Fe ion fractions and Fe II level populations to ≤ 7 % (nebular), ≤ 2.5 % otherwise; the same atoms with super-levels and the CMFGEN-matched C, N, O, Ca, Fe set pass the LTE and charge-exchange checks of `test_nltex` |
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
2b. ~~Super-levels (CMFGEN's F_TO_S)~~ (stage 3); ~~low-temperature dielectronic recombination~~ (CMFGEN's `die_data` as cross-section resonances); high-temperature (Burgess) dielectronic recombination at T ≳ 5×10⁴ K for the ions without such data, and explicit autoionizing levels (CMFGEN's `auto_data`); non-thermal ionization and excitation of the metals (Spencer–Fano).
3. ~~Time-dependent ionization~~ (`in_nltex_prev`: statistical equilibrium and release of the ionization energy over model sequences; the time-dependent thermal energy of the gas comes with the hydro L(r) of the `lumg` source).
4. Non-thermal deposition of ⁵⁶Co γ-rays (Kozma & Fransson 1992) with the same Shull & van Steenberg channels.
5. Shock structure: adiabatic vs radiative from the cooling time; ~~cooling-function spectra (line cooling at 10⁵–10⁷ K) instead of bremsstrahlung shapes~~ (CIE spectra, `in_shock_cie`; non-equilibrium ionization behind the shock still missing); T_e ≠ T_i behind fast shocks.
6. ~~NLTE macro-atoms for Fe II, Ca II, O I, Mg II~~ (stage 3; the ε treatment of stage 2 remains for the LTE elements and the lines beyond the model atoms).
7. DDMC: the NLTE lines (Sobolev resonances exist only in IMC, so a DDMC group carries no NLTE line opacity) and thermal Compton broadening in optically thick cells. ~~Macro-atom re-emission and non-thermal deposition~~: with the NLTE macro-atoms every absorption in a DDMC group (the cell's lump of thick groups, by the rate κ S of each group) goes to the macro-atom or the thermal pool as on IMC steps, and X-ray photons make fast photoelectrons. A re-emission in a thick group of the cell keeps the packet in DDMC; any other leaves as an IMC packet from a random point of the cell, on the line it was emitted in. Before, DDMC re-emitted the absorbed energy with the continuum emissivity κ_LTE B + η of the cell, without the lines: the energy absorbed in the thick continua (the H Lyman and He continua of a wind are DDMC groups in every cell) never cascaded through the recombination lines, and the H/He benchmark lost three quarters of its Hα (§4.4). ~~NLTE source functions~~: with the NLTE macro-atoms DDMC weights its groups with the source function S = (κ_LTE B + η)/κ of the absorption opacity instead of the Planck function at T_e (η: the recombination emissivity of the NLTE levels from the Milne relation; `ddmc_spectra`). In LTE S = B. In a thick NLTE continuum the Planck weights forced J → B(T_e): the Lyman continuum of a wind, where S = B/b₁ with b₁ ≪ 1. For line work in dense CSM keep DDMC for the thick continua only: `in_trn_tauddmc` = `in_nltex_ldb` (e.g. 30), so that the optical and UV lines stay IMC (A4/A7 have only 4 cells with Δτ_es > 5, τ_es ≈ 50–65 in total) and the Lyman and He continua are DDMC groups in local detailed balance. IMC everywhere (`in_puretran = t`) with local detailed balance does not converge in the H/He benchmark and over-ionizes its wind (§4.4).

**Stage 3 — dynamics and geometry**
8. Time-dependent radiation transport (SuperNu's IMC with census) with light-travel time, for interaction phases where t_diff ≈ t_dyn.
9. 2D/3D CSM (disks, clumps, asymmetric shells) with the SuperNu 2D/3D transport.
9a. Continuum absorption of trapped line photons (Hummer & Rybicki 1985 line-photon destruction: the same trapped path as `in_nltex_vturb`, with the absorption opacity at the line, to the bound-free macro-atom or the thermal pool, and the photoionization it drives in the rates). In the H/He benchmark it matters only below 1.07×10¹⁵ cm, in the near-LTE layer above the photosphere (trapped Lyα photoionizing n = 2: β_abs/β > 1), and is below 0.2 % of β in the wind. Partial redistribution and pressure broadening of the Lorentz wings (only radiative damping now blocks the escape of scattered photons, as for Lyα).
9b. Formal integral as CMFGEN's `cmf_flux`: ~~Doppler line profiles~~ (`in_fi_doppler`), ~~electron scattering with the J of the formal solution~~ (`in_fi_nes`, Ng's extrapolation every third pass). Still: an accelerated lambda operator for the scattering passes (the inner wind of a IIn, τ_es ≈ 5, converges slowly), the same scattering J in the statistical equilibrium (CMFGEN's incoherent scattering in its iterations), and the depth of the narrow P-Cygni trough of the H/He benchmark (0.7 of the continuum against CMFGEN's 0.25: CMFGEN's intensities on those rays fall below its own outer-wind source function, §4.4).

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

| CMFGEN model | lines | T_e | state |
|---|---|---|---|
| m002: H, He, C, N, O, Ca, Fe | comoving frame | held at its starting structure (≈ 5000 K beyond 2×10¹⁵ cm) | inner wind converged; beyond ≈ 2×10¹⁵ cm the H I ground state changes by factors of several per iteration |
| H/He | Sobolev | the same | outer wind oscillates by ≈ 10 %, ionized |
| H/He | Sobolev | radiative equilibrium (released after 25 iterations; ≈ 7400 K in the outer wind) | converged (< 10⁻⁴) |
| H/He | comoving frame | that T_e, fixed | converged |
| H/He | comoving frame, incoherent electron scattering in the iterations | that T_e, fixed | converged |

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
* m002 with all elements: being rerun with these changes.

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
* Remaining: Pβ 14 % weaker: b₅ is 10–30 % below CMFGEN's where it
  forms (1.1–2.2×10¹⁵ cm; b₅ < b₆ there, unlike CMFGEN's monotonic
  b_n), so its source function is 15–35 % lower; the incident J of the
  infrared lines of n = 5 is not starved (0.1–0.25 B); the cause is the
  grid (the n = 4 edge inside a group, see `Tools/make_wlgrid.py` in §2:
  with the edge as a boundary b₅ agrees within 5 %). The P-Cygni trough
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

### 4.5 StaNdaRT toy01: SN Ia-like ejecta at 19.5 d

The metal NLTE in a supernova rather than a wind: the StaNdaRT toy01
model (Blondin et al. 2022; `Tools/standart2superlite.py`, 100 zones,
the energy deposition of CMFGEN's model as the `lumg` source scaled to
its L_bol = 2.40×10⁴² erg/s), NLTE Si II–IV, S II–IV, Ca II–IV, Fe I–V,
Co II–V and Ni II–V, T_e from the thermal balance (`in_tesolve`), DDMC
for the thick continua with local detailed balance, 12 iterations.

| v [km/s] | 2000 | 4000 | 6000 | 8000 | 10 000 | 15 000 | 20 000 | 25 000 |
|---|---|---|---|---|---|---|---|---|
| T_e SuperLite [K] | 18 962 | 11 145 | 7477 | 6111 | 5169 | 4194 | 4303 | 4895 |
| T_e CMFGEN [K] | 18 683 | 10 734 | 7466 | 6223 | 5423 | 4335 | 5156 | 6027 |
| Fe II SuperLite | 1.5e-10 | 4.4e-4 | 1.1e-4 | 0.010 | 0.33 | 1.00 | 0.094 | 7.6e-5 |
| Fe II CMFGEN | 6.7e-11 | 1.6e-4 | 0.16 | 0.83 | 0.97 | 0.99 | 0.97 | 0.97 |

n_e agrees within 6 %. Band luminosities of the last iteration against
CMFGEN's (Sedona and SuperNu, LTE, in brackets): 2500–3500 Å 0.73 (0.88,
0.44), 3500–4500 Å 0.79 (0.99, 0.67), 4500–5500 Å 0.83 (1.05, 1.09),
5500–6500 Å 0.84 (1.48, 1.38), 6500–8000 Å 0.96 (1.18, 1.30), 8000–10 000 Å
1.03 (1.89, 1.75).
* Fe, Co and Ni stay doubly ionized out to 8000–10 000 km/s, where CMFGEN
  has them singly ionized from 6000 km/s on, and beyond 20 000 km/s, where
  SuperLite's T_e is 15–20 % below CMFGEN's, they are doubly ionized again.
  Not yet traced (the same direction as the Fe II deficit of the inner IIn
  wind, ×2–3, §4.4, where H I and O I are low by the same factor: there
  CMFGEN's ground state of H is 1.5–4 times SuperLite's near the
  photosphere while n₂ and n₃ agree within 10 %).
* CMFGEN's higher stages (Ca IV, Co IV, Ni IV) have floors of 10⁻⁶–10⁻³
  in the outer ejecta from the non-thermal ionization by the γ-ray
  deposition, which SuperLite does not have for the metals (§3, items 2b
  and 4).

## 5. Running

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
Outputs: `output.flx_luminos` (spectrum per iteration), `output.nltex` (ion fractions
of every NLTE element and n_e per cell), `output.nlteh` (H populations,
departure coefficients, Sobolev optical depths, escape probabilities per cell),
`output.tebal_<rank>` with `in_tesolve` (heating and cooling terms per cell at the last
iteration).  Unit tests: `Tests/unit/run_tests.sh`.
