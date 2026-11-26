# D22U and D22S: twins of published CMFGEN Ibn dense-shell models

SuperLite models of `pwr2e42_he4p0_CDS_3e15_unicomp` (D22U) and
`pwr2e42_he4p0_CDS_3e15` (D22S) of Dessart, Hillier & Kuncarayakti (2022,
A&A 658, A130, arXiv:2111.13360; spectra on Zenodo, record 8069692), from
the setup of its Sect. 4.1 and its Fig. 4:

- all of the he4p0 ejecta (1.62 Msun; Ertl et al. 2020) in a dense shell
  whose density is a Gaussian in velocity, centre 2000 km/s, sigma 70 km/s;
  homologous at 173.6 d, so that the shell is at 3e15 cm; next to no mass
  below or above it (1e-4 Msun of inner ejecta, a 1e-5 Msun/yr CSM);
- 2e42 erg/s deposited as heat, Gaussian profile exp(-((v - 2000)/98.99)^2)
  per unit volume (`in_nt_lumg = 0`);
- composition: D22S the stratified profile of their Fig. 4 against velocity
  (`../D22S/he4p0_strat.prof`, read from the vector paths of the figure in the
  PDF: mass fractions of He, C, N, O, Mg, Si, Ca, Fe and the velocity of each
  mass coordinate; the elements not plotted filled from the deficit:
  Ne = 0.4 O, S = 0.5 Si, the rest as Fe, i.e. 56Co, 56Fe and Ni; Na = 1e-3 O);
  D22U its mass average (`he4p0_uni.comp`). The averages agree with the
  he4p0 totals of Dessart et al. (2021, Table 1): He 0.558 (0.568), C 0.060
  (0.062), N 0.0039 (0.0040), O 0.198 (0.191), Si 0.032 (0.029), Ca 0.0009
  (0.0008), Mg 0.024 (0.018); Fe-group 0.029 (their Fe 0.016 plus 56Co);
- the formal-solution J of 5-912 A (`in_nltex_fij = 1`);
- the decay power of their 56Co (half of the 0.0445 Msun of 56Ni,
  6.69e9 erg/g/s), 2.96e41 erg/s, fully trapped, deposited as heat in
  proportion to the Fe-group mass density (`--decay`): with it the twins
  inject 2.30e42 erg/s, as the published spectra hold (2.37e42 erg/s
  over 1000 A - 2.5 um for the stratified model, 2.50e42 for unicomp).
  CMFGEN deposits decay energy through gamma rays and fast electrons.

Not in the twins: their depth grid, model atoms and CMFGEN settings. The comparison is with their published spectra
(1000 A - 2.5 um), their He I profiles (Fig. 6) and the region of formation
(Fig. 7).

Made with (P = the directory of the composition files)

    Tools/ibn_cds_reference.py <dir> --t 173.611111 --mcds 1.62 --mej 1e-4 \
      --mdot 1e-5 --power 2e42 --decay 2.9607e41 --comp $P/he4p0_uni.comp --comp-csm          # D22U
    Tools/ibn_cds_reference.py <dir> --t 173.611111 --mcds 1.62 --mej 1e-4 \
      --mdot 1e-5 --power 2e42 --decay 2.9607e41 --comp-profile $P/he4p0_strat.prof --comp-csm  # D22S

`SN_HYDRO_DATA` is the same structure for CMFGEN.

## The deposited power as fast electrons (D22UN, D22UF, D22SF)

D22U and D22S put the deposited power into the gas as heat. Their shells
came out with He neutral (He I 1.00 in D22U; He I 0.6–0.9 in the He-rich
outer layers of D22S at 11–13 kK), while Dessart+22 find He "partially
ionized" (he4p0) and "nearly once ionized" in the He-rich shells: the He I
lines of CMFGEN's spectra (3889, 5876, 6678, 7065 A) are missing, and in
D22S the hot He-rich layers cool instead through C II] 2326, N II] 2143,
Si III] 1892, C III] 1909, C II 1335 and Mg II 2800 (Sobolev estimates from
the state, cells at 1990–2130 km/s), UV lines that CMFGEN's spectrum does
not have.

- D22UN: D22U with the decay share of the power (0.1289) as fast
  electrons (`in_nt_lumg`; exact in the uniform model, whose decay and
  shock deposition profiles are the same): He II 0.16 of the shell's He
  after five iterations.
- D22UF, D22SF: D22U and D22S with all of the deposited power as fast
  electrons (`in_nt_lumg = 1`), as the energy of the shock reaches the
  gas in CMFGEN's treatment of Dessart's injected power (inferred from the
  ionization the paper reports; the paper does not say).

Results (SuperLite / CMFGEN luminosity in 1000–3000, 3000–5000,
5000–7000, 7000–10000 A; 12 NLTE iterations each, final median change
0.24 for D22UF and 0.18 for D22SF):

| model | deposition | UV | blue | green-red | red |
|---|---|---|---|---|---|
| D22U  | heat | 6.98 | 0.98 | 0.55 | 0.47 |
| D22UN | decay share fast | 6.17 | 0.96 | 0.66 | 0.57 |
| D22UF | all fast | 7.44 | 0.92 | 0.66 | 0.58 |
| D22S  | heat | 2.40 | 0.70 | 0.79 | 0.66 |
| D22SF | all fast | 1.75 | 0.71 | 0.90 | 0.81 |

- D22SF: the He-rich outer layers become mostly once ionized (He II
  0.28–0.82, outward), 0.3–1.2 kK cooler, and He I 3889, 5876, 6678 and
  7065 appear at nearly CMFGEN's strength. The UV intercombination lines
  halve but remain (Mg II 2800, C II] 2326, N II] 2143, Si III] 1892), and
  3000–5000 A stays 30% low.
- D22UF: the shell's He II rises to 0.3–0.5 and weak He I lines appear,
  but the UV excess does not shrink (Mg II 2800 and C II] 2326 strengthen)
  and the red stays 35–40% low; in the uniform model the deposition is not
  the main cause of the difference.

## CMFGEN's shock-power width (D22SG, D22UG)

The width of the power profile is not in the text of Dessart+22, and the
twins took dv = sqrt(2) sigma = 99 km/s, which makes the power per unit
mass uniform through the shell. The dashed "normalized power deposition
profile" of their Fig. 7 (digitized from the PDF at 600 dpi,
`../D22SG/fig7_edep_digitized.txt`; a Gaussian fit gives centre 1999.7 km/s,
dv 50.2 km/s, the same below and above the peak, rms 0.035 dex) is
exp(-((v - 2000)/dv)^2) with dv = 50 km/s, i.e. DVLOC_SHOCK_POWER = 50 km/s
of CMFGEN's Gaussian shock power per unit volume (which CMFGEN passes to its
non-thermal solver with the decay energy; if the plotted curve were per
unit mass instead, the volume profile would have dv = 45 km/s). With dv = 50 the He-rich outer
layers of the stratified shell (2070-2260 km/s) receive 3-100 times less
power per unit mass than in D22SF (the decay floor remains), and the inner
shell (<1950 km/s) 3 times less; the uniform shell's wings 8 times less.
D22SG and D22UG are D22SF and D22UF with `--dvpow 50`, warm-started from
their states. In D22SF those outer layers are the ones at 10.5-13 kK with
Fe 99% Fe III, cooling through the UV intercombination lines, where
Dessart+22 find a large part of the Fe II emission.
