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
