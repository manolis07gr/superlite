# D22A: twin of a published CMFGEN Ibn dense-shell model

SuperLite model of `pwr2e42_he2p6_CDS_3e15` of Dessart, Hillier &
Kuncarayakti (2022, A&A 658, A130; spectra on Zenodo, record 8069692), from
the setup of its Sect. 4.1:

- all of the ejecta (0.79 Msun, he2p6) in a dense shell whose density is a
  Gaussian in velocity, centre 2000 km/s, sigma 70 km/s; homologous at
  173.6 d, so that the shell is at 3e15 cm; next to no mass below or above
  it (here 1e-4 Msun of inner ejecta and a 1e-5 Msun/yr CSM);
- 2e42 erg/s deposited as heat with a Gaussian profile
  exp(-((v - 2000)/98.99)^2) per unit volume (`in_nt_lumg = 0`);
- the mixed composition of the He-star explosion he2p6 (Ertl et al. 2020,
  mixed as in Dessart et al. 2021): `he2p6_mix.comp`, the mass-weighted mean
  over the sub-shells of Table A.3 of Dessart et al. (2021, arXiv:2109.12350).
  With 56Ni decayed to half 56Co and half 56Fe; Co folded into Fe (no model
  atoms); Ar, K, Ti, Cr and Ni left out (0.6 % of the mass). Element totals
  as their Table 1 (He 0.90, O 0.031, Si 0.0046, Ca 3.4e-4 by mass);
- the formal-solution J of 5-912 A (`in_nltex_fij = 1`): the thin outer
  layers are starved of EUV packets otherwise.

Not in the twin: their radioactive decay power (a few per cent of the shock
power), their depth grid, model atoms and CMFGEN settings. The comparison is
with their published spectrum (1000 A - 2.5 um, `in_fi_wl1/2`).

Made with

    Tools/ibn_cds_reference.py <dir> --t 173.611111 --mcds 0.79 --mej 1e-4 \
      --mdot 1e-5 --power 2e42 --comp he2p6_mix.comp --comp-csm

`SN_HYDRO_DATA` is the same structure for CMFGEN.
