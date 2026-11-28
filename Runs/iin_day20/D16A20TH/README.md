# D16A20TH: model A at 20 d with T held fixed as in CMFGEN

D16A20 with `in_lumg_thermal = t`: each cell emits its own thermal
emissivity at the fixed T (the total rate of the k-packet outcomes: LTE
opacities, NLTE recombination, collisional cooling, lines in super-levels)
times its volume, not normalized to in_L_bol, and a packet absorbed into the
thermal pool is lost. This is what "T held fixed" means in CMFGEN (Dessart
et al. 2016, Sect. 2), so the two codes solve the same problem.

Why: with the lumg source, 42% of the luminosity (the L(r) jump of the
hydro model) is emitted by two thin cells outside the CDS (cells 61-62,
34.5 and 18.6 kK, 30-70 times their own thermal emissivity). Their EUV
overpopulates He II n=2 in the CSM by 1e3-1e7 against a CMFGEN run of the
same structure, and the CSM becomes He III (CMFGEN: He II).
Approximation: thermal energy reached inside a macro-atom chain (collisional
de-excitation) is still re-emitted.
