# D16A20R: model A at 20 d, outer CSM restarted as He II

D16A20T (cells 60-61 capped at 25.7 kK), warm-started from D16A20E with the
He III of cells 89-120 (the outer CSM, 7.6-11 kK) moved into the He II
ground state (Tools/nltex_state_he3to2.py).

Why: in D16A20/D16A20E the outer CSM is fully He III at 8 kK behind a He II
layer (cells 73-88) that absorbs the shell's photons below 228 A. The rate
diagnostic (Tools/cell_jtest.f90, cell 100: n_e 3.9e4 cm^-3) shows the He
ionization sustained by photons below 228 A made in the cell itself
(photoionization of excited He II by the optical/UV field, then recombination
to the ground state): a self-consistent state of very low density gas. A gas
that starts as He II absorbs those photons, so there may be a second, He II
solution (the one with CMFGEN's narrow P Cygni absorptions).
