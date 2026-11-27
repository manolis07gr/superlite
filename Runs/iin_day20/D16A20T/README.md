# D16A20T: model A at 20 d with the hottest shocked-shell cells capped

Diagnostic twin of D16A20E (same lumg emission, fixed T). Cells 60 and 61
(50.6 and 34.5 kK in the input) are set to 25.7 kK, the temperature of
cell 59, the hottest neighbouring shell cell.

Why: in D16A20 and D16A20E the CSM is doubly ionized in He (He III) right
ahead of the shell (cells 62-72, 18 kK) and in the cold outer CSM (cells
91-120, 7.6-11 kK), although the formal-solution EUV J is on
(in_nltex_fij = 1). That removes the narrow P Cygni absorptions of the
CMFGEN spectrum (Dessart et al. 2016, Figs 5-6). If the He II continuum of
these two cells keeps the CSM doubly ionized, capping them recombines it.
