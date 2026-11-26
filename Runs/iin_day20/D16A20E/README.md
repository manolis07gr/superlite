# D16A20E: D16A20 with the shock-spike emission moved into the CDS

A diagnostic variant of D16A20 (Dessart et al. 2016 model A at 20 d). In
D16A20 the lum column puts 1.27e43 erg/s (half of L_bol) into cells 59-61:
the thin post-shock spike at 26-50 kK (their Fig. 3, T held fixed), where the
'lumg' source emits it with the local emissivity of 50 kK gas. The result is
EUV-rich. D16A20's CSM (cells 62-120) ends up He III throughout, even at
8 kK. Its spectrum therefore shows a strong broad He II 4686 that CMFGEN's
20 d spectrum (their Fig. 5) lacks. The narrow CSM P Cygni absorptions of
Halpha and He I 6678 (their Fig. 6) are missing for the same reason.

CMFGEN, with T fixed, has no such normalization: each cell emits its thermal
emissivity, and for the optically thin, fully ionized 50 kK layer that is
about 1e42 erg/s, not 8e42. Here the emission of cells 59-61 is added to the
CDS cell 58 (20 kK); every other cell, T and the total are unchanged.
