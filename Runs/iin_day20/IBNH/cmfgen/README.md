# CMFGEN run of IBNH

The CMFGEN model of the IBNH structure (`Tools/ibn_cds_reference.py`,
`Tools/cmfgen_ibn_setup.py`).  `VADAT` here is the corrected one.

## Diagnosis of the first run (33 iterations, not converging)

1. **The shock power was counted twice as heat, and drove no non-thermal
   ionization.**  With `TRT_NON_TE = T` and `INC_SHOCK_POWER = T`, the
   CMFGEN version used (cur_cmf, main program of 18-Jun-2025) adds the
   shock power to the electron energy equation in its non-thermal routine
   and again in its shock-power routine.  The non-thermal ionization and
   excitation rates are scaled by the radioactive decay energy alone
   (zero here), so none of the shock power ionizes or excites.  The
   shock-power routine itself states that non-thermal processes are off.
   In effect the first run heated the shell with 2P = 1.2e43 erg/s.
   SuperLite IBNH (`in_nt_lumg = 1`) instead passes P to fast electrons
   and splits it with the Spencer-Fano shares.  The two runs were not the
   same model.
2. **LAMBDA and linearization (BA) iterations alternated.**  `MAX_LIN = 10`
   lets a linearization change a population tenfold, which exceeds
   `LAM_VAL = 800` %.  So every BA iteration was followed by a LAMBDA one
   with T fixed, and that one undid it.  The variables swinging most were
   trace high ions of the thin CSM at depths 5-26 (2200-5800 km/s):
   O V/O VI, Fe VII and Mg IV fell by about 97 %, and N IV rose 5-7 fold.
   The shell at depths 51-54 swung by 100-250 %.  Variables changing by
   more than 100 % per depth: 32 at iteration 25, 48 at 29, 54 at 33.
   Luminosity 1.81-1.88e9 Lsun, against 1.567e9 Lsun deposited.

## Changes in VADAT

| Key | Was | Now | Why |
|---|---|---|---|
| `TRT_NON_TE` | T | F | shock power counted once, all as heat |
| `FIX_IMP` | F | T | CMFGEN fixes negligible (impurity) species until near convergence |
| `LAM_VAL` | 800 | 1e4 | no LAMBDA iteration after every linearization |

The matching SuperLite run is `Runs/iin_day20/IBNHT`: IBNH with
`in_nt_lumg = 0`, where the deposited power is heat (k-packets).  That
makes the comparison like for like.  IBNH (fast electrons) stays as the
physics run.  A CMFGEN run with non-thermal shock power would need its
non-thermal routine changed to scale its rates by the shock power too,
and the separate shock heating skipped.  That is a change to CMFGEN to
report to its author, not something done here.

## Restart

The corrected run can start from the saved populations of the first run
(SCRTEMP, POINT1, POINT2 at iteration 33): copy VADAT into the model
directory and run `run_cmfgen.sh` again.  T will move a lot at first,
because the heating halves.
