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

## Running the corrected model

A restart from the first run's populations (iteration 33) did not
converge: each restart read the stored BA matrix (`BAMAT`) of the old
state and took a stale step, and the populations, built with twice the
heating, kept swinging (at iteration 38, over 1000 of about 1350
variables per depth changed by more than 100 % at depths 1-21).  Two
rules for running it:

1. **Start fresh** (no `SCRTEMP`, `POINT1`, `POINT2`): LTE start,
   `LTE_EST = T`.  The first 15-17 iterations are LAMBDA iterations
   with T fixed; the variables changing by more than 100 % per depth
   fell from 818 to 0 by iteration 15.
2. **Delete `BAMAT` and `BAMATPNT` before every restart.**  CMFGEN reads
   them at restart and skips computing BA.

Then **set `LAM_VAL` to 1e15** (restart, it is read at start).  In the
cool inner layers (depths 65-80) a few trace levels (O V) fall almost to
zero and come back in LAMBDA iterations, where impurity species are not
fixed.  The maximum change stays at 1e5-1e11 %, above `LAM_VAL = 1e4`,
so CMFGEN never leaves LAMBDA mode and T is never freed
(`DO_T_AUTO` frees it when the maximum change falls below 50 %).  At
iteration 17 the luminosity was 2.22e9 Lsun, against 1.567e9 Lsun
deposited.  With LAMBDA iterations off, the BA iterations fix the
impurity species (`FIX_IMP = T`) and T can vary.

Each iteration that computes BA takes about 80 minutes on 2 threads
(80 depths, NT = 1596); one that reuses it, a few minutes.
