# CMFGEN inputs of model A at 20 d (the same structure as D16A20)

SN_HYDRO_DATA from `Tools/dessart/cmfgen_modelA.py ../input.str <dir> 20.0`
(the cell centres of input.str: T held fixed, the non-monotonic V, rho, the
mass fractions); atomic-data links as Tools/cmfgen_ibn_setup.py with the
table Runs/iin_day20/table_sn5.txt (plus HI) and the Fe II F_TO_S of the Ibn
CMFGEN run.

The configuration of Dessart et al. (2016, Sect. 2): T fixed ([FIX_T]),
Sobolev rates ([GLOBAL_LINE] SOB), the continuum without velocity terms
([NOV_CONT], which allows the arbitrary velocity field; it needs the
non-relativistic moments, so [USE_J_REL], [INCL_REL], [INCL_ADV_TRANS] F
and [STOP_IF_BP] F to pass the SN-model checks), diffusion inner boundary,
no thick outer boundary, no inserted boundary points, and their model atom
(H I, He I-II, C I-III, N I-III, O I-II, Na I, Si II, Ca II, Fe II-V; no
Sc II, Ti II here). Convergence from the LTE start: full linearization at
every iteration ([LAM_VAL] 1e30, [N_FIX_BA] 0), corrections limited to a
factor 2 ([MAX_LIN] 2); [NUM_BNDS] 1 (disk: the BA matrix of 1627 levels
and 3 bands is 7 GB).
