#!/bin/bash
# A4 and A7 at day 20 with the IIn CSM physics of CSM_PHYSICS.md (#83):
# NLTE model atoms, CIE shock spectra, a warm start from the saved state in
# <model>/init, six iterations, then the emergent spectrum as a formal
# integral with 12 electron-scattering passes (as CMFGEN's cmf_flux).
#
# usage: Runs/iin_day20/run.sh <CMFGEN atomic dir> <work dir> [ranks] [models]
#   CMFGEN atomic dir: atomic/ of atomic_data_21jun23 (CMFGEN web page);
#     "-" when ATOMS gives the model atoms already made
#   ranks: MPI ranks per model (default 4; ~3 GB of memory each)
#   models: default "A4 A7" (run one after the other)
# environment (all optional):
#   SL_BUILD  a CMake build of this branch (with MPI); default: build one in
#             <work dir>/build
#   ATOMS     a directory with the model atoms; default: make them in
#             <work dir>/atoms (Tools/make_nlte_atoms.sh with table_sn3.txt
#             and CMFGEN's line cut GF_CUT = 1e-3; identical to the atoms of
#             the saved states)
#   CIE       the CIE table of the shocks; default: make <work dir>/shock_cie.dat
#             with Tools/chianti_cie_spectra.py (ChiantiPy, XUVTOP set to the
#             CHIANTI database; about an hour)
#   JMC       a directory with A4/output.nltex_jmc and A7/output.nltex_jmc, the
#             Monte Carlo estimators of the saved states: without them the
#             first iteration builds its macro-atom tables with J = B(T_R)
#   NITER     iterations (default 6; the run stops before when the NLTE
#             populations have converged, in_nltex_tol)
#   INIT      a work dir of an earlier run: each model starts from its
#             <INIT>/<model>/output.nltex_state and output.nltex_jmc (to
#             continue that run) instead of <model>/init
#   EXTRA     namelist settings added to input.par of the iterations, e.g.
#             EXTRA="in_nltex_overlap = 0d0" (tests; use another work dir)
#   NOFI      if set, no formal integral (tests)
#   MPIRUN    the MPI launcher (default mpirun; e.g. "srun" on a cluster, or
#             "mpirun --mca btl self,sm" on one machine without the TCP warnings)
# Each model runs in <work dir>/<model> (run.log, output.*), its spectrum in
# <work dir>/<model>_fi (output.fi_spectrum).  A model whose run.log ends
# with "exit 0" is skipped, so the script can be rerun after a stop; an
# interrupted model starts again from init/ (set its in_nltex_init to its
# own output.nltex_state by hand to continue from the last iteration).
set -e
[ $# -ge 2 ] || { sed -n '2,30p' "$0"; exit 1; }
src=$(cd "$(dirname "$0")/../.." && pwd)
kit=$src/Runs/iin_day20
if [ "$1" = "-" ]; then
  atomic=-
  [ -n "$ATOMS" ] || { echo "atomic dir '-' needs ATOMS"; exit 1; }
else
  atomic=$(cd "$1" && pwd)
fi
mkdir -p "$2"; work=$(cd "$2" && pwd)
np=${3:-4}
models=${4:-"A4 A7"}
niter=${NITER:-6}
mpirun=${MPIRUN:-mpirun}

# 1. SuperLite
if [ -z "$SL_BUILD" ]; then
  SL_BUILD=$work/build
  mkdir -p "$SL_BUILD"
  (cd "$SL_BUILD" && cmake -DUSE_MPI=ON "$src" > cmake.log && make -j "$np" superlite > make.log)
fi
exe=$SL_BUILD/superlite
[ -x "$exe" ] || { echo "no $exe"; exit 1; }

# 2. model atoms
if [ -z "$ATOMS" ]; then
  ATOMS=$work/atoms
  if [ ! -f "$ATOMS/fe6.atom" ]; then
    echo "$(date +%H:%M) model atoms -> $ATOMS"
    SL_TABLE=$kit/table_sn3.txt GFCUT=1e-3 bash "$src/Tools/make_nlte_atoms.sh" "$atomic" "$ATOMS" > "$work/atoms.log" 2>&1
  fi
fi

# 3. CIE spectra of the radiative shocks
if [ -z "$CIE" ]; then
  CIE=$work/shock_cie.dat
  if [ ! -f "$CIE" ]; then
    echo "$(date +%H:%M) CIE table -> $CIE"
    python3 -I "$src/Tools/chianti_cie_spectra.py" -o "$CIE" --logt 4.5,8.0,0.1 \
      --cache "$work/cie_cache" > "$work/cie.log" 2>&1
  fi
fi

setup() {   # setup <run dir> <model> <init state>
  local d=$1 m=$2 init=$3 f
  mkdir -p "$d"
  cp "$kit/$m/input.str" "$kit/$m/input.wlgrid" "$kit/$m/input.fluxwl" "$d/"
  for f in Atoms Tinit.dat data.bf_verner data.ff_sutherland data.ion data.zeta; do
    ln -sfn "$src/Data/$f" "$d/$f"
  done
  sed -e "s|@ATOMS@|$ATOMS|" -e "s|@CIE@|$CIE|" -e "s|@INIT@|$init|" \
      -e "s|^\( *in_niter *= *\).*|\1$niter|" -e '$d' "$kit/$m/input.par" > "$d/input.par"
  [ -n "$EXTRA" ] && printf ' %s\n' "$EXTRA" >> "$d/input.par"
  [ -n "$NOFI" ] && printf ' in_fi = f\n' >> "$d/input.par"
  printf '/\n' >> "$d/input.par"
}
done_ok() {   # done_ok <run dir>: the run ended with exit status 0
  [ -f "$1/run.log" ] && [ "$(tail -n 1 "$1/run.log" | tr -d ' ')" = "exit0" ]
}
run() {   # run <run dir>
  (cd "$1" && OPENBLAS_NUM_THREADS=1 $mpirun -np "$np" "$exe" > run.log 2>&1; printf '\nexit %d\n' $? >> run.log)
  done_ok "$1" || { echo "$(date +%H:%M) $1 failed, see run.log"; exit 1; }
}

for m in $models; do
  d=$work/$m
  if done_ok "$d"; then
    echo "$(date +%H:%M) $m done already"
  else
    mkdir -p "$d/init"
    if [ -n "$INIT" ] && [ -f "$INIT/$m/output.nltex_state" ]; then
      cp "$INIT/$m/output.nltex_state" "$d/init/"
      [ -f "$INIT/$m/output.nltex_jmc" ] && cp "$INIT/$m/output.nltex_jmc" "$d/init/"
    else
      cp "$kit/$m/init/output.nltex_state" "$d/init/"
      [ -n "$JMC" ] && [ -f "$JMC/$m/output.nltex_jmc" ] && cp "$JMC/$m/output.nltex_jmc" "$d/init/"
    fi
    setup "$d" $m "$d/init/output.nltex_state"
    echo "$(date +%H:%M) $m: $niter iterations on $np ranks"
    run "$d"
  fi
# the emergent spectrum: formal integral of the last state
  [ -n "$NOFI" ] && continue
  f=${d}_fi
  if done_ok "$f"; then continue; fi
  setup "$f" $m "$d/output.nltex_state"
  sed -e "s|^\( *in_niter *= *\).*|\11|" -e '$d' "$f/input.par" > "$f/input.tmp"
  printf ' in_fi_only = t\n in_fi_nes = 12\n in_fi_nesdv = 60d0\n/\n' >> "$f/input.tmp"
  mv "$f/input.tmp" "$f/input.par"
  echo "$(date +%H:%M) $m: formal integral"
  run "$f"
done
echo "$(date +%H:%M) all done: <work>/<model>_fi/output.fi_spectrum"
