#!/bin/bash
# Uniqueness of the statistical equilibrium of chosen cells of a run, as
# the run is and without charge exchange: Tools/cell_jtest with SL_UNIQ=t
# (the solution with the run's own J from its saved state, from LTE and
# from itself), the ion window and zones kept within a solution as the
# iin_day20 kit.  The run dir needs output.nltex_state and output.nltex_jmc.
#   usage: Tools/cell_uniq_test.sh <build dir with cell_jtest> <run dir> <cells, e.g. 12,40> <out dir>
#   writes <out dir>/uniq_base.txt and <out dir>/uniq_nocx.txt, and prints the summary lines
set -e
b=$(cd "$1" && pwd); run=$(cd "$2" && pwd); cells=$3
mkdir -p "$4"; out=$(cd "$4" && pwd)
mkdir -p "$out/j"
for c in $(echo "$cells" | tr ',' ' '); do printf '1e14 1e-5\n1e16 1e-5\n' > "$out/j/cmfj_$c.dat"; done
# the run without charge exchange: links to its files, its input.par with in_nltex_chgexch = 'none'
rm -rf "$out/nocx"; mkdir "$out/nocx"
for f in "$run"/*; do ln -s "$f" "$out/nocx/"; done
rm "$out/nocx/input.par"
awk '/^ *\/ *$/ && !d {print " in_nltex_chgexch = '"'"'none'"'"'"; d=1} {print}' "$run/input.par" > "$out/nocx/input.par"
cd "$out"
export SL_UNIQ=t SL_IWIN=1d-12 SL_OVLFIX=t OPENBLAS_NUM_THREADS=1
"$b/cell_jtest" "$run" "$cells" j > uniq_base.txt 2>&1 || echo "(base: cell_jtest exit $?)"
"$b/cell_jtest" "$out/nocx" "$cells" j > uniq_nocx.txt 2>&1 || echo "(nocx: cell_jtest exit $?)"
for v in base nocx; do
  echo "== $v"
  grep -a "uniq cell\|ierr" "uniq_$v.txt" | grep -v "^ *#"
  grep -a "ions" "uniq_$v.txt" | head -3
done
