#!/bin/bash
# Build and run the unit tests of the CSM-physics branch.
#   usage: Tests/unit/run_tests.sh [build_dir]     (needs gfortran, LAPACK/BLAS)
set -e
src=$(cd "$(dirname "$0")/../.." && pwd)
b=${1:-$(mktemp -d)}
mkdir -p "$b"
cd "$b"
FC="gfortran -O2 -fallow-argument-mismatch"
LAPACK="-llapack -lblas"
build() {   # build <test> <sources...>
  local t=$1; shift
  $FC "$@" "$src/Tests/unit/$t.f90" -o "$t" $LAPACK 2>&1 | grep -v "^$" || true
}
build test_doppler   "$src/TRANSPORT1/dopplermod.f90"
build test_compton   "$src/randommod.f" "$src/TRANSPORT1/comptonmod.f90"
build test_nlteh     "$src/nltehmod.f90"
build test_macroatom "$src/physconstmod.f" "$src/miscmod.f" "$src/groupmod.f" "$src/gridmod.f" \
                     "$src/randommod.f" "$src/nltehmod.f90" "$src/nltehtrmod.f90"
build test_bbsob     "$src/physconstmod.f" "$src/bbsobmod.f90"
build test_tebal     "$src/nltehmod.f90" "$src/tebalmod.f90"
fail=0
for t in test_doppler test_compton test_nlteh test_macroatom test_bbsob test_tebal; do
  echo "=== $t"
  ./$t > $t.log 2>&1 || { fail=1; echo "$t FAILED"; }
  grep -v IEEE $t.log || true
done
[ $fail -eq 0 ] && echo "all unit tests passed"
exit $fail
