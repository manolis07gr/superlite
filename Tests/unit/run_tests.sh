#!/bin/bash
# Build and run the unit tests of the CSM-physics branch.
#   usage: [SL_BUILD=<cmake build dir>] [SL_ATOMS=<model-atom dir>] \
#          [SL_ATOMS_SL=<super-level atoms> SL_ATOMS_FULL=<the same without>] \
#          Tests/unit/run_tests.sh [work_dir]
#   (needs gfortran, LAPACK/BLAS; test_nebular links against the objects of
#   a SuperLite CMake build given by SL_BUILD, with mpif90 for MPI builds;
#   the He tests of test_nltex and test_tebal need he1.atom and he2.atom in
#   SL_ATOMS, made with Tools/make_nlte_atoms.sh; test_superlevel compares
#   atoms with super-levels against the same atoms without, see its header)
set -e
src=$(cd "$(dirname "$0")/../.." && pwd)
b=${1:-$(mktemp -d)}
mkdir -p "$b"
cd "$b"
FC="gfortran -O2 -ffree-line-length-none -fallow-argument-mismatch"
LAPACK="-llapack -lblas"
build() {   # build <test> <sources...>
  local t=$1; shift
  $FC "$@" "$src/Tests/unit/$t.f90" -o "$t" $LAPACK 2>&1 | grep -v "^$" || true
}
build test_doppler   "$src/TRANSPORT1/dopplermod.f90"
build test_compton   "$src/randommod.f" "$src/TRANSPORT1/comptonmod.f90"
NLTEX="$src/physconstmod.f $src/miscmod.f $src/groupmod.f $src/gridmod.f $src/randommod.f \
       $src/nltehmod.f90 $src/nlteatommod.f90 $src/nltexmod.f90 $src/nltextrmod.f90"
build test_nlteh     "$src/nltehmod.f90"
build test_nltex     $NLTEX
build test_bbsob     "$src/physconstmod.f" "$src/bbsobmod.f90"
build test_tebal     $NLTEX "$src/tebalmod.f90"
build test_superlevel $NLTEX
build test_shockcie  "$src/shockciemod.f90"
tests="test_doppler test_compton test_nlteh test_nltex test_bbsob test_tebal test_superlevel test_shockcie"
if [ -n "$SL_BUILD" ]; then
  objs=$(find "$SL_BUILD/CMakeFiles/superlite.dir" -name "*.o" ! -name "superlite.f90.o")
  LD=gfortran; command -v mpif90 >/dev/null && grep -q "USE_MPI:BOOL=ON" "$SL_BUILD/CMakeCache.txt" && LD=mpif90
  $LD -O2 -fallow-argument-mismatch -I"$SL_BUILD" "$src/Tests/unit/test_nebular.f90" $objs $LAPACK -o test_nebular
  ln -sf "$src/Data/data.ion" "$src/Data/data.zeta" .
  tests="$tests test_nebular"
else
  echo "(SL_BUILD not set: test_nebular skipped)"
fi
fail=0
for t in $tests; do
  echo "=== $t"
  args=""
  case $t in
    test_nltex|test_tebal) [ -n "$SL_ATOMS" ] && args="$SL_ATOMS";;
    test_superlevel) [ -n "$SL_ATOMS_SL" ] && args="$SL_ATOMS_SL $SL_ATOMS_FULL";;
  esac
  ./$t $args > $t.log 2>&1 || { fail=1; echo "$t FAILED"; }
  grep -v IEEE $t.log || true
done
[ $fail -eq 0 ] && echo "all unit tests passed"
exit $fail
