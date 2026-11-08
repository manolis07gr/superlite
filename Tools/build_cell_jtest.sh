#!/bin/bash
# Build Tools/cell_jtest against the objects of a SuperLite CMake build,
# with the compiler and libraries of that build's own link line.
#   usage: Tools/build_cell_jtest.sh <CMake build dir>   ->  <build dir>/cell_jtest
set -e
src=$(cd "$(dirname "$0")/.." && pwd)
b=$(cd "$1" && pwd)
link=$(cat "$b/CMakeFiles/superlite.dir/link.txt")
fc=${link%% *}
libs=""
set -- $link
while [ $# -gt 0 ]; do
  case "$1" in
    -framework) libs="$libs $1 $2"; shift;;
    -l*|-L*|-Wl,*|*.dylib|*.so|*.so.*|*.a|*.tbd) libs="$libs $1";;
  esac
  shift
done
objs=$(find "$b/CMakeFiles/superlite.dir" -name '*.o' ! -name 'superlite.f90.o')
$fc -O2 -fallow-argument-mismatch -ffree-line-length-none -I"$b" "$src/Tools/cell_jtest.f90" $objs $libs -o "$b/cell_jtest"
echo "built $b/cell_jtest"
