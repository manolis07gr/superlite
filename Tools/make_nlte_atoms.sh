#!/bin/bash
# Make the SuperLite NLTE model atoms from the CMFGEN atomic data
# (https://sites.pitt.edu/~hillier/web/CMFGEN.htm, atomic_data_21jun23.tar.gz).
#   usage: Tools/make_nlte_atoms.sh <CMFGEN atomic dir> <output dir> [ions...]
# ions: names as in the table below (default: all), e.g.  he1 he2 c2 c3
# Each line of the table: name Z stage directory nlev (levels kept, lowest
# first) and, for older data sets, the osc/phot/col file names.  nlev may be
# N:file, with file the CMFGEN super-level assignment (F_TO_S) of the data
# set: the atom then has super-levels.  The environment variables
#   SL_TABLE  a file with a table of the same format replacing the one below
#   GFMIN     keep only the lines with g_l f >= GFMIN
#   GFCUT     CMFGEN's line cut (GF_CUT; with GFLEV = GF_LEV_CUT, default 9,
#             MINTRANS = MIN_TRANS, default 10, for Z >= ATCUT, default 10)
# select other atoms.
set -e
here=$(cd "$(dirname "$0")" && pwd)
A=${1:?CMFGEN atomic data directory (the one containing HYD, HE, CARB, ...)}
out=${2:?output directory}
shift 2
mkdir -p "$out"
table="
he1 2 1 HE/I/15jul15     69 hei_osc phot_hei heicol.dat
he2 2 2 HE/II/5dec96     30 he2_osc.dat he2phot.dat he2col.dat
c1  6 1 CARB/I/19apr23   26
c2  6 2 CARB/II/19apr23  40
c3  6 3 CARB/III/19apr23 50
c4  6 4 CARB/IV/19apr23  30
n1  7 1 NIT/I/19apr23    40
n2  7 2 NIT/II/19apr23   45
n3  7 3 NIT/III/19apr23  40
n4  7 4 NIT/IV/19apr23   30
o1  8 1 OXY/I/19apr23    51
o2  8 2 OXY/II/19apr23   50
o3  8 3 OXY/III/19apr23  50
o4  8 4 OXY/IV/19apr23   30
ne1 10 1 NEON/I/19apr23  40
ne2 10 2 NEON/II/19apr23 40
ne3 10 3 NEON/III/19apr23 30
na1 11 1 NA/I/19apr23    40
na2 11 2 NA/II/19apr23   20
mg1 12 1 MG/I/19apr23    40
mg2 12 2 MG/II/19apr23   40
mg3 12 3 MG/III/19apr23  20
si1 14 1 SIL/I/19apr23   40
si2 14 2 SIL/II/19apr23  40
si3 14 3 SIL/III/19apr23 40
si4 14 4 SIL/IV/19apr23  30
s1  16 1 SUL/I/19apr23   40
s2  16 2 SUL/II/19apr23  50
s3  16 3 SUL/III/19apr23 40
s4  16 4 SUL/IV/19apr23  30
ca1 20 1 CA/I/19apr23    40
ca2 20 2 CA/II/19apr23   40
ca3 20 3 CA/III/19apr23  30
fe1 26 1 FE/I/19apr23    60
fe2 26 2 FE/II/19apr23   100
fe3 26 3 FE/III/19apr23  80
fe4 26 4 FE/IV/19apr23   50
fe5 26 5 FE/V/19apr23    120
fe6 26 6 FE/VI/19apr23   80
co2 27 2 COB/II/19apr23  80
co3 27 3 COB/III/19apr23 80
co4 27 4 COB/IV/19apr23  60
co5 27 5 COB/V/19apr23   140
co6 27 6 COB/VI/19apr23  120
ni2 28 2 NICK/II/19apr23 60
ni3 28 3 NICK/III/19apr23 60
ni4 28 4 NICK/IV/19apr23 60
ni5 28 5 NICK/V/19apr23  120
ni6 28 6 NICK/VI/19apr23 140
"
[ -n "$SL_TABLE" ] && table=$(cat "$SL_TABLE")
want=" $* "
roman="I II III IV V VI VII VIII IX X"
echo "$table" | while read -r name z stage dir nlev osc phot col; do
  [ -z "$name" ] && continue
  if [ $# -gt 0 ] && [[ "$want" != *" $name "* ]]; then continue; fi
  d="$A/$dir"
  osc=${osc:-osc_data}; phot=${phot:-phot_data_A}; col=${col:-col_data}
  sl=()
  if [[ "$nlev" == *:* ]]; then
    sl=(--ftos "$d/${nlev#*:}")
    nlev=${nlev%%:*}
  fi
  [ -n "$GFMIN" ] && sl+=(--gfmin "$GFMIN")
  if [ -n "$GFCUT" ] && [ "$z" -ge "${ATCUT:-10}" ]; then
    sl+=(--gfcut "$GFCUT" --gflev "${GFLEV:-9}" --mintrans "${MINTRANS:-10}")
  fi
# further photoionization routes (to excited terms of the next ion) and the
# levels of the next ion (same data date if present, else the latest)
  extra=()
  for r in B C D E; do
    [ -f "$d/phot_data_$r" ] && extra+=(--phot-b "$d/phot_data_$r")
  done
  if [ ${#extra[@]} -gt 0 ]; then
    el=${dir%%/*}; date=${dir##*/}
    nr=$(echo $roman | cut -d' ' -f$((stage+1)))
    nosc="$A/$el/$nr/$date/osc_data"
    [ -f "$nosc" ] || nosc=$(ls -d "$A/$el/$nr"/*/osc_data 2>/dev/null | tail -1)
    if [ -f "$nosc" ]; then extra+=(--next-osc "$nosc"); else extra=(); fi
  fi
  python3 -I "$here/cmfgen2superlite_atom.py" --z "$z" --ion "$stage" \
    --osc "$d/$osc" --phot "$d/$phot" --col "$d/$col" --hyd "$A/HYD/I/5dec96" \
    "${extra[@]}" "${sl[@]}" --nlev "$nlev" -o "$out/$name.atom" 2>&1 | grep -v "^note: [0-9]* photo\|^note: no "
done
# charge-exchange reactions among the NLTE elements (SuperLite reads
# chgexch.dat from the model-atom directory, see in_nltex_chgexch)
cx="$A/misc/chg_exch_data.30oct14"
[ -f "$cx" ] && python3 -I "$here/cmfgen2superlite_chgexch.py" "$cx" "$out/chgexch.dat"
