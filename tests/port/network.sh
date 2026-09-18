#!/bin/zsh
# Comparison harness for legacy/network.cc (auxin/PIN network models).
#
# Run from a directory holding s.rk5, euler.rk5 and seeded inits built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
# Cell variables 0-7 and wall variables 1-10 carry non-zero values there.
#
# Forms that store the membrane PIN into a wall variable write it from derivs,
# so they differ from legacy at the t=0 print only - see NOTES.md.
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-twoSquareT.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }

P12="0.1
0.02
0.05
0.3
1.5
0.2
0.4
0.1
0.6
0.15
0.25
0.08"

printf "%-46s " "AuxinModelSimple1"
mk "AuxinModelSimple1 12 1 4
$P12
0 1 2 3"; run
printf "%-46s " "  ... 267 cells"; run s.rk5 meristemT.init

printf "%-46s " "AuxinModel1 (geometry weighted)"
mk "AuxinModel1 12 1 4
$P12
0 1 2 3"; run
printf "%-46s " "  ... 267 cells"; run s.rk5 meristemT.init

printf "%-46s " "AuxinModelSimple1Wall"
mk "AuxinModelSimple1Wall 7 1 3
0.05
0.03
0.4
1.2
0.1
0.3
0.07
0 1 2"; run
printf "%-46s " "  ... 267 cells"; run s.rk5 meristemT.init

printf "%-46s " "AuxinTransportCellCellNoGeometry"
mk "AuxinTransportCellCellNoGeometry 7 1 4
0.1
0.8
0.2
0.6
0.4
2.0
0.3
0 1 2 3"; run
printf "%-46s " "  ... 267 cells"; run s.rk5 meristemT.init

echo
echo "Storing the membrane PIN (MISMATCH at t=0 only; see NOTES.md):"
printf "%-46s " "AuxinModelSimple1 + PIN store"
mk "AuxinModelSimple1 12 2 4 1
$P12
0 1 2 3
5"; run euler.rk5

echo
echo "The SimpleROPModel family:"
# ROP7 gates on a wall marker, so the fixture needs real 0/1 flags rather than
# seeded noise - without them neither branch fires and the test is vacuous.
python3 - <<'PYEOF'
lines = open("twoSquareT.init").read().split("\n")
out, inwall, seen = [], False, 0
for l in lines:
    t = l.split()
    if len(t) == 3 and t[0] == "7" and t[1] == "1" and t[2] == "10":
        inwall = True; out.append(l); continue
    if inwall and len(t) == 11:
        t[7] = "1" if seen % 2 == 0 else "0"
        seen += 1; out.append(" ".join(t)); continue
    out.append(l)
open("twoSquareM.init", "w").write("\n".join(out))
PYEOF
P13="0.05
0.02
0.3
0.4
0.6
0.15
0.1
0.03
0.5
0.25
0.4
2.0
0.07"
for n in "" 2 3 4 6; do
  printf "%-46s " "SimpleROPModel$n"
  mk "SimpleROPModel$n 13 2 2 2
$P13
0 1
3 1"; run euler.rk5
done
printf "%-46s " "SimpleROPModel5"
mk "SimpleROPModel5 16 2 2 2
$P13
0.2
1.5
0.35
0 1
3 1"; run euler.rk5
printf "%-46s " "SimpleROPModel7 (marked walls)"
mk "SimpleROPModel7 15 2 2 3
$P13
0.12
0.09
0 1
3 1 7"; run euler.rk5 twoSquareM.init

echo
echo "The gradient models:"
P15="0.05
0.02
0.3
0.4
0.6
0.15
0.1
0.03
0.5
0.25
0.4
2.0
0.07
0.12
0.09"
for n in UpInternalGradientModel DownInternalGradientModel UpExternalGradientModel \
         DownInternalGradientModelSingleCell DownInternalGradientModelGeometric; do
  printf "%-46s " "$n"
  mk "$n 15 2 2 3
$P15
0 1
3 1 7"; run euler.rk5 twoSquareM.init
done
