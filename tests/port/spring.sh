#!/bin/zsh
# Comparison harness for the wall-spring family in legacy/mechanicalSpring.cc.
#
# Run from a directory holding euler.rk5 and seeded inits built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
# Fixture layout: wall variable 1 = a per-wall signal, 5 = a spare to save the
# force into; cell variable 1 = the concentration the Hill variant reads.
#
# Any form that saves the force into a wall variable is expected to MISMATCH
# at t=0 only: that write happens inside derivs, so it sees legacy's extra
# pre-print derivative evaluation (README item 5). The trajectory and every
# later print agree exactly. See NOTES.md.
T=/Users/ross/projects/tissue
mk() { printf '2 0 0\n\nWallGrowth::Constant 2 1 1\n0.05\n1\n0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${1:-twoSquareT.init} euler.rk5 --tol=1e-9 2>&1 | tail -1; }

printf "%-52s " "WallMechanics::Spring"
mk "WallMechanics::Spring 2 1 1
0.4
0.8
0"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

printf "%-52s " "WallMechanics::SpringEpidermal"
mk "WallMechanics::SpringEpidermal 2 1 1
0.4
0.8
0"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

printf "%-52s " "WallMechanics::SpringEpidermalCell"
mk "WallMechanics::SpringEpidermalCell 2 1 1
0.4
0.8
0"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

printf "%-52s " "WallMechanics::SpringConcentrationHill"
mk "WallMechanics::SpringConcentrationHill 5 1 2
0.2
0.9
0.5
2.0
0.8
0 1"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

printf "%-52s " "VertexFromWallBoundarySpring"
mk "VertexFromWallBoundarySpring 2 1 1
0.4
0.8
0"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

# Fibre (microtubule) direction read from cell variables: components at the
# given index, then a flag. Index 2 has the flag set in both fixture cells,
# index 5 in only one - the unoriented path matters, and differs between these
# two reactions.
printf "%-52s " "VertexFromWallSpringMT (both cells oriented)"
mk "VertexFromWallSpringMT 3 1 2
0.2
0.6
0.8
0 2"; run
printf "%-52s " "  ... one cell unoriented"
mk "VertexFromWallSpringMT 3 1 2
0.2
0.6
0.8
0 5"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init

printf "%-52s " "VertexFromWallSpringMTConcentrationHill"
mk "VertexFromWallSpringMTConcentrationHill 6 2 2 1
1.0
0.7
0.4
0.5
2.0
0.8
0 5
1"; run
printf "%-52s " "  ... 267 cells"
run meristemT.init
printf "%-52s " "  ... both cells oriented"
mk "VertexFromWallSpringMTConcentrationHill 6 2 2 1
1.0
0.7
0.4
0.5
2.0
0.8
0 2
1"; run

echo
echo "Force-save forms (MISMATCH at t=0 only; see NOTES.md):"
printf "%-52s " "WallMechanics::Spring + K_force2 + save"
mk "WallMechanics::Spring 3 3 1 1 1
0.4
0.8
1.2
0
5
1"; run
printf "%-52s " "WallMechanics::SpringEpidermal + save"
mk "WallMechanics::SpringEpidermal 2 2 1 1
0.4
0.8
0
5"; run
