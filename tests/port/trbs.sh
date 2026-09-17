#!/bin/zsh
# Comparison harness for the TRBS elasticity and the center-triangulation
# growth rules that sit on it.
#
# Run from a directory holding twoSquare3D.init (examples/tutorials/mesh),
# tests/port/tri3D.init, euler3.rk5 and rk3.rk5:
#   printf 'Euler\n0 2\n0 6\n0.001\n'          > euler3.rk5
#   printf 'RK5Adaptive\n0 2\n0 6\n0.1 1e-6\n' > rk3.rk5
#
# Anything that also drives *wall* variables must be compared under Euler:
# legacy's RK5 mis-integrates those (README item 1).
T=/Users/ross/projects/tissue
run() { python3 $T/tools/port/compare.py $1 $2 $3 --tol=1e-9 2>&1 | tail -1; }

cat > _trbs.model <<'M'
4 0 0

CenterTriangulation::Initiate 0 1 1
10

WallGrowth::Constant 2 1 1
0.1
1
0

CenterTriangulation::WallGrowth::Constant 2 1 1
0.1
1
10

VertexFromTRBScenterTriangulation 2 2 1 1
10.0
0.3
0
10
M
printf "%-52s " "VertexFromTRBScenterTriangulation"
run _trbs.model twoSquare3D.init euler3.rk5

# Same without wall dynamics, so RK5Adaptive is comparable too.
cat > _trbs_nowall.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
10

CenterTriangulation::WallGrowth::Constant 2 1 1
0.1
1
10

VertexFromTRBScenterTriangulation 2 2 1 1
10.0
0.3
0
10
M
printf "%-52s " "  ... no wall dynamics, RK5Adaptive"
run _trbs_nowall.model twoSquare3D.init rk3.rk5

cat > _trbs_hill.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
10

CenterTriangulation::WallGrowth::Constant 2 1 1
0.1
1
10

VertexFromTRBScenterTriangulationConcentrationHill 5 2 2 1
1.0
20.0
0.3
0.5
2.0
0 9
10
M
printf "%-52s " "VertexFromTRBScenterTriangulationConcentrationHill"
run _trbs_hill.model twoSquare3D.init euler3.rk5
printf "%-52s " "  ... under RK5Adaptive"
run _trbs_hill.model twoSquare3D.init rk3.rk5

# VertexFromTRBS is defined only for triangular cells.
cat > _trbs_flat.model <<'M'
2 0 0

WallGrowth::Constant 2 1 1
0.1
1
0

VertexFromTRBS 2 1 1
10.0
0.3
0
M
printf "%-52s " "VertexFromTRBS (triangular cells)"
run _trbs_flat.model tri3D.init euler3.rk5
printf '1 0 0\n\nVertexFromTRBS 2 1 1\n10.0\n0.3\n0\n' > _trbs_flat2.model
printf "%-52s " "  ... no wall dynamics, RK5Adaptive"
run _trbs_flat2.model tri3D.init rk3.rk5

for flag in 1 0; do
cat > _ct_stress.model <<M
3 0 0

CenterTriangulation::Initiate 0 1 1
10

VertexFromTRBScenterTriangulation 2 2 1 1
10.0
0.3
0
10

WallGrowth::CenterTriangulation::Stress 4 1 1
0.5
0.01
1
$flag
10
M
printf "%-52s " "CenterTriangulation::WallGrowth::Stress (linear=$flag)"
run _ct_stress.model twoSquare3D.init euler3.rk5
done

cat > _ct_hill.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
10

VertexFromTRBScenterTriangulation 2 2 1 1
10.0
0.3
0
10

WallGrowth::CenterTriangulation::StressConcentrationHill 7 1 2
0.1
0.9
0.5
2.0
0.01
1
1
10 9
M
printf "%-52s " "CenterTriangulation::...::StressConcentrationHill"
run _ct_hill.model twoSquare3D.init euler3.rk5
printf "%-52s " "  ... under RK5Adaptive"
run _ct_hill.model twoSquare3D.init rk3.rk5
