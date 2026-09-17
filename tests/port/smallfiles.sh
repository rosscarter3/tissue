#!/bin/zsh
# Comparison harness for the twelve classes finishing off pressure2D.cc,
# bending.cc, sisterVertex.cc, turgorGrowth.cc, cellTime.cc and dilution.cc.
#
# Run from a directory holding twoSquare.init, meristem.init, euler.rk5 and
# s.rk5. Three reactions are expected to MISMATCH - they fix wrong-index bugs
# in legacy; see NOTES.md and tests/port/bending_refcheck.py.
T=/Users/ross/projects/tissue
mk() { printf '%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-twoSquare.init} $1 --tol=1e-9 2>&1 | tail -1; }

printf "%-52s " "CellTimeDerivative"
mk "1 0 0

CellTimeDerivative 0 1 1
5"; run s.rk5

printf "%-52s " "Dilution::FromVertexDerivs"
mk "2 0 0

Pressure2D::AreaPotential 2 0
0.5
0

DilutionFromVertexDerivs 0 1 2
1 2"; run euler.rk5

printf "%-52s " "WaterVolumeFromTurgor"
mk "1 0 0

WaterVolumeFromTurgor 5 2 1 1
0.1
0.5
2.0
0
0
1
5"; run s.rk5

printf "%-52s " "WaterVolumeFromTurgor (denyShrink, negative P)"
mk "1 0 0

WaterVolumeFromTurgor 5 2 1 1
0.4
0.05
2.0
1
1
1
5"; run s.rk5

printf "%-52s " "Pressure2D::AreaPotentialTri (V-normalised)"
mk "1 0 0

Pressure2D::AreaPotentialTri 2 0
0.5
1"; run euler.rk5

printf "%-52s " "Pressure2D::AreaPotentialTri (internal only, 267)"
mk "1 0 0

Pressure2D::AreaPotentialTri 3 0
0.5
0
1"; run euler.rk5 meristem.init

printf "%-52s " "Pressure2D::AreaPotentialTriSpatialThreshold (267)"
mk "1 0 0

Pressure2D::AreaPotentialTriSpatialThreshold 2 1 1
0.5
3.0
1"; run euler.rk5 meristem.init

printf "%-52s " "Pressure2D::AreaPotentialTargetArea (inflating)"
mk "2 0 0

Creation::Zero 1 1 1
0.5
3

Pressure2D::AreaPotentialTargetArea 2 1 1
0.5
0
3"; run euler.rk5

printf "%-52s " "Pressure2D::AreaPotentialTargetArea (contracting)"
mk "2 0 0

Degradation::One 1 1 1
0.5
3

Pressure2D::AreaPotentialTargetArea 2 1 1
0.5
0
3"; run euler.rk5

printf "%-52s " "Pressure2D::AreaPotentialTargetArea (no contraction)"
mk "2 0 0

Degradation::One 1 1 1
0.5
3

Pressure2D::AreaPotentialTargetArea 2 1 1
0.5
1
3"; run euler.rk5

printf "%-52s " "Bending::AngleInitiate + AngleRelax (267)"
mk "3 0 0

Bending::AngleInitiate 0 1 1
1

Pressure2D::AreaPotential 2 0
0.5
0

Bending::AngleRelax 1 1 1
0.5
1"; run euler.rk5 meristem.init

printf "%-52s " "Bending::Angle (fixed; MISMATCH expected)"
mk "1 0 0

Bending::Angle 1 1 1
0.5
1"; run euler.rk5

printf "%-52s " "Bending::NeighborCenter (fixed; MISMATCH expected)"
mk "1 0 0

Bending::NeighborCenter 1 1 1
0.5
0"; run euler.rk5

printf '1\n0 5\n' > sister
printf "%-52s " "SisterVertex::InitiateFromFile + SpringCellConc"
mk "2 0 0

SisterVertex::InitiateFromFile 0 0

SisterVertex::SpringCellConc 1 1 1
0.5
9"; run euler.rk5

printf '2\n0 5\n1 4\n' > sister
printf "%-52s " "SpringCellConc with BreakLength"
mk "2 0 0

SisterVertex::InitiateFromFile 0 0

Pressure2D::AreaPotential 2 0
2.0
0

SisterVertex::SpringCellConc 2 1 1
0.5
1.1
9"; run euler.rk5

printf '1\n5 0\n' > sister
printf "%-52s " "SpringCellConc reversed pair (MISMATCH expected)"
mk "2 0 0

SisterVertex::InitiateFromFile 0 0

SisterVertex::SpringCellConc 1 1 1
0.5
9"; run euler.rk5

echo
echo "Independent check of the two reactions that cannot be compared:"
python3 $T/tests/port/bending_refcheck.py
