#!/bin/zsh
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model twoSquare.init euler.rk5 --tol=1e-9 2>&1 | tail -1; }
printf "%-42s " "StressSpatial";        mk "WallGrowth::StressSpatial 6 2 2 1
0.3
0.01
1.5
2.0
1.0
0.0
0 1
1"; run
printf "%-42s " "StressSpatialSingle";  mk "WallGrowth::StressSpatialSingle 6 2 2 1
0.3
0.01
1.5
2.0
1.0
0.0
0 1
1"; run
printf "%-42s " "StressConcentrationHill"; mk "WallGrowth::StressConcentrationHill 7 2 2 1
0.2
0.5
0.8
2.0
0.01
1.0
0.0
0 6
1"; run
printf "%-42s " "ConstantStressEpidermalAsymmetric"; mk "WallGrowth::ConstantStressEpidermalAsymmetric 2 1 1
0.3
2.0
0"; run
printf "%-42s " "Force";                mk "WallGrowth::Force 2 2 1 1
0.4
0.1
0
1"; run
