#!/bin/zsh
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model twoSquare.init s.rk5 --tol=1e-9 2>&1 | tail -1; }
printf "%-26s " "SpatialCylinder";   mk "Creation::SpatialCylinder 4 1 1
0.5
1.2
2.0
1.0
3"; run
printf "%-26s " "SpatialRing";       mk "Creation::SpatialRing 5 1 1
0.5
1.2
0.8
2.0
1.0
3"; run
printf "%-26s " "SpatialCoordinate"; mk "Creation::SpatialCoordinate 4 2 1 1
0.5
1.2
2.0
1.0
3
0"; run
printf "%-26s " "SpatialPlane";      mk "Creation::SpatialPlane 3 2 1 1
0.4
0.5
1.0
3
0"; run
printf "%-26s " "FromList";          mk "Creation::FromList 1 2 1 2
0.6
3
0 1"; run
printf "%-26s " "OneGeometric";      mk "Creation::OneGeometric 1 2 1 1
0.5
3
1"; run
printf "%-26s " "Sinus";             mk "Creation::Sinus 3 1 1
0.5
2.0
0.25
3"; run
