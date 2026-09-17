#!/bin/zsh
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model twoSquare.init s.rk5 --tol=1e-9 2>&1 | tail -1; }
printf "%-28s " "Hill";          mk "Degradation::Hill 3 2 1 1
0.7
0.4
2.0
0
1"; run
printf "%-28s " "HillN";         mk "Degradation::HillN 5 3 1 1 1
0.5
0.3
2.0
0.6
1.5
2
3
4"; run
printf "%-28s " "TwoGeometric";  mk "Degradation::TwoGeometric 1 2 1 1
0.2
5
1"; run
printf "%-28s " "OneWall";       mk "Degradation::OneWall 1 1 1
0.3
0"; run
printf "%-28s " "OneBoundary";   mk "Degradation::OneBoundary 1 1 1
0.25
6"; run
printf "%-28s " "OneFromList";   mk "Degradation::OneFromList 1 2 1 2
0.4
7
0 1"; run
