#!/bin/zsh
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model twoSquare.init "$2" --tol=1e-9 2>&1 | tail -1; }
printf "%-28s " "OneToTwo";          mk "MassAction::OneToTwo 1 2 3 0
0.4
1 2 3"; run x s.rk5
printf "%-28s " "HillSimple";        mk "MassAction::HillSimple 3 2 1 1
0.4
0.8
2.0
1
2"; run x s.rk5
printf "%-28s " "GeneralEnzymatic";  mk "MassAction::GeneralEnzymatic 1 3 1 1 1
0.3
1
2
3"; run x s.rk5
printf "%-28s " "GeneralWall (Euler)";  mk "MassAction::GeneralWall 1 2 1 1
0.3
1
2"; run x euler.rk5
printf "%-28s " "OneToTwoWall (Euler)"; mk "MassAction::OneToTwoWall 1 2 0 3
0.3
0 1 2"; run x euler.rk5
printf "%-28s " "TwoToOneWall (Euler)"; mk "MassAction::TwoToOneWall 1 2 0 3
0.3
0 1 2"; run x euler.rk5
