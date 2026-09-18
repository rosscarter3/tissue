#!/bin/zsh
# Comparison harness for legacy/directionReaction.cc (MT direction updates).
#
# Run from a directory holding s.rk5 and the fixtures built with
#   python3 tests/port/cellvar_init.py twoSquareT.init dir2d.init \
#     "1.0 0.0 -0.3 0.9 0.10 2.0 1.00 0.30" \
#     "0.0 1.0  0.6 0.8 0.20 2.0 1.98 0.75"
#   python3 tests/port/cellvar_init.py tests/port/tri3D.init dir3d.init \
#     "1.0 0.0 0.0 -0.3 0.9  0.2 0.10 2.0 1.00 0.30" \
#     "0.0 0.0 1.0  0.6 0.8 -0.2 0.20 2.0 1.98 0.75"
#
# dir2d: target 0-1, MT 2-3, velocity 4, stresses 5-6, concentration 7.
# dir3d: target 0-2, MT 3-5, velocity 6, stresses 7-8, concentration 9.
# Cell 0 is set up to fail the axis-alignment test (its MT points away from
# the target) and to pass both the velocity and stress gates; cell 1 passes
# neither, so a gate that is ignored shows up as a difference.
#
# UpdateMTDirection and UpdateMTDirectionEquilibrium initiate the MT direction
# *to* the target, so on their own they have nothing left to follow and are
# permanent no-ops. Every case for those two therefore pairs them with
# RotatingDirection driving the target, which is also what exercises it.
T=/Users/ross/projects/tissue
mk() { printf '%d 0 0\n\n%s\n' "$1" "$2" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-dir2d.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }

ROT2="RotatingDirection 1 1 1
0.6

0"
ROT3="RotatingDirection 1 1 1
0.6

0"

printf "%-46s " "RotatingDirection (2D)"
mk 1 "$ROT2"; run
printf "%-46s " "RotatingDirection (3D)"
mk 1 "$ROT3"; run s.rk5 dir3d.init

printf "%-46s " "ContinousMTDirection (2D)"
mk 1 "ContinousMTDirection 1 2 1 1
0.8

0
2"; run

printf "%-46s " "ContinousMTDirection + rotating target"
mk 2 "$ROT2

ContinousMTDirection 1 2 1 1
0.8

0
2"; run

printf "%-46s " "ContinousMTDirection3d"
mk 1 "ContinousMTDirection3d 1 2 1 1
0.8

0
3"; run s.rk5 dir3d.init

printf "%-46s " "ContinousMTDirection3d + rotating target"
mk 2 "$ROT3

ContinousMTDirection3d 1 2 1 1
0.8

0
3"; run s.rk5 dir3d.init

printf "%-46s " "UpdateMTDirection (2D, rotating target)"
mk 2 "$ROT2

UpdateMTDirection 1 2 1 1
0.9

0
2"; run

printf "%-46s " "UpdateMTDirection (3D, rotating target)"
mk 2 "$ROT3

UpdateMTDirection 1 2 1 1
0.9

0
3"; run s.rk5 dir3d.init

printf "%-46s " "UpdateMTDirection k_rate=0 (no-op)"
mk 2 "$ROT2

UpdateMTDirection 1 2 1 1
0

0
2"; run

# --- Equilibrium: one, two and three parameters ------------------------------
printf "%-46s " "MTDirEquilibrium 1 param (ungated)"
mk 2 "$ROT2

UpdateMTDirectionEquilibrium 1 2 1 1
0.9

0
2"; run

for thresh in 0.15 0.5; do
  printf "%-46s " "MTDirEquilibrium velocity gate $thresh"
  mk 2 "$ROT2

UpdateMTDirectionEquilibrium 2 3 1 1 1
0.9
$thresh

0
2
4"; run
done

for sthresh in 0.1 1.0; do
  printf "%-46s " "MTDirEquilibrium +stress gate $sthresh"
  mk 2 "$ROT2

UpdateMTDirectionEquilibrium 3 4 1 1 1 2
0.9
0.5
$sthresh

0
2
4
5 6"; run
done

printf "%-46s " "MTDirEquilibrium 3 param (3D)"
mk 2 "$ROT3

UpdateMTDirectionEquilibrium 3 4 1 1 1 2
0.9
0.5
0.1

0
3
6
7 8"; run s.rk5 dir3d.init

printf "%-46s " "MTDirEquilibrium k_rate=0 (no-op)"
mk 2 "$ROT2

UpdateMTDirectionEquilibrium 2 3 1 1 1
0
0.5

0
2
4"; run

# --- ConcenHill --------------------------------------------------------------
# Writes the direction straight into cellData from derivs, so it runs once per
# derivative evaluation rather than once per step; see NOTES.md.
printf "%-46s " "UpdateMTDirectionConcenHill (3D)"
mk 1 "UpdateMTDirectionConcenHill 3 3 1 1 1
0.7
0.4
2.0

0
3
9"; run s.rk5 dir3d.init

printf "%-46s " "  ... with rotating target"
mk 2 "$ROT3

UpdateMTDirectionConcenHill 3 3 1 1 1
0.7
0.4
2.0

0
3
9"; run s.rk5 dir3d.init

printf "%-46s " "UpdateMTDirectionConcenHill k_rate=0"
mk 1 "UpdateMTDirectionConcenHill 3 3 1 1 1
0
0.4
2.0

0
3
9"; run s.rk5 dir3d.init
