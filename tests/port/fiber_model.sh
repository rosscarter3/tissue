#!/bin/zsh
# Comparison harness for legacy/fiberModel.cc (mechanical-property feedback).
#
# Run from a directory holding s.rk5 and the fixtures built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
#   python3 tests/port/cellvar_init.py twoSquareT.init fiber.init \
#     "0.30 10.0 1.0 0.50 0.25 0.10 0.4 0.6" \
#     "0.75 16.0 2.0 0.60 0.40 0.20 0.5 0.7"
#   python3 tests/port/cellvar_init.py twoSquareT.init fiberlow.init \
#     "0.30 10.0 1.0 0.50 0.25 0.10 0.4 0.6" \
#     "0.75  5.0 2.0 0.60 0.40 0.20 0.5 0.7"
#
# fiber.init puts cell 0 inside Deposition's hard-coded [8,14] stress window
# and cell 1 above it; fiberlow.init puts cell 1 below the floor instead.
# Cell variables are: 0 anisotropy, 1 max stress, 2 area, 3 fiber,
# 4 longitudinal fiber, 5 velocity.
#
# These reactions do all their work in update(), so a model with only a
# FiberModel rule never moves a vertex and the comparison is purely of the
# cell-variable table.
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-fiber.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }

# --- General -----------------------------------------------------------------
# Gradual branches: gate on velocity (cell 0 at 0.10, cell 1 at 0.20) and stop
# at Y_M+Y_F = 0.8. A threshold of 0.15 lets only cell 0 through, which is what
# makes the gate observable rather than a no-op.
for flag in 0 1; do
  for thresh in 0.15 0.5; do
    printf "%-46s " "General flag=$flag thresh=$thresh"
    mk "FiberModel 8 3 1 1 1
0.7
$thresh
$flag
0.4
2.0
0.2
0.6
0
0
3
5"; run
  done
done

# Initiation: 1 zeroes the anisotropy and sets Y_L isotropic, 2 starts Y_L from
# the anisotropy already in the cell. k_rate 0 then isolates initiate() from
# any subsequent update.
for init in 1 2; do
  for flag in 0 1; do
    printf "%-46s " "General init=$init flag=$flag (k_rate 0)"
    mk "FiberModel 8 3 1 1 1
0
0.5
$flag
0.4
2.0
0.2
0.6
$init
0
3
5"; run
  done
done

# Direct (flag 2) needs the second level-1 index, which legacy reads and never
# writes; the update is ungated and omits the matrix term.
printf "%-46s " "General flag=2 (direct, 2 indices)"
mk "FiberModel 8 3 1 2 1
0.7
0.5
2
0.4
2.0
0.2
0.6
0
0
3 4
5"; run

printf "%-46s " "General flag=2 init=2"
mk "FiberModel 8 3 1 2 1
0.7
0.5
2
0.4
2.0
0.2
0.6
2
0
3 4
5"; run

# k_rate 0 with no initiation must be a complete no-op.
printf "%-46s " "General k_rate=0 init=0 (no-op)"
mk "FiberModel 8 3 1 1 1
0
0.5
1
0.4
2.0
0.2
0.6
0
0
3
5"; run

# --- Deposition --------------------------------------------------------------
# The velocity gate is tissue-wide: one cell above the threshold stops every
# cell. 0.15 blocks (cell 1 is at 0.20), 0.5 lets the redistribution run.
for thresh in 0.15 0.5; do
  printf "%-46s " "Deposition thresh=$thresh"
  mk "FiberModel::Deposition 7 1 6
0.3
0.1
$thresh
0
0.4
2.0
1.5
0 1 2 3 4 5"; run
done

# init_flag 3 also drives the longitudinal fiber through the Hill function
# rather than setting it to half the fiber.
printf "%-46s " "Deposition init=3 (Hill longitudinal)"
mk "FiberModel::Deposition 7 1 6
0.3
0.1
0.5
3
0.4
2.0
1.5
0 1 2 3 4 5"; run

printf "%-46s " "Deposition init=1 (uniform initiate)"
mk "FiberModel::Deposition 7 1 6
0.3
0.1
0.5
1
0.4
2.0
1.5
0 1 2 3 4 5"; run

# A cell below the floor takes no share, and the cells above it split the whole
# amount - the branch that stays at zero while the others are rescaled.
printf "%-46s " "Deposition, one cell below floor"
mk "FiberModel::Deposition 7 1 6
0.3
0.1
0.5
0
0.4
2.0
1.5
0 1 2 3 4 5"; run s.rk5 fiberlow.init

printf "%-46s " "Deposition k_rate=0 (no-op)"
mk "FiberModel::Deposition 7 1 6
0
0.1
0.5
0
0.4
2.0
1.5
0 1 2 3 4 5"; run
