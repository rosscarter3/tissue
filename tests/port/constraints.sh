#!/bin/zsh
# Comparison harness for the VertexNoUpdate* clamps in legacy/adhocReaction.cc.
#
# Run from a directory holding s.rk5 and the fixtures built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
#   python3 tests/port/cellvar_init.py tests/port/tri3D.init dir3d.init \
#     "1.0 0.0 0.0 -0.3 0.9  0.2 0.10 2.0 1.00 0.30" \
#     "0.0 0.0 1.0  0.6 0.8 -0.2 0.20 2.0 1.98 0.75"
#
# A clamp only zeroes derivatives, so on its own it does nothing observable:
# every case pairs it with a wall spring, listed first because these reactions
# cancel what has accumulated so far.
#
# The spring is pointed at wall variable 1, not 0. Variable 0 holds each
# wall's true resting length, so a spring reading it finds the seeded mesh
# already at equilibrium, exerts no force, and nothing moves - every clamp
# then "matches" legacy while clamping a tissue that was not going anywhere.
# Variable 1 holds a seeded value unrelated to the geometry, so the springs
# are out of equilibrium and the vertices actually move. Verified by
# comparing each clamped run against the unclamped one; see NOTES.md.
T=/Users/ross/projects/tissue
mk() { printf '%d 0 0\n\n%s\n' "$1" "$2" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-twoSquareT.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }
# One case below sits on a print-rounding boundary rather than differing in
# behaviour: with this solver's 1e-6 error tolerance a single wall variable
# prints as 0.131243 in one binary and 0.131242 in the other, one unit in the
# last of six significant digits, and no other value in the run differs. The
# same comparison is exact (0.000e+00) at solver tolerances 1e-5, 1e-7 and
# 1e-9, which is what rules out a real divergence.
roundrun() { python3 $T/tools/port/compare.py single.model ${2:-twoSquareT.init} ${1:-s.rk5} --tol=2e-6 2>&1 | tail -1; }

SPRING="VertexFromWallSpring 2 1 1
6.0
0.6

1"

printf "%-46s " "VertexNoUpdateFromIndex (already ported)"
mk 2 "$SPRING

VertexNoUpdateFromIndex 0 1 2
0 3"; run

for axis in X Y; do
  printf "%-46s " "VertexNoUpdateFromIndexHold$axis"
  mk 2 "$SPRING

VertexNoUpdateFromIndexHold$axis 0 1 2
0 3"; run
done

printf "%-46s " "  ... HoldZ (3D)"
mk 2 "$SPRING

VertexNoUpdateFromIndexHoldZ 0 1 2
0 3"; run s.rk5 dir3d.init

# The threshold splits the two-square mesh: vertices sit at x = 0, 1 and 2.
for dir in 1 -1; do
  for thresh in 0.5 1.5; do
    printf "%-46s " "VertexNoUpdateFromPosition dir=$dir t=$thresh"
    mk 2 "$SPRING

VertexNoUpdateFromPosition 2 1 1
$thresh
$dir

0"
    # dir=1 t=1.5 is the print-rounding case described above.
    if [ "$dir" = "1" ] && [ "$thresh" = "1.5" ]; then roundrun; else run; fi
  done
done

printf "%-46s " "  ... on y instead of x"
mk 2 "$SPRING

VertexNoUpdateFromPosition 2 1 1
0.5
1

1"; run

printf "%-46s " "VertexNoUpdateFromList"
mk 2 "$SPRING

VertexNoUpdateFromList 0 0"; run

printf "%-46s " "  ... 3D"
mk 2 "$SPRING

VertexNoUpdateFromList 0 0"; run s.rk5 dir3d.init

printf "%-46s " "VertexNoUpdateBoundary (all directions)"
mk 2 "$SPRING

VertexNoUpdateBoundary 0 0"; run

for ax in 0 1; do
  printf "%-46s " "VertexNoUpdateBoundary (axis $ax only)"
  mk 2 "$SPRING

VertexNoUpdateBoundary 0 1 1
$ax"; run
done

printf "%-46s " "  ... two axes"
mk 2 "$SPRING

VertexNoUpdateBoundary 0 1 2
0 1"; run

printf "%-46s " "  ... 3D, all directions"
mk 2 "$SPRING

VertexNoUpdateBoundary 0 0"; run s.rk5 dir3d.init
