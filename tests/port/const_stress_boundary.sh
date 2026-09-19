#!/bin/zsh
# Comparison harness for VertexFromConstStressBoundary (adhocReaction.cc).
#
# Run from a directory holding s.rk5 and stretch3d.init, built with
#   awk 'BEGIN{split("0.60 0.70 0.80 0.90 0.65 0.75 0.85",L," ")}
#        NR>=20 && NR<=26 { sub(/^1\.0/, L[NR-19]) } { print }' \
#     examples/tutorials/mesh/twoSquare3D.init > stretch3d.init
#
# That is the mesh tutorial's 1 x 2 rectangle in the z = 0 plane, so the four
# boundary planes are right_x = 1, left_x = 0, top_y = 2, bottom_y = 0, with
# the walls given *differing* resting lengths against unit-length edges. Every
# vertex lies on at least one side, which is the normal case here: the
# reaction exists to stretch a rectangular template at fixed stress.
#
# Two fixture traps, both of which made this file test nothing at all:
#
# - The reaction *replaces* the velocity of each boundary vertex with its
#   side's mean, so it needs a force to absorb. On the unmodified init the
#   resting lengths equal the edge lengths, the springs are already balanced,
#   and there is nothing to average.
# - One resting length for every wall is not enough either. The template is
#   symmetric, so every vertex on a side then gets the same velocity, the mean
#   equals each of them, and the averaging - the whole point of the reaction -
#   is still a no-op. The seven differing lengths break that symmetry.
# - Pointing the spring at a different wall variable to get that imbalance
#   does not work here: variables 1 onwards are all zero in this init, so the
#   spring divides by a zero resting length and the whole run is NaN. Both
#   binaries produced NaN and every case "matched". The line range above
#   edits only the wall block; a plain sed also rewrites a vertex position on
#   line 13 and quietly moves the template's corner off x = 1.
T=/Users/ross/projects/tissue
SPRING="VertexFromWallSpring 2 1 1
6.0
0.6

0"
mk() { printf '2 0 0\n\n%s\n\n%s\n' "$SPRING" "$1" > single.model; }
# One fixed Euler step rather than an adaptive run. This reaction *assigns*
# velocities, so under RK5Adaptive the two binaries' step sequences part
# company and the comparison reports the solver's own error (~1e-6 at the
# harness's 1e-6 tolerance) instead of anything about the port. One step of
# Euler removes the solver from the question entirely and compares the
# derivative directly, which is what caught the ordering bug below.
printf 'Euler\n0 0.001\n0 2\n0.001\n' > one_step.rk5
run() { python3 $T/tools/port/compare.py single.model stretch3d.init one_step.rk5 --tol=1e-12 2>&1 | tail -1; }

printf "%-46s " "no stress (sides averaged only)"
mk "VertexFromConstStressBoundary 7 0
0
0
0.1
1.0
0.0
2.0
0.0"; run

printf "%-46s " "x stress only"
mk "VertexFromConstStressBoundary 7 0
0.5
0
0.1
1.0
0.0
2.0
0.0"; run

printf "%-46s " "y stress only"
mk "VertexFromConstStressBoundary 7 0
0
0.5
0.1
1.0
0.0
2.0
0.0"; run

printf "%-46s " "both, equal"
mk "VertexFromConstStressBoundary 7 0
0.5
0.5
0.1
1.0
0.0
2.0
0.0"; run

printf "%-46s " "both, unequal"
mk "VertexFromConstStressBoundary 7 0
0.8
0.2
0.1
1.0
0.0
2.0
0.0"; run

printf "%-46s " "compressive (negative stress)"
mk "VertexFromConstStressBoundary 7 0
-0.3
-0.3
0.1
1.0
0.0
2.0
0.0"; run

# A wider sensitivity puts the y = 1 vertices on the top and bottom sides too,
# which changes which vertices each side averages over.
# Sides overlap here: every vertex is within 1.5 of every boundary plane, so
# each vertex belongs to all four. That is the only configuration in which the
# order of the four assignments is observable, and it is what exposed the port
# taking each side's mean *after* writing the previous one instead of taking
# all four first. A rectangular template never does this, which is exactly why
# it needs a test.
printf "%-46s " "wide sensitivity (sides overlap)"
mk "VertexFromConstStressBoundary 7 0
0.5
0.5
1.5
1.0
0.0
2.0
0.0"; run
