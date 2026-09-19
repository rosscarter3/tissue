#!/bin/zsh
# Comparison harness for Hypocotyl3D::limitZdis (legacy hypocotyl3D.cc).
#
# Run from a directory holding one_step.rk5 (a single fixed Euler step) and a
# copy of the hypocotyl template:
#   cp <bench>/3Dhypocotyl/Fig4/3Dhypocotyl.init .
#   python3 tests/port/restlength_init.py 3Dhypocotyl.init hypo_relaxed.init
#
# This reaction is hard-coded to that template - cell variables 37 and 38 as
# the layer and boundary labels, and z = -50 as the plane between the
# cylinder's two ends - so there is no smaller fixture that exercises it. The
# shipped template has 770 cells with 39 variables, 34 of them qualifying,
# and vertices on both sides of the plane.
#
# The shipped hypocotyl models cannot be used directly: they set
# doubleEdge_flag=2, which legacy itself now rejects. So this pairs the
# reaction with a wall spring instead, which is all it needs - it replaces
# the z velocity of each end with that end's mean, so it only needs a force
# to average.
#
# The resting lengths have to be perturbed for that to be true. In the shipped
# template they equal the actual lengths, so the springs start balanced,
# nothing moves, and the reaction averages a set of zeros - it "matches"
# legacy while doing nothing at all. restlength_init.py also varies them wall
# to wall, which matters here: with a uniform scale the template is regular
# enough that each end's vertices share a velocity anyway, so the mean equals
# each of them and the averaging is still a no-op.
#
# One fixed Euler step rather than an adaptive run, for the same reason as in
# const_stress_boundary.sh: this reaction assigns velocities, so under an
# adaptive solver the comparison measures the two step sequences rather than
# the port.
T=/Users/ross/projects/tissue
run() { python3 $T/tools/port/compare.py single.model hypo_relaxed.init one_step.rk5 --tol=1e-12 2>&1 | tail -1; }

printf "%-46s " "limitZdis with a wall spring"
printf '2 0 0\n\nVertexFromWallSpring 2 1 1\n6.0\n0.6\n\n0\n\nHypocotyl3D::limitZdis 0 1 2\n0 1\n' > single.model; run

printf "%-46s " "  ... stronger spring"
printf '2 0 0\n\nVertexFromWallSpring 2 1 1\n25.0\n0.6\n\n0\n\nHypocotyl3D::limitZdis 0 1 2\n0 1\n' > single.model; run

printf "%-46s " "  ... indices are ignored (7 9 vs 0 1)"
printf '2 0 0\n\nVertexFromWallSpring 2 1 1\n6.0\n0.6\n\n0\n\nHypocotyl3D::limitZdis 0 1 2\n7 9\n' > single.model; run

printf "%-46s " "  ... with an axial force as well"
printf '3 0 0\n\nVertexFromWallSpring 2 1 1\n6.0\n0.6\n\n0\n\nGrowthForce::Radial 2 0\n0.01\n1\n\nHypocotyl3D::limitZdis 0 1 2\n0 1\n' > single.model; run
