#!/bin/zsh
# Comparison harness for the centre-triangulation pressure pair in
# legacy/mechanical.cc.
#
# Run from a directory holding s.rk5 and a copy of the mesh tutorial's
# twoSquare3D.init. These reactions read the cell centre from the
# centre-triangulation cell variables, so every model here runs
# CenterTriangulation::Initiate first to create them - without it the centre
# reads as whatever happens to be in those columns and the force is nonsense.
#
# 3D only, hence a 3D init with 10 cell variables, so the CT block starts at
# index 10. Build it with
#   python3 tests/port/cellvar_init.py \
#     examples/tutorials/mesh/twoSquare3D.init ct3d.init \
#     "0.0 1.0 1.0 0.4 0.0 0 0 0 0 1" \
#     "0.0 1.0 1.0 0.9 0.0 0 0 0 0 0"
# rather than using twoSquare3D.init directly: there, cell variable 3 is 1.0
# in both cells, so the concentration-scaled case multiplies by one and is
# indistinguishable from the constant one. 0.4 and 0.9 make it bite.
#
# Two cases need a tighter solver than s.rk5 (error tolerance 1e-6), and both
# for the same reason rather than any difference in the port: at 1e-6 the
# answer is only good to about 1e-6, so once the two binaries' adaptive step
# sequences part company the results differ by roughly that much. Tightening
# to 1e-9 makes them agree to 1e-12, which is what shows the arithmetic is
# identical. A local tight.rk5 is written below for those.
T=/Users/ross/projects/tissue
mk() { printf '%d 0 0\n\nCenterTriangulation::Initiate 0 1 1\n10\n\n%s\n' "$1" "$2" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-ct3d.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }
printf 'RK5Adaptive\n0 2\n0 2\n0.1 1e-9\n' > tight.rk5
tightrun() { python3 $T/tools/port/compare.py single.model ct3d.init tight.rk5 --tol=1e-12 2>&1 | tail -1; }

# --- constant pressure -------------------------------------------------------
# Concentration index 0 means "no concentration" to legacy, so this is the
# plain constant-pressure form.
printf "%-50s " "VertexFromCellPressure (no conc, no vol norm)"
mk 2 "CenterTriangulation::VertexFromCellPressure 2 1 2
0.8
0

10 0"; run

printf "%-50s " "  ... volume normalised"
mk 2 "CenterTriangulation::VertexFromCellPressure 2 1 2
0.8
1

10 0"; run

printf "%-50s " "  ... scaled by cell variable 3"
mk 2 "CenterTriangulation::VertexFromCellPressure 2 1 2
0.8
0

10 3"; run

printf "%-50s " "  ... scaled and volume normalised"
mk 2 "CenterTriangulation::VertexFromCellPressure 2 1 2
0.8
1

10 3"; run

printf "%-50s " "  ... alias VertexFromCellPressurecenterTriang."
mk 2 "VertexFromCellPressurecenterTriangulation 2 1 2
0.8
0

10 0"; run

# Deflation is the ill-conditioned direction: the cell shrinks towards its own
# centre, so the centre-to-midpoint vector this force is built from shrinks
# towards zero and normalising it amplifies everything. At -0.8 the cell has
# collapsed by t=2 and v2 needs 176k accepted and 63k *rejected* steps to get
# there; the two binaries agree exactly (1e-12) up to t=1 and then part. -0.1
# keeps the geometry valid for the whole run and still tests the sign.
printf "%-50s " "  ... negative K (deflating)"
mk 2 "CenterTriangulation::VertexFromCellPressure 2 1 2
-0.1
0

10 0"; run

# --- linear ramp -------------------------------------------------------------
# deltaT against the solver's end time of 2 decides whether the ramp finishes:
# 0.5 completes early, 4 is still climbing at the end, so both the ramping and
# the saturated regimes are covered.
for dt in 0.5 4; do
  printf "%-50s " "VertexFromCellPressureLinear deltaT=$dt"
  mk 2 "CenterTriangulation::VertexFromCellPressureLinear 3 1 1
0.8
0
$dt

10"; run
done

printf "%-50s " "  ... volume normalised, deltaT=1"
mk 2 "CenterTriangulation::VertexFromCellPressureLinear 3 1 1
0.8
1
1

10"; run

printf "%-50s " "  ... alias ...centerTriangulationLinear"
mk 2 "VertexFromCellPressurecenterTriangulationLinear 3 1 1
0.8
0
1

10"; run

# With a wall spring resisting it, so the tissue reaches a balance rather than
# expanding freely - closer to how these are used.
printf "%-50s " "  ... against a wall spring"
printf '3 0 0\n\nCenterTriangulation::Initiate 0 1 1\n10\n\nVertexFromWallSpring 2 1 1\n6.0\n0.6\n\n0\n\nCenterTriangulation::VertexFromCellPressureLinear 3 1 1\n0.8\n0\n1\n\n10\n' > single.model; tightrun
