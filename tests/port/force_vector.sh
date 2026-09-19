#!/bin/zsh
# Comparison harness for Force::VectorLinear (legacy force.cc).
#
# Run from a directory holding s.rk5, twoSquareT.init (2D) and dir3d.init (3D).
#
# The force is applied to a named list of vertices and ramped over deltaT, so
# the cases below vary the number of components, the ramp against the run
# length, and which vertices are pulled. A wall spring is included so the
# pulled vertices have something to pull against - without it the whole
# template just translates and several configurations become
# indistinguishable.
T=/Users/ross/projects/tissue
SPRING="VertexFromWallSpring 2 1 1
6.0
0.6

0"
mk() { printf '2 0 0\n\n%s\n\n%s\n' "$SPRING" "$1" > single.model; }
printf 'Euler\n0 2\n0 2\n0.001\n' > force_euler.rk5
run() { python3 $T/tools/port/compare.py single.model ${2:-twoSquareT.init} ${1:-force_euler.rk5} --tol=1e-12 2>&1 | tail -1; }

# deltaT against the solver's end time of 2: 0.5 finishes the ramp early, 8
# is still climbing at the end.
for dt in 0.5 8; do
  printf "%-46s " "2 params (F_x, deltaT=$dt), 2D"
  mk "Force::VectorLinear 2 1 2
0.4
$dt

0 3"; run
done

printf "%-46s " "3 params (F_x F_y deltaT), 2D"
mk "Force::VectorLinear 3 1 2
0.4
-0.2
1.0

0 3"; run

printf "%-46s " "  ... alias VertexFromForceLinear"
mk "VertexFromForceLinear 3 1 2
0.4
-0.2
1.0

0 3"; run

printf "%-46s " "  ... a different vertex list"
mk "Force::VectorLinear 3 1 3
0.4
-0.2
1.0

1 2 4"; run

printf "%-46s " "  ... one vertex only"
mk "Force::VectorLinear 3 1 1
0.4
-0.2
1.0

2"; run

# In 3D. The 4-parameter form is the only one whose component count matches
# the tissue dimension, so it is also the only one legacy does not feed its
# deltaT in as an extra force component - see the note in force_vector.cpp.
printf "%-46s " "4 params (F_x F_y F_z deltaT), 3D"
mk "Force::VectorLinear 4 1 2
0.4
-0.2
0.3
1.0

0 3"; run force_euler.rk5 dir3d.init

printf "%-46s " "  ... 3 params in 3D (deltaT leaks into F_z)"
mk "Force::VectorLinear 3 1 2
0.4
-0.2
1.0

0 3"; run force_euler.rk5 dir3d.init

printf "%-46s " "  ... 2 params in 3D (deltaT leaks into F_y)"
mk "Force::VectorLinear 2 1 2
0.4
1.0

0 3"; run force_euler.rk5 dir3d.init
