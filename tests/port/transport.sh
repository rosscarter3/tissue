#!/bin/zsh
# Comparison harness for legacy/transport.cc.
#
# Run from a directory holding euler.rk5, s.rk5 and seeded inits built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
# Stock inits leave wall variables at zero, which makes every PIN- or
# AUX-driven reaction here transport nothing and pass vacuously.
#
# Wall variable layout in the seeded fixture: 1,2 = PIN pair; 3,4 = wall auxin
# pair; 5,6 = flux pair; 7 = conductivity; 8,9 = membrane species pair.
# Cell variables: 0 = auxin, 1 = AUX, 2 = gate.
#
# The two "+flux" forms are expected to MISMATCH: they write a diagnostic into
# wallData during derivs, so they see legacy's extra derivative evaluations.
# See NOTES.md.
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-twoSquareT.init} $1 --tol=1e-9 2>&1 | tail -1; }

printf "%-46s " "Diffusion::MembraneSimple"
mk "Diffusion::MembraneSimple 1 1 1
0.3
8"; run euler.rk5

printf "%-46s " "Diffusion::SimpleOne"
mk "Diffusion::SimpleOne 1 2 1 1
0.4
0
2"; run s.rk5

printf "%-46s " "Diffusion::ConductiveSimple"
mk "Diffusion::ConductiveSimple 5 2 1 1
0.5
0.2
1.0
0.5
0.1
0
7"; run euler.rk5

printf "%-46s " "Diffusion::2D (static mesh)"
mk "Diffusion::2D 1 1 1
0.3
0"; run s.rk5

printf "%-46s " "Diffusion::2D (267 cells)"
run s.rk5 meristemT.init

printf "%-46s " "ActiveTransportCellEfflux"
mk "ActiveTransportCellEfflux 1 2 1 1
0.5
0
1"; run s.rk5

printf "%-46s " "ActiveTransportCellEffluxMM"
mk "ActiveTransportCellEffluxMM 2 2 1 1
0.5
0.3
0
1"; run s.rk5

printf "%-46s " "DiffusionActiveTransportCell"
mk "DiffusionActiveTransportCell 2 2 1 1
0.2
0.5
0
1"; run s.rk5

printf "%-46s " "DiffusionActiveTransportCell (267 cells)"
run s.rk5 meristemT.init

printf "%-46s " "ActiveTransportWall"
mk "ActiveTransportWall 5 2 2 2
0.1
0.4
0.2
0.6
0.3
0 1
3 1"; run euler.rk5

printf "%-46s " "ActiveTransportWall (267 cells)"
run euler.rk5 meristemT.init

printf "%-46s " "InfluxActiveTransportCell"
mk "InfluxActiveTransportCell 2 2 1 1
0.5
0.2
0
1"; run s.rk5

printf "%-46s " "InfluxActiveTransportCell (267 cells)"
run s.rk5 meristemT.init

echo
echo "Expected mismatches (diagnostic flux field only; see NOTES.md):"
printf "%-46s " "DiffusionActiveTransportCell (+flux)"
mk "DiffusionActiveTransportCell 2 3 1 1 1
0.2
0.5
0
1
5"; run euler.rk5

printf "%-46s " "InfluxActiveTransportCell (+flux)"
mk "InfluxActiveTransportCell 2 3 1 1 1
0.5
0.2
0
1
5"; run euler.rk5

echo
echo "Diffusion::2D cached-position fix (static exact, moving diverges):"
printf "%-46s " "Diffusion::2D + pressure (moving mesh)"
printf '2 0 0\n\nPressure2D::AreaPotential 2 0\n0.6\n0\n\nDiffusion::2D 1 1 1\n0.3\n0\n' > single.model
run euler.rk5
