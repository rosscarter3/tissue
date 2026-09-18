#!/bin/zsh
# Comparison harness for legacy/calculate.cc (diagnostic quantities).
#
# Run from a directory holding s.rk5 and the fixtures built with
#   python3 tests/port/cellvar_init.py twoSquareT.init dir2d.init \
#     "1.0 0.0 -0.3 0.9 0.10 2.0 1.00 0.30" \
#     "0.0 1.0  0.6 0.8 0.20 2.0 1.98 0.75"
#   python3 tests/port/cellvar_init.py tests/port/tri3D.init dir3d.init \
#     "1.0 0.0 0.0 -0.3 0.9  0.2 0.10 2.0 1.00 0.30" \
#     "0.0 0.0 1.0  0.6 0.8 -0.2 0.20 2.0 1.98 0.75"
#
# dir3d: vector A at 0-2, vector B at 3-5, spare scalars at 6-9.
#
# Every reaction here writes cellData from derivs() rather than adding to
# cellDerivs, so each differs from legacy at the t=0 print only - legacy does
# one extra derivative evaluation before that first print. That is the known
# class described in README item 5, and it is why the lines below report a
# difference at t=0 and none after.
#
# VertexVelocity and TissueVolumeChange read vertexDerivs, so a model with
# only those in it measures a tissue at rest and reports zero: every case
# pairs them with a wall spring that actually moves vertices. The spring is
# listed first because these read the derivatives accumulated so far.
T=/Users/ross/projects/tissue
mk() { printf '%d 0 0\n\n%s\n' "$1" "$2" > single.model; }
run() { python3 $T/tools/port/compare.py single.model ${2:-dir3d.init} ${1:-s.rk5} --tol=1e-9 2>&1 | tail -1; }
# TissueVolumeChange stores a difference against its own previous value, so it
# reports not the state but the gap between two derivative evaluations. Legacy
# evaluates once more before printing than this build does, so that gap is
# measured over a different pair and the stored value lags by one evaluation.
# Measured: with Euler the offset is exactly linear in h (rel/h constant at
# 2.28 over h from 5e-4 to 4e-3), which is a one-step lag and not an
# arithmetic difference - the other two quantities this reaction stores, which
# are functions of the state alone, match exactly at every print.
lagrun() { python3 $T/tools/port/compare.py single.model ${2:-dir3d.init} ${1:-s.rk5} --tol=1e-8 2>&1 | tail -1; }

SPRING="VertexFromWallSpring 2 1 1
6.0
0.6

0"

printf "%-46s " "Calculate::AngleVectors"
mk 1 "Calculate::AngleVectors 0 2 2 1
0 3

6"; run

printf "%-46s " "  ... alias CalculateAngleVectors"
mk 1 "CalculateAngleVectors 0 2 2 1
0 3

6"; run

printf "%-46s " "  ... with the tissue moving"
mk 2 "$SPRING

Calculate::AngleVectors 0 2 2 1
0 3

6"; run

printf "%-46s " "Calculate::AngleVectorXYplane"
mk 1 "Calculate::AngleVectorXYplane 0 2 1 1
0

6"; run

printf "%-46s " "  ... vector along z (xy degenerate)"
mk 1 "Calculate::AngleVectorXYplane 0 2 1 1
3

6"; run

printf "%-46s " "Calculate::AngleVector (axis 0)"
mk 1 "Calculate::AngleVector 1 2 1 1
0

0
6"; run

printf "%-46s " "  ... alias AngleVector, other vector"
mk 1 "AngleVector 1 2 1 1
0

3
6"; run

printf "%-46s " "Calculate::VertexVelocity (at rest)"
mk 1 "Calculate::VertexVelocity 0 1 1
6"; run

printf "%-46s " "Calculate::VertexVelocity (moving)"
mk 2 "$SPRING

Calculate::VertexVelocity 0 1 1
6"; run

printf "%-46s " "  ... 2D"
mk 2 "$SPRING

Calculate::VertexVelocity 0 1 1
4"; run s.rk5 dir2d.init

printf "%-46s " "  ... velocity gating an MT update"
mk 3 "$SPRING

Calculate::VertexVelocity 0 1 1
6

UpdateMTDirectionEquilibrium 2 3 1 1 1
0.9
0.5

0
3
6"; run

printf "%-46s " "TissueVolumeChange (moving) [lag 1e-8]"
mk 2 "$SPRING

Calculate::TissueVolumeChange 0 1 6
0 6 0 7 0 8"; lagrun

printf "%-46s " "  ... alias TemplateVolumeChange [lag 1e-8]"
mk 2 "$SPRING

TemplateVolumeChange 0 1 6
0 6 1 6 0 7"; lagrun

printf "%-46s " "  ... at rest"
mk 1 "Calculate::TissueVolumeChange 0 1 6
0 6 0 7 0 8"; run
