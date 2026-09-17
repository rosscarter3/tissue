#!/bin/zsh
# Comparison harness for legacy/membraneCycling.cc and membraneCyclingAll.cc.
#
# Run from a directory holding euler.rk5 and a seeded init built with
#   python3 tests/port/seed_init.py <tutorial>.init twoSquareT.init 10 8
# Stock inits leave wall variables at zero, and the carrier these reactions
# cycle lives in paired wall variables.
#
# Fixture layout: wall 1,2 = carrier pair; 3,4 = signal pair.
# Cell 0 = auxin, 1 = carrier pool.
#
# The last four are expected to MISMATCH - they fix typos that made legacy
# depend on which of a wall's cells the init listed first (README item 12).
# membrane_swapcheck.py tests them on that symmetry instead.
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model twoSquareT.init euler.rk5 --tol=1e-9 2>&1 | tail -1; }
swap() { python3 $T/tests/port/membrane_swapcheck.py single.model twoSquareT.init euler.rk5; }

CONSTANT="MembraneCycling::Constant 2 2 1 1
0.3
0.4
1
1"
ALLCONSTANT="MembraneCyclingAll::Constant 2 2 1 1
0.3
0.4
1
1"
CROSS="MembraneCycling::CrossMembraneNonLinear 4 2 1 1
0.3
0.4
0.5
2.0
1
1"
LWFNL="MembraneCycling::LocalWallFeedbackNonLinear 4 2 1 2
0.3
0.4
0.5
2.0
1
3 1"
ALLLWFNL="MembraneCyclingAll::LocalWallFeedbackNonLinear 4 2 1 2
0.3
0.4
0.5
2.0
1
3 1"
UTGNL="MembraneCycling::CellUpTheGradientNonLinear 4 2 2 1
0.3
0.4
0.5
2.0
0 1
1"
UTGL="MembraneCycling::CellUpTheGradientLinear 2 2 2 1
0.3
0.4
0 1
1"
INTNL="MembraneCycling::InternalCellNonLinear 4 2 2 1
0.3
0.4
0.5
2.0
0 1
1"
INTL="MembraneCycling::InternalCellLinear 2 2 2 1
0.3
0.4
0 1
1"
FLUX="MembraneCycling::CellFluxExocytosis 3 2 2 1
0.2
0.5
0.6
0 1
1"
LWFL="MembraneCycling::LocalWallFeedbackLinear 2 2 1 2
0.3
0.4
1
3 1"
PINL="MembraneCycling::PINFeedbackLinear 2 2 1 1
0.3
0.4
1
1"
PINNL="MembraneCycling::PINFeedbackNonLinear 4 2 1 1
0.3
0.4
0.5
2.0
1
1"
INHIB="MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition 4 2 1 2
0.3
0.4
0.5
2.0
1
3 1"

for pair in "Constant|$CONSTANT" "All::Constant|$ALLCONSTANT" \
            "CrossMembraneNonLinear|$CROSS" "LocalWallFeedbackNonLinear|$LWFNL" \
            "All::LocalWallFeedbackNonLinear|$ALLLWFNL" \
            "CellUpTheGradientNonLinear|$UTGNL" "CellUpTheGradientLinear|$UTGL" \
            "InternalCellNonLinear|$INTNL" "InternalCellLinear|$INTL" \
            "CellFluxExocytosis|$FLUX"; do
  printf "%-42s " "${pair%%|*}"; mk "${pair#*|}"; run
done

echo
echo "Fixed (expect MISMATCH against legacy, invariance under wall flip):"
for pair in "LocalWallFeedbackLinear|$LWFL" "PINFeedbackLinear|$PINL" \
            "PINFeedbackNonLinear|$PINNL" \
            "All::LWFNonLinearInhibition|$INHIB"; do
  printf "%-42s " "${pair%%|*}"; mk "${pair#*|}"; run; swap
done

echo
echo "Control - the swap check on reactions that needed no fix:"
printf "%-42s\n" "CrossMembraneNonLinear"; mk "$CROSS"; swap
printf "%-42s\n" "LocalWallFeedbackNonLinear"; mk "$LWFNL"; swap
