#!/bin/zsh
# Comparison harness for legacy/boolean.cc (flag gates and counters).
#
# Run from a directory holding s.rk5 and the fixture built with
#   python3 tests/port/cellvar_init.py twoSquareT.init bool.init \
#     "1 1 0 1 0 0.7 0.2 0" \
#     "1 0 0 0 1 0.3 0.9 0"
#
# Cell variables 0-4 are exact 0/1 flags, 5 and 6 are concentrations for the
# threshold gates, 7 is a spare output. The two cells differ in every input,
# so a gate that ignores one of its inputs shows up.
#
# These compare against 1 and 0 with == on doubles, so the flags must be
# exactly 0 or 1; seeded values in (0.05, 0.95) would make every gate false
# and every case pass while testing nothing.
T=/Users/ross/projects/tissue
mk() { printf '1 0 0\n\n%s\n' "$1" > single.model; }
run() { python3 $T/tools/port/compare.py single.model bool.init s.rk5 --tol=1e-12 2>&1 | tail -1; }

# gatetype 0 clears the output when the condition is false; anything else
# latches. Both are covered for the two gates that offer the choice.
for gt in 0 1; do
  printf "%-46s " "Boolean::AndGate gatetype=$gt"
  mk "Boolean::AndGate 1 2 2 1
$gt
0 1
7"; run
  printf "%-46s " "Boolean::AndNotGate gatetype=$gt"
  mk "Boolean::AndNotGate 1 2 2 1
$gt
0 1
7"; run
done

printf "%-46s " "  ... AndGate alias, output over a live flag"
mk "AndGate 1 2 2 1
0
0 1
3"; run

# The two three-input gates get the *same* indices, so the only thing between
# them is the condition (!input2 against input2) and the clearing. With
# different indices they happen to agree on this fixture, which would make the
# pair look interchangeable.
printf "%-46s " "Boolean::AndSpecialGate (always clears)"
mk "Boolean::AndSpecialGate 0 2 3 1
0 1 2
7"; run

printf "%-46s " "Boolean::AndSpecialGate2 (latches)"
mk "Boolean::AndSpecialGate2 0 2 3 1
0 1 2
7"; run

for thr in 0.5 0.9; do
  printf "%-46s " "Boolean::AndSpecialGate3 threshold=$thr"
  mk "Boolean::AndSpecialGate3 1 2 3 1
$thr
5 1 2
7"; run
done

for t1 in 0.5 0.9; do
  printf "%-46s " "Boolean::AndThresholdGate t1=$t1"
  mk "Boolean::AndThresholdGate 2 2 2 1
$t1
0.1
5 6
7"; run
done

printf "%-46s " "  ... alias AndThresholdsGate (with the s)"
mk "AndThresholdsGate 2 2 2 1
0.5
0.1
5 6
7"; run

# The counters add one per update() - per solver step, not per unit time.
printf "%-46s " "Boolean::Count"
mk "Boolean::Count 0 1 1
7"; run

printf "%-46s " "Boolean::FlagCount"
mk "Boolean::FlagCount 0 2 1 1
1
7"; run

printf "%-46s " "Boolean::AndGateCount"
mk "Boolean::AndGateCount 0 2 2 1
0 1
7"; run

printf "%-46s " "Boolean::OrGateCount"
mk "Boolean::OrGateCount 0 2 2 1
1 4
7"; run

printf "%-46s " "Boolean::OrSpecialGateCount"
mk "Boolean::OrSpecialGateCount 0 2 2 1
2 6
7"; run

# --- threshold and flag bookkeeping (adhocReaction.cc) ------------------------
# Same fixture: flags in 0-4, concentrations in 5-6, output in 7.
# resetFlag only decides what happens to cells *below* the threshold, and only
# 0 writes anything there. Writing into variable 7, which starts at zero, both
# settings leave a zero behind and the parameter looks inert. Variable 4 is 1
# in the below-threshold cell, so the two settings are distinguishable.
for rf in 0 1; do
  printf "%-46s " "ThresholdSwitch resetFlag=$rf"
  mk "ThresholdSwitch 2 2 1 1
0.5
$rf
5
4"; run
  printf "%-46s " "ThresholdReset resetFlag=$rf"
  mk "ThresholdReset 2 2 1 1
0.5
$rf
5
4"; run
done

printf "%-46s " "  ... ThresholdSwitch over a live flag"
mk "ThresholdSwitch 2 2 1 1
0.5
0
6
1"; run

printf "%-46s " "FlagAddValue"
mk "FlagAddValue 1 2 1 1
0.25
1
7"; run

printf "%-46s " "  ... negative value, different flag"
mk "FlagAddValue 1 2 1 1
-2.5
4
7"; run

printf "%-46s " "CopyVariable"
mk "CopyVariable 0 2 1 1
5
7"; run

printf "%-46s " "  ... copying a flag onto a concentration"
mk "CopyVariable 0 2 1 1
0
6"; run
