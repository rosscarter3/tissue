#!/bin/zsh
# Comparison harness for the transversely isotropic (microtubule) TRBS.
#
# Run from a directory holding euler3.rk5 and a seeded 3D fixture:
#   printf 'Euler\n0 2\n0 6\n0.001\n' > euler3.rk5
#   python3 tests/port/seed_init.py examples/tutorials/mesh/twoSquare3D.init \
#           twoSquare3D_MT.init 6 60
#
# Level-0 layout used here: wall length 0; MT direction 1-3; strain anisotropy
# 4; stress anisotropy 5; area ratio 6; transverse modulus / iso energy 7;
# anisotropic energy 8; longitudinal modulus 9; Mises stress 10; stress tensor
# 20-25; cell normal 26-28. Cell variables 13 and 40 are left clear because
# legacy reads them by hard-coded number in some MF modes.
T=/Users/ross/projects/tissue
gen() { # $1=MF $2=neighbourWeight $3=wallGrowth $4=planeStress $5=MTflag
       #  $6=Y_fibre (default 2.0)  $7="1" to add the 12th (loosening) index
python3 - "$1" "$2" "$3" "$4" "$5" "${6:-2.0}" "${7:-0}" <<'PYEOF'
import sys
mf, nw, wg, ps, mt, yf, loosen = sys.argv[1:8]
block = "\nWallGrowth::Constant 2 1 1\n0.05\n1\n0\n" if wg == "1" else ""
idx = "0 1 4 5 6 7 8 9 10 20 26" + (" 30" if loosen == "1" else "")
open("_mt.model", "w").write(f"""{4 if wg=='1' else 3} 0 0

CenterTriangulation::Initiate 0 1 1
60
{block}
CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

VertexFromTRBScenterTriangulationMT 11 2 {12 if loosen=="1" else 11} 1
1.0
{yf}
0.3
0.2
{mf}
{nw}
1.0
{ps}
0.7
{mt}
0
{idx}
60
""")
PYEOF
}
run() { python3 $T/tools/port/compare.py _mt.model twoSquare3D_MT.init euler3.rk5 --tol=1e-9 ${1} 2>&1 | tail -1; }

echo "Material (MF) flags, forces and stored state:"
for mf in 0 2 4 6 7 8 9; do
  printf "  %-42s " "MF flag $mf"; gen $mf 0.0 1 1 0; run
done
# MF 3 reads Y_fibre as an anisotropy and divides by (2 - it), so 2.0 is a
# pole; MF 10 needs the 12th index, the loosening compound.
printf "  %-42s " "MF flag 3 (Y_fibre 0.5)"; gen 3 0.0 1 1 0 0.5; run
printf "  %-42s " "MF flag 10 (loosening compound)"; gen 10 0.0 1 1 0 2.0 1; run
printf "  %-42s " "MF flag 0, plane strain"; gen 0 0.0 1 0 0; run
printf "  %-42s " "MF flag 0, no wall dynamics"; gen 0 0.0 0 1 0; run

echo
echo "Forces only (MF flags whose derivs writes a cell variable, so the"
echo "t=0 print differs - see NOTES.md; the trajectory is identical):"
printf "  %-42s " "MF flag 1 (FiberModel)"; gen 1 0.0 1 1 0; run --block=vertex
printf "  %-42s " "MF flag 5"; gen 5 0.0 1 1 0; run --block=vertex
printf "  %-42s " "MF flag 0, MT update flag 1 (TETA)"; gen 0 0.0 1 1 1; run --block=vertex

echo
echo "Neighbour weighting - expected to MISMATCH: legacy's guard is"
echo "'size_t > -1', always false, so its averaging never happens (README 13)."
for w in 0.2 0.4; do
  printf "  %-42s " "neighbour weight $w"; gen 0 $w 1 1 0; run
done
echo
echo "  Stored stress anisotropy as the weight rises. The fixture's two cells"
echo "  agree to <1%, so a real average must barely move; legacy tracks (1-w)."
printf "    %-10s %-16s %-16s\n" "weight" "legacy" "v2 (fixed)"
printf 'Euler\n0 0.01\n0 2\n0.01\n' > _one.rk5
for w in 0.0001 0.2 0.4 0.8; do
  gen 0 $w 0 1 0
  L=$($T/bin/simulator _mt.model twoSquare3D_MT.init _one.rk5 2>/dev/null | awk '$1==4 && NF>60 {v=$11} END{print v}')
  N=$($T/build/simulator _mt.model twoSquare3D_MT.init _one.rk5 2>/dev/null | awk '$1==4 && NF>60 {v=$11} END{print v}')
  printf "    %-10s %-16s %-16s\n" "$w" "$L" "$N"
done

echo
echo "VertexFromTRBScenterTriangulationConcentrationHillMT:"
cat > _chmt.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

VertexFromTRBScenterTriangulationConcentrationHillMT 8 2 3 1
1.0
3.0
0.3
0.5
2.0
0.2
0.5
2.0
0 4 1
60
M
printf "  %-42s " "two levels (no stored directions)"
python3 $T/tools/port/compare.py _chmt.model twoSquare3D_MT.init euler3.rk5 --tol=1e-9 2>&1 | tail -1
# With storage, derivs writes cell variables, so only the t=0 print differs -
# the same derivs-side-effect timing as elsewhere; every later print agrees.
cat > _chmt6.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

VertexFromTRBScenterTriangulationConcentrationHillMT 8 6 3 1 1 1 1 1
1.0
3.0
0.3
0.5
2.0
0.2
0.5
2.0
0 4 1
60
15
31
35
41
M
printf "  %-42s " "six levels, forces"
python3 $T/tools/port/compare.py _chmt6.model twoSquare3D_MT.init euler3.rk5 --tol=1e-9 --block=vertex 2>&1 | tail -1
