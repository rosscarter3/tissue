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

echo
echo "VertexFromTRBSMT (same material, no center triangulation):"
cat > _tmt.model <<'M'
2 0 0

WallGrowth::Constant 2 1 1
0.05
1
0

VertexFromTRBSMT 10 1 10
1.0
2.0
0.3
0.2
0
0.0
0.0
1
0.7
0
0 1 4 5 6 7 8 9 10 20
M
printf "  %-42s " "MF flag 0, forces"
python3 $T/tools/port/compare.py _tmt.model tri3D_MT.init euler3.rk5 --tol=1e-9 --block=vertex 2>&1 | tail -1
printf "  %-42s " "MF flag 1, forces"
python3 - <<'PYEOF'
s = open("_tmt.model").read().replace("0.2\n0\n0.0", "0.2\n1\n0.0", 1)
open("_tmt1.model", "w").write(s)
PYEOF
python3 $T/tools/port/compare.py _tmt1.model tri3D_MT.init euler3.rk5 --tol=1e-9 --block=vertex 2>&1 | tail -1
echo "  The stored diagnostics are written from derivs, so legacy samples them"
echo "  one derivative evaluation later. The gap is exactly linear in the step:"
printf "    %-10s %-14s %-14s %s\n" "h" "legacy" "v2" "difference"
for h in 0.001 0.0005 0.00025; do
  printf 'Euler\n0 2\n0 6\n%s\n' $h > _eh.rk5
  A=$($T/bin/simulator _tmt.model tri3D_MT.init _eh.rk5 2>/dev/null | awk '$1==3 && NF>45 {n++; if(n==3) print $11}')
  B=$($T/build/simulator _tmt.model tri3D_MT.init _eh.rk5 2>/dev/null | awk '$1==3 && NF>45 {n++; if(n==3) print $11}')
  python3 -c "print(f'    {$h:<10} {$A:<14} {$B:<14} {abs($A-$B):.3e}')"
done

echo
echo "VertexFromTRLScenterTriangulationMT (orthotropic linear elements):"
trls() { python3 - "$1" "$2" "$3" "$4" <<'PYEOF'
import sys
ym, yf, mf, nw = sys.argv[1:5]
open("_trls.model", "w").write(f"""3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

VertexFromTRLScenterTriangulationMT 11 2 11 1
{ym}
{yf}
0.3
0.2
{mf}
{nw}
1.0
1
0.7
0
0
0 1 4 5 6 7 8 9 10 20 26
60
""")
PYEOF
}
runv() { python3 $T/tools/port/compare.py _trls.model twoSquare3D_MT.init euler3.rk5 --tol=1e-9 --block=vertex 2>&1 | tail -1; }
for mf in 0 1 5; do printf "  %-42s " "MF flag $mf, forces"; trls 1.0 2.0 $mf 0.0; runv; done
# MF 2 reads Y_matrix as a total, so it needs Y_matrix > Y_fibre.
printf "  %-42s " "MF flag 2 (Y_total 5, Y_fibre 2)"; trls 5.0 2.0 2 0.0; runv
printf "  %-42s " "neighbour weight 0.4"; trls 1.0 2.0 0 0.4; runv
echo "  Its stored diagnostics are sampled one evaluation later by legacy,"
echo "  the same as VertexFromTRBSMT; the gap is linear in the step:"
printf "    %-10s %-14s %-14s %s\n" "h" "legacy" "v2" "difference"
trls 1.0 2.0 0 0.0
for h in 0.001 0.0005 0.00025; do
  printf 'Euler\n0 2\n0 6\n%s\n' $h > _eh.rk5
  A=$($T/bin/simulator _trls.model twoSquare3D_MT.init _eh.rk5 2>/dev/null | awk '$1==4 && NF>60 {n++; if(n==3) print $12}')
  B=$($T/build/simulator _trls.model twoSquare3D_MT.init _eh.rk5 2>/dev/null | awk '$1==4 && NF>60 {n++; if(n==3) print $12}')
  python3 -c "print(f'    {$h:<10} {$A:<14} {$B:<14} {abs($A-$B):.3e}')"
done

echo
echo "Hypocotyl3D::VertexFromTRBScenterTriangulationMT:"
# Its MF flag -1 is a tissue-layer material keyed on cell variable 37, so the
# fixture needs real layer codes rather than seeded noise.
python3 - <<'PYEOF'
lines = open("twoSquare3D_MT.init").read().split("\n")
out, incell, seen = [], False, 0
for l in lines:
    t = l.split()
    if len(t) == 2 and t[0] == "2" and t[1] == "60":
        incell = True; out.append(l); continue
    if incell and len(t) == 60:
        t[37] = "-1" if seen == 0 else "-2"   # epidermis, then inner
        seen += 1; out.append(" ".join(t)); continue
    out.append(l)
open("twoSquare3D_HYP.init", "w").write("\n".join(out))
PYEOF
hyp() { python3 - "$1" "$2" "$3" "$4" <<'PYEOF'
import sys
mf, nw, p2, p3 = sys.argv[1:5]
open("_hyp.model", "w").write(f"""3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

Hypocotyl3D::VertexFromTRBScenterTriangulationMT 11 2 11 1
2.0
6.0
{p2}
{p3}
{mf}
{nw}
1.0
1
0.7
0
0
0 1 4 5 6 7 8 9 10 20 26
60
""")
PYEOF
}
runh() { python3 $T/tools/port/compare.py _hyp.model ${1:-twoSquare3D_MT.init} euler3.rk5 --tol=1e-9 2>&1 | tail -1; }
for mf in 0 1; do printf "  %-42s " "MF flag $mf"; hyp $mf 0.0 0.3 0.2; runh; done
# MF 2 reads Y_matrix as a total stiffness, so it needs Y_matrix > Y_fibre;
# below that the transverse modulus goes negative and legacy spins forever in
# its unbounded Jacobi loop (this port caps every one of them at 50 sweeps).
printf "  %-42s " "MF flag 2 (needs Y_total > Y_fibre)"
python3 - <<'PYEOF'
open("_hyp.model", "w").write("""3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.05
1
60

Hypocotyl3D::VertexFromTRBScenterTriangulationMT 11 2 11 1
9.0
6.0
0.3
0.2
2
0.0
1.0
1
0.7
0
0
0 1 4 5 6 7 8 9 10 20 26
60
""")
PYEOF
runh
printf "  %-42s " "MF flag -1 (layer material)"; hyp -1 0.0 0.5 10.0; runh twoSquare3D_HYP.init
echo "  (parameters 2 and 3 are layer scalings under MF -1, not Poisson"
echo "   ratios, so legacy skips its range check there and so does this.)"
printf "  %-42s " "MF -1, neighbour 0.4 (MISMATCH expected)"; hyp -1 0.4 0.5 10.0; runh twoSquare3D_HYP.init

echo
echo "VertexFromTRBScenterTriangulationMTOpt (an energy probe - it applies"
echo "no force; the annealing search it is named for was never finished):"
cat > _opt.model <<'M'
3 0 0

CenterTriangulation::Initiate 0 1 1
60

CenterTriangulation::WallGrowth::Constant 2 1 1
0.2
1
60

VertexFromTRBScenterTriangulationMTOpt 11 2 9 1
1.0
2.0
0.3
0.2
0.5
0
0.1
0.1
1.0
0.99
0.99
0 1 4 5 6 7 8 9 10
60
M
printf 'Euler\n0 0.2\n0 3\n0.01\n' > _opt.rk5
printf "  %-42s " "aniso flag 0 (energies on stdout)"
python3 $T/tools/port/compare.py _opt.model twoSquare3D_MT.init _opt.rk5 --tol=1e-9 2>&1 | tail -1
python3 - <<'PYEOF'
open("_opt1.model","w").write(open("_opt.model").read().replace("0.5\n0\n0.1\n","0.5\n1\n0.1\n",1))
PYEOF
printf "  %-42s " "aniso flag 1"
python3 $T/tools/port/compare.py _opt1.model twoSquare3D_MT.init _opt.rk5 --tol=1e-9 2>&1 | tail -1
# It must leave the simulation untouched.
python3 - <<'PYEOF'
s = open("_opt.model").read()
head = s.split("VertexFromTRBScenterTriangulationMTOpt")[0].rstrip()
open("_opt_none.model","w").write(head.replace("3 0 0","2 0 0",1)+"\n")
PYEOF
$T/build/simulator _opt.model twoSquare3D_MT.init _opt.rk5 2>/dev/null \
  | grep -vE "^-?[0-9.e+-]+  [0-9.e+-]" > _with.txt
$T/build/simulator _opt_none.model twoSquare3D_MT.init _opt.rk5 2>/dev/null > _without.txt
printf "  %-42s " "leaves the tissue untouched"
diff -q _with.txt _without.txt >/dev/null && echo "yes" || echo "NO"
