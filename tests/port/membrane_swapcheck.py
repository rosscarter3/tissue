#!/usr/bin/env python3
"""Check membrane-cycling reactions are invariant to a wall's cell1/cell2 order.

Which of a wall's two cells an init file stores as cell1 is arbitrary: the
paired wall variables are stored in that order too, so relabelling both
together must leave the dynamics untouched. Four reactions in legacy's
membraneCycling.cc and membraneCyclingAll.cc break that symmetry through index
and parameter typos, which is why they cannot be compared against legacy
directly (see NOTES.md). This tests the property instead.

    python3 membrane_swapcheck.py MODEL INIT SOLVER

Reports, for each simulator, whether the cell variables come out the same with
the wall orientation flipped.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BINS = [("legacy", os.path.join(ROOT, "bin", "simulator")),
        ("v2", os.path.join(ROOT, "build", "simulator"))]


def swap_init(src, dst):
    """Flip cell1/cell2 on every wall and swap each paired wall variable."""
    toks = []
    for line in open(src):
        toks += line.split("#", 1)[0].split()
    it = iter(toks)
    nxt = lambda: next(it)
    nC, nW, nV = int(nxt()), int(nxt()), int(nxt())
    conn = []
    for _ in range(nW):
        label, c1, c2, v1, v2 = (nxt() for _ in range(5))
        conn.append([label, c2, c1, v1, v2])          # the flip
    nV2, dim = int(nxt()), int(nxt())
    pos = [[nxt() for _ in range(dim)] for _ in range(nV2)]
    nE, nL, nWV = int(nxt()), int(nxt()), int(nxt())
    walls = []
    for _ in range(nE):
        length = [nxt() for _ in range(nL)]
        v = [nxt() for _ in range(nWV)]
        for k in range(0, nWV - 1, 2):                # swap each (k, k+1) pair
            v[k], v[k + 1] = v[k + 1], v[k]
        walls.append(length + v)
    nC2, nCV = int(nxt()), int(nxt())
    cells = [[nxt() for _ in range(nCV)] for _ in range(nC2)]

    out = [f"{nC} {nW} {nV}"] + [" ".join(c) for c in conn]
    out += ["", f"{nV2} {dim}"] + [" ".join(p) for p in pos]
    out += ["", f"{nE} {nL} {nWV}"] + [" ".join(w) for w in walls]
    out += ["", f"{nC2} {nCV}"] + [" ".join(c) for c in cells]
    open(dst, "w").write("\n".join(out) + "\n")
    return nC, nCV


def cell_values(binary, model, init, solver, nC, nCV):
    """The cell variables from the last printed state."""
    p = subprocess.run([binary, model, init, solver], capture_output=True,
                       text=True)
    if p.returncode != 0:
        sys.exit(f"{binary} failed:\n{p.stderr[-600:]}")
    lines = [l.split() for l in p.stdout.splitlines() if l.split()]
    blocks = [n for n, l in enumerate(lines)
              if len(l) == 2 and l[0] == str(nC) and l[1].isdigit()]
    start = blocks[-1] + 1
    vals = []
    for row in lines[start:start + nC]:
        nv = int(row[0])
        vals += [float(x) for x in row[1 + nv:1 + nv + nCV]]
    return vals


def main():
    model, init, solver = sys.argv[1:4]
    swapped = "_swapped.init"
    nC, nCV = swap_init(init, swapped)
    for name, binary in BINS:
        a = cell_values(binary, model, init, solver, nC, nCV)
        b = cell_values(binary, model, swapped, solver, nC, nCV)
        worst = max(abs(x - y) / max(1.0, abs(x), abs(y))
                    for x, y in zip(a, b)) if a else 0.0
        verdict = "invariant" if worst < 1e-9 else "ORIENTATION-DEPENDENT"
        print(f"  {name:<7} {verdict:<22} (worst {worst:.3e})")


main()
