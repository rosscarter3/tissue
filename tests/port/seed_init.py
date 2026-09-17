#!/usr/bin/env python3
"""Re-emit a tissue init file with seeded cell and wall variables.

The transport reactions read paired wall variables (PIN, AUX/LAX, wall auxin),
and every init shipped with the tutorials leaves wall variables at zero - so a
comparison against legacy passes trivially with no transport happening at all.
This writes a copy with deterministic non-zero values in both tables.

    python3 seed_init.py IN.init OUT.init [numWallVar] [numCellVar]
"""
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    nwv = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    ncv = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    toks = []
    for line in open(src):
        toks += line.split("#", 1)[0].split()
    it = iter(toks)
    nxt = lambda: next(it)

    nC, nW, nV = int(nxt()), int(nxt()), int(nxt())
    conn = [[nxt() for _ in range(5)] for _ in range(nW)]
    nV2, dim = int(nxt()), int(nxt())
    pos = [[nxt() for _ in range(dim)] for _ in range(nV2)]
    nE, nL, oldWV = int(nxt()), int(nxt()), int(nxt())
    lengths = []
    for _ in range(nE):
        lengths.append([nxt() for _ in range(nL)])
        for _ in range(oldWV):
            nxt()                         # discard the old wall variables

    out = [f"{nC} {nW} {nV}"]
    out += [" ".join(c) for c in conn]
    out += ["", f"{nV2} {dim}"]
    out += [" ".join(p) for p in pos]
    out += ["", f"{nE} {nL} {nwv}"]
    for w in range(nE):
        # a cheap deterministic spread in (0.05, 0.95), different per column
        vals = [f"{0.05 + 0.9 * (((w + 1) * (k + 3) * 2654435761) % 997) / 997:.6f}"
                for k in range(nwv)]
        out.append(" ".join(lengths[w] + vals))
    out += ["", f"{nC} {ncv}"]
    for c in range(nC):
        vals = [f"{0.05 + 0.9 * (((c + 2) * (k + 5) * 40503) % 991) / 991:.6f}"
                for k in range(ncv)]
        out.append(" ".join(vals))
    open(dst, "w").write("\n".join(out) + "\n")
    print(f"{dst}: {nC} cells x {ncv} vars, {nE} walls x {nwv} vars")


main()
