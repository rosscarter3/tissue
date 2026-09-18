#!/usr/bin/env python3
"""Re-emit a tissue init file with explicit cell-variable values.

seed_init.py spreads its values over (0.05, 0.95), which is the right range
for concentrations but makes some reactions vacuous: FiberModel::Deposition
compares a cell's maximal stress against a hard-coded window of [8, 14], so
every seeded cell falls below the floor, every target share is zero and a
comparison against legacy passes without the redistribution ever running.
This writes a copy whose cell variables are whatever the caller asks for, so
a fixture can put cells above, inside and below such a window on purpose.

    python3 cellvar_init.py IN.init OUT.init "0.3 10 1 .5 .25 .1" "0.75 16 2 ..."

One quoted row of values per cell; every row must be the same length. Wall
variables are carried through unchanged.
"""
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    rows = [r.split() for r in sys.argv[3:]]
    if not rows or len({len(r) for r in rows}) != 1:
        sys.exit("give one quoted row of values per cell, all the same length")

    toks = []
    for line in open(src):
        toks += line.split("#", 1)[0].split()
    it = iter(toks)
    nxt = lambda: next(it)

    nC, nW, nV = int(nxt()), int(nxt()), int(nxt())
    if len(rows) != nC:
        sys.exit(f"{src} has {nC} cells but {len(rows)} rows were given")
    conn = [[nxt() for _ in range(5)] for _ in range(nW)]
    nV2, dim = int(nxt()), int(nxt())
    pos = [[nxt() for _ in range(dim)] for _ in range(nV2)]
    nE, nL, nwv = int(nxt()), int(nxt()), int(nxt())
    walls = [[nxt() for _ in range(nL + nwv)] for _ in range(nE)]

    out = [f"{nC} {nW} {nV}"]
    out += [" ".join(c) for c in conn]
    out += ["", f"{nV2} {dim}"]
    out += [" ".join(p) for p in pos]
    out += ["", f"{nE} {nL} {nwv}"]
    out += [" ".join(w) for w in walls]
    out += ["", f"{nC} {len(rows[0])}"]
    out += [" ".join(r) for r in rows]
    open(dst, "w").write("\n".join(out) + "\n")
    print(f"{dst}: {nC} cells x {len(rows[0])} vars, {nE} walls x {nwv} vars")


if __name__ == "__main__":
    main()
