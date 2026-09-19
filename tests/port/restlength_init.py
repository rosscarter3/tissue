#!/usr/bin/env python3
"""Re-emit a tissue init file with the walls' resting lengths perturbed.

Nearly every shipped init has each wall's resting length equal to its actual
length, so a model whose only force is a spring starts at equilibrium and
never moves. That makes a whole class of reactions untestable: the clamps,
the boundary conditions and the velocity diagnostics all act on vertex
derivatives, so with nothing moving they "match" legacy while doing nothing.

    python3 restlength_init.py IN.init OUT.init [scale] [spread]

`scale` multiplies every resting length (default 0.8, so the springs pull
inwards). `spread` additionally varies them wall to wall by up to that
fraction (default 0.1), which matters whenever a reaction averages over a
group of vertices: with uniform lengths a symmetric template gives every
vertex in the group the same velocity, the mean equals each of them, and the
averaging is still a no-op.
"""
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    scale = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
    spread = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1

    toks = []
    for line in open(src, errors="ignore"):
        toks += line.split("#", 1)[0].split()
    it = iter(toks)
    nxt = lambda: next(it)

    nC, nW, nV = int(nxt()), int(nxt()), int(nxt())
    conn = [[nxt() for _ in range(5)] for _ in range(nW)]
    nV2, dim = int(nxt()), int(nxt())
    pos = [[nxt() for _ in range(dim)] for _ in range(nV2)]
    nE, nL, nwv = int(nxt()), int(nxt()), int(nxt())
    walls = []
    for w in range(nE):
        lengths = [float(nxt()) for _ in range(nL)]
        # A deterministic wobble in [-spread, +spread], different per wall.
        wobble = 1.0 + spread * (2.0 * (((w + 1) * 2654435761) % 1000) / 999.0
                                 - 1.0)
        lengths = [f"{x * scale * wobble:.10g}" for x in lengths]
        walls.append(lengths + [nxt() for _ in range(nwv)])
    nC2, ncv = int(nxt()), int(nxt())
    cells = [[nxt() for _ in range(ncv)] for _ in range(nC2)]

    out = [f"{nC} {nW} {nV}"]
    out += [" ".join(c) for c in conn]
    out += ["", f"{nV2} {dim}"]
    out += [" ".join(p) for p in pos]
    out += ["", f"{nE} {nL} {nwv}"]
    out += [" ".join(w) for w in walls]
    out += ["", f"{nC2} {ncv}"]
    out += [" ".join(c) for c in cells]
    open(dst, "w").write("\n".join(out) + "\n")
    print(f"{dst}: {nE} walls, resting lengths x{scale} +/-{spread:.0%}")


if __name__ == "__main__":
    main()
