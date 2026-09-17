#!/usr/bin/env python3
"""Independent check of the two Bending reactions whose legacy sources have a
wrong-vertex bug, so tools/port/compare.py cannot be used.

Re-derives the force law from the legacy algebra in Python, applies one
explicit Euler step of size h, and checks the simulator lands in the same
place. Also asserts the properties the fixes are supposed to restore:
Bending::Angle's three per-turn contributions must sum to zero (it is a
gradient of an energy), which the legacy version violates.
"""
import math, re, subprocess, sys

SIM = "/Users/ross/projects/tissue/build/simulator"
PI = 3.14159


def read_init(path):
    toks = []
    for line in open(path):
        line = line.split("#", 1)[0]
        toks += line.split()
    it = iter(toks)
    nxt = lambda: next(it)
    nC, nW, nV = int(nxt()), int(nxt()), int(nxt())
    walls = []
    for _ in range(nW):
        nxt()                                   # edge label
        nxt(); nxt()                            # cell1 cell2
        walls.append((int(nxt()), int(nxt())))
    nV2, dim = int(nxt()), int(nxt())
    pos = [[float(nxt()) for _ in range(dim)] for _ in range(nV2)]
    nE, nL, nWV = int(nxt()), int(nxt()), int(nxt())
    wdata = [[float(nxt()) for _ in range(nL + nWV)] for _ in range(nE)]
    return walls, pos, wdata


def cycles(walls, cells):
    """Pair each cell's vertex cycle with the walls joining consecutive pairs."""
    out = []
    for verts in cells:
        n = len(verts)
        ws = []
        for k in range(n):
            a, b = verts[k], verts[(k + 1) % n]
            ws.append(next(i for i, (u, v) in enumerate(walls)
                           if {u, v} == {a, b}))
        out.append((verts, ws))
    return out


def neighbor_center(k_bend, Li, cyc, pos, wdata, dim):
    F = [[0.0] * dim for _ in pos]
    for verts, ws in cyc:
        n = len(verts)
        for k in range(n):
            j, jp, jm = verts[k], verts[(k + 1) % n], verts[k - 1]
            ep, em = ws[k], ws[k - 1]
            Lem, Lep = wdata[em][Li], wdata[ep][Li]
            for d in range(dim):
                target = (pos[jp][d] * Lem + pos[jm][d] * Lep) / (Lem + Lep)
                F[j][d] -= k_bend * (pos[j][d] - target)
    return F


def angle(k_bend, Ti, cyc, pos, wdata, dim):
    F = [[0.0] * dim for _ in pos]
    worst_sum = 0.0
    for verts, ws in cyc:
        n = len(verts)
        for k in range(n):
            j, jp, jm = verts[k], verts[(k + 1) % n], verts[k - 1]
            em = ws[k - 1]
            dm = [pos[j][d] - pos[jm][d] for d in range(dim)]
            dp = [pos[jp][d] - pos[j][d] for d in range(dim)]
            f = sum(a * b for a, b in zip(dm, dp))
            Lm = math.sqrt(sum(a * a for a in dm))
            Lp = math.sqrt(sum(a * a for a in dp))
            g = 1.0 / (Lp * Lm)
            Fc = min(0.999, max(-0.999, f * g))
            theta = math.acos(Fc) - PI
            f0 = k_bend * (theta - wdata[em][Ti]) * g / math.sqrt(1 - Fc * Fc)
            f1 = f * g * Lm / Lp
            f2 = f * g * Lp / Lm
            for d in range(dim):
                a = f0 * (dm[d] * f2 - dp[d])
                b = f0 * (dp[d] * (1 + f1) - dm[d] * (1 + f2))
                c = f0 * (dm[d] - dp[d] * f1)
                worst_sum = max(worst_sum, abs(a + b + c))
                F[jm][d] += a
                F[j][d] += b
                F[jp][d] += c
    return F, worst_sum


def sim_positions(model, init, h, nvert, dim):
    open("tiny.rk5", "w").write(f"Euler\n0 {h}\n0 2\n{h}\n")
    out = subprocess.run([SIM, model, init, "tiny.rk5"], capture_output=True,
                         text=True).stdout
    # Last printed vertex block: "<nvert> <dim>" then nvert rows.
    m = [x.start() for x in re.finditer(rf"^{nvert} {dim}$", out, re.M)]
    toks = out[m[-1]:].split()[2:2 + nvert * dim]
    return [[float(t) for t in toks[i * dim:(i + 1) * dim]] for i in range(nvert)]


def check(name, model, F, pos, h, dim):
    got = sim_positions(model, "twoSquare.init", h, len(pos), dim)
    worst = 0.0
    for v in range(len(pos)):
        for d in range(dim):
            want = pos[v][d] + h * F[v][d]
            worst = max(worst, abs(want - got[v][d]) / max(1.0, abs(want)))
    print(f"{name:<44} {'agrees' if worst < 2e-6 else 'DISAGREES'} "
          f"(worst {worst:.2e} vs independent Python force law)")
    return worst < 2e-6


walls, pos, wdata = read_init("twoSquare.init")
CELLS = [[0, 1, 3, 2], [2, 3, 5, 4]]   # sorted cycles, as printed by the simulator
cyc = cycles(walls, CELLS)
dim = len(pos[0])
ok = True

open("nc.model", "w").write("1 0 0\n\nBending::NeighborCenter 1 1 1\n0.5\n0\n")
ok &= check("Bending::NeighborCenter (fixed)", "nc.model",
            neighbor_center(0.5, 0, cyc, pos, wdata, dim), pos, 1.0, dim)

open("an.model", "w").write("1 0 0\n\nBending::Angle 1 1 1\n0.5\n1\n")
F, worst_sum = angle(0.5, 1, cyc, pos, wdata, dim)
ok &= check("Bending::Angle (fixed)", "an.model", F, pos, 1.0, dim)
print(f"{'Bending::Angle force triples sum to zero':<44} "
      f"{'yes' if worst_sum < 1e-12 else 'NO'} (worst |sum| {worst_sum:.1e})")
sys.exit(0 if ok else 1)
