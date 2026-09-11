#!/usr/bin/env python3
"""Analyze 3D hook simulation output (print flag 0): hook angle from the
tube-axis (ring centroids) and inner/outer surface arc lengths in the hook."""
import math
import sys

NCIRC = 12
N_BASAL, N_HOOK, N_APICAL = 3, 24, 2
NRINGS = N_BASAL + N_HOOK + N_APICAL + 1
T_END = 10.0


def vid(i, j):
    return i * NCIRC + (j % NCIRC)


def parse_flag0(path):
    tok = open(path).read().split()
    pos = 1
    frames = []
    while pos < len(tok):
        nv, dim = int(tok[pos]), int(tok[pos + 1])
        pos += 2
        coords = [float(x) for x in tok[pos:pos + nv * dim]]
        pos += nv * dim
        verts = [tuple(coords[dim * i:dim * i + 3]) for i in range(nv)]
        nc, ncol = int(tok[pos]), int(tok[pos + 1])
        pos += 2
        for _ in range(nc):
            nvert = int(tok[pos])
            pos += 1 + nvert + ncol
        nw, wcol = int(tok[pos]), int(tok[pos + 1])
        pos += 2 + nw * (wcol + 2)
        frames.append(verts)
    return frames


def ring_center(verts, i):
    xs = [verts[vid(i, j)] for j in range(NCIRC)]
    return tuple(sum(p[d] for p in xs) / NCIRC for d in range(3))


def angle_between(u, v):
    dot = sum(a * b for a, b in zip(u, v))
    nu = math.sqrt(sum(a * a for a in u))
    nv = math.sqrt(sum(b * b for b in v))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (nu * nv)))))


def analyze(path, t_end=T_END):
    frames = parse_flag0(path)
    n = len(frames)
    rows = []
    for fi, verts in enumerate(frames):
        t = t_end * fi / (n - 1) if n > 1 else 0.0
        b1, b2 = ring_center(verts, 2), ring_center(verts, 6)
        a1, a2 = ring_center(verts, NRINGS - 3), ring_center(verts, NRINGS - 1)
        basal = tuple(q - p for p, q in zip(b1, b2))
        apical = tuple(q - p for p, q in zip(a1, a2))
        hook_angle = angle_between(basal, apical)

        def arc(j):
            s = 0.0
            for i in range(N_BASAL, N_BASAL + N_HOOK):
                s += math.dist(verts[vid(i, j)], verts[vid(i + 1, j)])
            return s
        rows.append((t, hook_angle, arc(0), arc(NCIRC // 2)))
    return rows


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "run3d_light.out"
    t_end = float(sys.argv[2]) if len(sys.argv) > 2 else T_END
    rows = analyze(path, t_end)
    print(f"{'t(h)':>6} {'angle':>7} {'inner':>8} {'outer':>8}")
    for (t, a, li, lo) in rows:
        print(f"{t:6.2f} {a:7.1f} {li:8.1f} {lo:8.1f}")
    r0, r1 = rows[0], rows[-1]
    print(f"\nangle {r0[1]:.1f} -> {r1[1]:.1f} deg")
    print(f"inner hook arc {r0[2]:.1f} -> {r1[2]:.1f} ({r1[2]/r0[2]:.2f}x)")
    print(f"outer hook arc {r0[3]:.1f} -> {r1[3]:.1f} ({r1[3]/r0[3]:.2f}x)")
