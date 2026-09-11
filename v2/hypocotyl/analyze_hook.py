#!/usr/bin/env python3
"""Analyze hook-opening simulation output (simulator print flag 0).

Computes per timepoint: hook angle (angle between basal and apical axis
directions; 160 deg = closed, ~0 = fully open), inner/outer epidermal arc
lengths in the hook region, and compares with the experimental opening
kinetics from Walia, Carter et al. (2024), Figure 1/Data S1.
"""
import math
import sys

# generator layout (must match make_hook_init.py)
NFILES = 10
NV_PER_LINE = NFILES + 1
N_BASAL, N_HOOK, N_APICAL = 8, 32, 3
NSLICES = N_BASAL + N_HOOK + N_APICAL
T_END = 10.0

# experimental hook opening (Walia et al. 2024, mean angle vs h after light)
EXP_T = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
EXP_A = [158.7, 155.6, 133.4, 102.9, 81.5, 63.2, 53.3, 46.8, 44.3, 40.4, 34.5]


def vid(l, fb):
    return l * NV_PER_LINE + fb


def parse_flag0(path):
    tok = open(path).read().split()
    pos = 0

    def take(n):
        nonlocal pos
        out = tok[pos:pos + n]
        pos += n
        return out

    num_print = int(take(1)[0])  # header (tCount==0 only)
    frames = []
    while pos < len(tok):
        nv, dim = int(tok[pos]), int(tok[pos + 1])
        pos += 2
        coords = [float(x) for x in take(nv * dim)]
        verts = [(coords[2 * i], coords[2 * i + 1]) for i in range(nv)]
        nc, ncol = int(tok[pos]), int(tok[pos + 1])
        pos += 2
        for _ in range(nc):
            nvert = int(take(1)[0])
            take(nvert + ncol)  # vertex ids + variables+volume+numWall
        nw, wcol = int(tok[pos]), int(tok[pos + 1])
        pos += 2
        take(nw * (wcol + 2))  # v1 v2 + wall columns
        frames.append(verts)
    return num_print, frames


def angle_between(u, v):
    dot = u[0] * v[0] + u[1] * v[1]
    nu = math.hypot(*u)
    nv = math.hypot(*v)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (nu * nv)))))


def analyze(path):
    _, frames = parse_flag0(path)
    n = len(frames)
    times = [T_END * k / (n - 1) for k in range(n)]
    rows = []
    for t, verts in zip(times, frames):
        def axis(l):
            return verts[vid(l, 5)]  # file boundary 5 = central axis
        def vec(a, b):
            return (b[0] - a[0], b[1] - a[1])
        basal = vec(axis(2), axis(6))
        apical = vec(axis(NSLICES - 2), axis(NSLICES))
        hook_angle = angle_between(basal, apical)
        def arc(fb):
            s = 0.0
            for l in range(N_BASAL, N_BASAL + N_HOOK):
                a, b = verts[vid(l, fb)], verts[vid(l + 1, fb)]
                s += math.hypot(b[0] - a[0], b[1] - a[1])
            return s
        rows.append((t, hook_angle, arc(NV_PER_LINE - 1), arc(0)))
    return rows


if __name__ == "__main__":
    rows = analyze(sys.argv[1] if len(sys.argv) > 1 else "run0.out")
    print(f"{'t(h)':>6} {'angle':>7} {'inner':>8} {'outer':>8} {'exp.angle':>9}")
    for (t, a, li, lo) in rows:
        exp = ""
        for et, ea in zip(EXP_T, EXP_A):
            if abs(et - t) < 0.13:
                exp = f"{ea:9.1f}"
        print(f"{t:6.2f} {a:7.1f} {li:8.1f} {lo:8.1f} {exp}")
    t0, a0, li0, lo0 = rows[0]
    t1, a1, li1, lo1 = rows[-1]
    print(f"\nangle {a0:.1f} -> {a1:.1f} deg (exp: 158.7 -> 34.5)")
    print(f"inner epidermal hook arc: {li0:.1f} -> {li1:.1f} um "
          f"({li1 / li0:.2f}x; exp inner fold ~2-4x)")
    print(f"outer epidermal hook arc: {lo0:.1f} -> {lo1:.1f} um "
          f"({lo1 / lo0:.2f}x; exp outer fold ~1.07x)")
