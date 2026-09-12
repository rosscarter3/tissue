#!/usr/bin/env python3
"""Validate the dynamic anisotropic material against the paper's predictions
(Walia, Carter et al. 2024, Figs 3E-F, 4E-G):

  1. Closed hook (dark / t=0): stress anisotropy > 0 with the principal
     stress circumferential everywhere (the model's geometric prediction that
     matches dark CMT orientations).
  2. During light opening the principal stress switches circumferential ->
     longitudinal, on the OUTER hook side before the INNER side.

Reads print-flag-0 output; signed anisotropy s = +a if the principal stress
is circumferential (angle to local tube axis > 45 deg), -a if longitudinal
(the paper's plotting convention)."""
import math
import sys

NCIRC = 12
N_BASAL, N_HOOK, N_APICAL = 3, 24, 2
NRINGS = N_BASAL + N_HOOK + N_APICAL + 1
ANISO_VAR = 7  # cell variables: [.,.,.,., auxin, type, theta, a, nx, ny, nz, s1]


def vid(i, j):
    return i * NCIRC + (j % NCIRC)


def parse(path):
    tok = open(path).read().split()
    pos = 1
    frames = []
    try:
        while pos < len(tok):
            nv, dim = int(tok[pos]), int(tok[pos + 1])
            pos += 2
            coords = [float(x) for x in tok[pos:pos + nv * dim]]
            if len(coords) < nv * dim:
                break
            pos += nv * dim
            verts = [tuple(coords[dim * i:dim * i + 3]) for i in range(nv)]
            nc, ncol = int(tok[pos]), int(tok[pos + 1])
            pos += 2
            cells = []
            for _ in range(nc):
                nvert = int(tok[pos])
                row = tok[pos + 1 + nvert:pos + 1 + nvert + ncol]
                cells.append([float(x) for x in row])
                pos += 1 + nvert + ncol
            nw, wcol = int(tok[pos]), int(tok[pos + 1])
            pos += 2 + nw * (wcol + 2)
            frames.append((verts, cells))
    except (ValueError, IndexError):
        pass
    return frames


def ring_center(verts, i):
    xs = [verts[vid(i, j)] for j in range(NCIRC)]
    return tuple(sum(p[d] for p in xs) / NCIRC for d in range(3))


def signed_aniso(verts, cells, i, j):
    """+a if principal stress circumferential, -a if longitudinal."""
    c = cells[i * NCIRC + j]
    a = c[ANISO_VAR]
    n = c[ANISO_VAR + 1:ANISO_VAR + 4]
    nn = math.sqrt(sum(x * x for x in n))
    if nn < 1e-12:
        return 0.0
    p, q = ring_center(verts, i), ring_center(verts, i + 1)
    ax = tuple(b - a2 for a2, b in zip(p, q))
    an = math.sqrt(sum(x * x for x in ax))
    cosang = abs(sum(x * y for x, y in zip(ax, n))) / (an * nn)
    longitudinal = cosang > math.cos(math.radians(45))
    return -a if longitudinal else a


def side_mean(verts, cells, sectors):
    vals = []
    mid = N_BASAL + N_HOOK // 2
    for i in range(mid - 5, mid + 5):  # central hook slices
        for j in sectors:
            vals.append(signed_aniso(verts, cells, i, j))
    return sum(vals) / len(vals)


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "run3d_light.out"
    t_end = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
    frames = parse(path)
    n = len(frames)
    inner_sectors = [0, NCIRC - 1]
    outer_sectors = [NCIRC // 2 - 1, NCIRC // 2]
    print(f"{'t(h)':>6} {'inner s':>8} {'outer s':>8}   "
          "(+ = circumferential principal stress, - = longitudinal)")
    switch = {}
    prev = {}
    for fi, (verts, cells) in enumerate(frames):
        t = t_end * fi / (n - 1) if n > 1 else 0.0
        si = side_mean(verts, cells, inner_sectors)
        so = side_mean(verts, cells, outer_sectors)
        for name, s in (("inner", si), ("outer", so)):
            if name in prev and prev[name] > 0 >= s and name not in switch:
                switch[name] = t
            prev[name] = s
        if fi % 4 == 0 or fi == n - 1:
            print(f"{t:6.2f} {si:+8.3f} {so:+8.3f}")
    print()
    for name in ("outer", "inner"):
        print(f"{name} switch (circ -> long): "
              f"{switch[name]:.2f} h" if name in switch else
              f"{name}: no switch within the run")
