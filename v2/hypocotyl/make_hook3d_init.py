#!/usr/bin/env python3
"""Generate a 3D apical hook: the epidermal shell of a bent hypocotyl tube
as a closed center-triangulated surface mesh (the toroidal pressure vessel
of Walia, Carter et al. 2024, cell-resolved).

Tube of radius r = 75.8 um bent through BEND_DEG about R = 2.88 r, with
straight basal/apical segments, closed by polygonal end caps. Cells are
surface quads (rings x sectors); sector theta=0 faces the hook center
(inner side).

Wall variables: [restLength, growFlag, bendFlag]
  growFlag = 1 on axial walls (longitudinal growth only)
  bendFlag = 1 on axial walls (chains for bending stiffness)
Cell variables: [0,0,0,0, auxin, cellType, thetaDeg] (CT data appended at 7)
  auxin: 0.6 baseline (dark maintenance), up to 1.0 on the inner hook side
  cellType: 0 = shell quad, 1 = end cap
"""
import math
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "hook3d.init"

R_HYP = 75.8
LAMBDA_R = 2.88
R_AXIS = LAMBDA_R * R_HYP
BEND_DEG = 175.0            # relaxes toward ~160 deg under turgor
BASAL_LEN = 50.0   # short arms: minimize rotational drag
APICAL_LEN = 40.0
N_BASAL, N_HOOK, N_APICAL = 3, 24, 2
NCIRC = 12

AUXIN_BASE = 0.6
AUXIN_MAX = 1.0

# --- axis path (bend in the xy plane) -----------------------------------------
bend = math.radians(BEND_DEG)
axis = []  # (point, normal-toward-center)
for k in range(N_BASAL + 1):
    y = -BASAL_LEN + BASAL_LEN * k / N_BASAL
    axis.append(((0.0, y), (1.0, 0.0)))
for k in range(1, N_HOOK + 1):
    phi = bend * k / N_HOOK
    p = (R_AXIS - R_AXIS * math.cos(phi), R_AXIS * math.sin(phi))
    n = (math.cos(phi), -math.sin(phi))
    axis.append((p, n))
end_p, _ = axis[-1]
tangent = (math.sin(bend), math.cos(bend))
for k in range(1, N_APICAL + 1):
    s = APICAL_LEN * k / N_APICAL
    axis.append(((end_p[0] + tangent[0] * s, end_p[1] + tangent[1] * s),
                 axis[-1][1]))

NRINGS = len(axis)
NAX = NRINGS - 1           # axial cell slices
NQUAD = NAX * NCIRC
NCELLS = NQUAD + 2         # + basal cap, apical cap
CAP_B, CAP_A = NQUAD, NQUAD + 1

# --- vertices: ring i, sector j  (theta=0 faces the curvature center) ----------
verts = []
for (p, n) in axis:
    for j in range(NCIRC):
        th = 2.0 * math.pi * j / NCIRC
        c, s = math.cos(th), math.sin(th)
        verts.append((p[0] + R_HYP * c * n[0],
                      p[1] + R_HYP * c * n[1],
                      R_HYP * s))

def vid(i, j):
    return i * NCIRC + (j % NCIRC)

def qid(i, j):
    return i * NCIRC + (j % NCIRC)

def dist(a, b):
    return math.dist(verts[a], verts[b])

# --- walls ---------------------------------------------------------------------
# (cell1, cell2, v1, v2, growFlag, bendFlag)
walls = []
# circumferential edges: ring i, sector j
for i in range(NRINGS):
    for j in range(NCIRC):
        below = qid(i - 1, j) if i > 0 else CAP_B
        above = qid(i, j) if i < NAX else CAP_A
        walls.append((below, above, vid(i, j), vid(i, j + 1), 0, 0))
# axial edges: between rings i,i+1 at sector j
for i in range(NAX):
    for j in range(NCIRC):
        walls.append((qid(i, j - 1), qid(i, j), vid(i, j), vid(i + 1, j), 1, 1))

# --- cells ----------------------------------------------------------------------
cells = []
for i in range(NAX):
    if N_BASAL <= i < N_BASAL + N_HOOK:
        ax_w = math.sin(math.pi * (i - N_BASAL + 0.5) / N_HOOK) ** 0.25
        phi_deg = (i - N_BASAL + 0.5) / N_HOOK * BEND_DEG
    else:
        ax_w, phi_deg = 0.0, 0.0
    for j in range(NCIRC):
        th = 2.0 * math.pi * (j + 0.5) / NCIRC
        circ_w = max(0.0, math.cos(th))  # inner half of the tube
        auxin = AUXIN_BASE + (AUXIN_MAX - AUXIN_BASE) * ax_w * circ_w
        cells.append((auxin, 0, phi_deg))
cells.append((AUXIN_BASE, 1, 0.0))  # basal cap
cells.append((AUXIN_BASE, 1, 0.0))  # apical cap

# --- write ------------------------------------------------------------------------
with open(OUT, "w") as f:
    f.write(f"{NCELLS} {len(walls)} {len(verts)}\n")
    for w, (c1, c2, v1, v2, _, _) in enumerate(walls):
        f.write(f"{w} {c1} {c2} {v1} {v2}\n")
    f.write(f"\n{len(verts)} 3\n")
    for (x, y, z) in verts:
        f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
    f.write(f"\n{len(walls)} 1 2\n")
    for (c1, c2, v1, v2, grow, bendf) in walls:
        f.write(f"{dist(v1, v2):.6f} {grow} {bendf}\n")
    f.write(f"\n{NCELLS} 7\n")
    for (auxin, ctype, phi) in cells:
        f.write(f"0 0 0 0 {auxin:.4f} {ctype} {phi:.2f}\n")

inner_len = dist(vid(N_BASAL, 0), vid(N_BASAL + 1, 0))
outer_len = dist(vid(N_BASAL, NCIRC // 2), vid(N_BASAL + 1, NCIRC // 2))
print(f"wrote {OUT}: {NCELLS} cells ({NQUAD} quads + 2 caps), "
      f"{len(walls)} walls, {len(verts)} vertices")
print(f"hook axial cell lengths: inner ~{inner_len:.1f} um, "
      f"outer ~{outer_len:.1f} um; circumferential width "
      f"{2*math.pi*R_HYP/NCIRC:.1f} um")
