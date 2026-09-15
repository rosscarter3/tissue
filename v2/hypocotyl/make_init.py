#!/usr/bin/env python3
"""Generate a biologically-parameterised 3D apical hook: the hypocotyl
epidermis as a closed, tapered, centre-triangulated surface with the auxin
and apoplastic-pH fields measured in Walia, Carter et al. (2024) Dev Cell
59:3245.

Geometry (paper's fitted toroid + hypocotyl anatomy)
  r = 75.8 um hypocotyl radius, hook curvature R = 2.88 r, hook angle ~160 deg
  16 epidermal cell files around the circumference (Arabidopsis has ~16-24)
  gentle taper over the apical segment where the cotyledons attach
  inner-flank cells are short and outer-flank cells long purely by geometry:
  on the arc the inner surface is at radius R-r and the outer at R+r

Chemistry (measured, Figs 2 and S3)
  auxin  cell var 5: DR5 maximum on the INNER hook side, baseline elsewhere.
                     Depleted after illumination (Degradation::One).
  acid   cell var 6: apoplastic acidification, 0 = dark/alkaline, 1 = fully
                     acidified. In the dark the INNER side is more alkaline
                     than the outer, so it starts lower. Illumination drives
                     both to 1 with tau = 0.5 h (measured: acidification at
                     both flanks by 30 min).

Wall variables: [restLength, growFlag, bendFlag, Lmax]
  Lmax is the cell's mature length. All hook epidermal cells mature to the
  SAME length; outer-flank cells are born near it (so they can only extend
  ~1.07x, as measured) while inner-flank cells are born short (so they reach
  ~2.2x, measured 2.15x). This single assumption reproduces both measured
  fold changes without a free inner/outer parameter.

Cell variables: [sdir x, sdir y, sdir z, anisotropy, sigma1,
                 auxin, acid, cellType, thetaDeg]   (CT data appended at 9)
The first four are written by the TRBS stress solver and land on VTK's native
"cell vector" / "cell vector length" slots, so the predicted CMT axis is a
first-class vector field in ParaView.
"""
import math
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "hook.init"

# --- geometry ---------------------------------------------------------------
R_HYP = 75.8                 # hypocotyl radius (um), fitted in the paper
LAMBDA_R = 2.88              # hook curvature ratio R/r, fitted in the paper
R_AXIS = LAMBDA_R * R_HYP
BEND_DEG = 172.0             # generated bend; relaxes to ~160 under turgor
BASAL_LEN, APICAL_LEN = 55.0, 45.0
N_BASAL, N_HOOK, N_APICAL = 3, 26, 3
NCIRC = 16
APEX_TAPER = 0.86            # radius factor at the cotyledon end

# --- chemistry (measured patterns) ------------------------------------------
AUXIN_BASE = 0.45            # transiting auxin, represses extension everywhere
AUXIN_PEAK = 1.00            # DR5 maximum, inner hook side
ACID_DARK_OUTER = 0.25       # outer flank is already somewhat acidic in dark
ACID_DARK_INNER = 0.05       # inner flank is more alkaline in the dark
MATURE_FACTOR = 1.05         # outer-flank extension, measured 1.055 at t=8 h
                             # (this alone fixes the inner cap at 2.17x,
                             #  measured 2.149 - not a free parameter)

bend = math.radians(BEND_DEG)
axis = []                                   # (point, inward normal)
for k in range(N_BASAL + 1):
    axis.append(((0.0, -BASAL_LEN + BASAL_LEN * k / N_BASAL), (1.0, 0.0)))
for k in range(1, N_HOOK + 1):
    phi = bend * k / N_HOOK
    axis.append(((R_AXIS - R_AXIS * math.cos(phi), R_AXIS * math.sin(phi)),
                 (math.cos(phi), -math.sin(phi))))
end_p = axis[-1][0]
tangent = (math.sin(bend), math.cos(bend))
for k in range(1, N_APICAL + 1):
    s = APICAL_LEN * k / N_APICAL
    axis.append(((end_p[0] + tangent[0] * s, end_p[1] + tangent[1] * s),
                 axis[-1][1]))

NRINGS = len(axis)
NAX = NRINGS - 1
NQUAD = NAX * NCIRC
NCELLS = NQUAD + 2
CAP_B, CAP_A = NQUAD, NQUAD + 1
APEX_START = N_BASAL + N_HOOK               # taper begins here

def ring_radius(i):
    if i <= APEX_START:
        return R_HYP
    f = (i - APEX_START) / max(1, NRINGS - 1 - APEX_START)
    return R_HYP * (1.0 - (1.0 - APEX_TAPER) * f)

verts = []
for i, (p, n) in enumerate(axis):
    rad = ring_radius(i)
    for j in range(NCIRC):
        th = 2.0 * math.pi * j / NCIRC       # theta = 0 faces the hook centre
        c, s = math.cos(th), math.sin(th)
        verts.append((p[0] + rad * c * n[0], p[1] + rad * c * n[1], rad * s))

def vid(i, j): return i * NCIRC + (j % NCIRC)
def qid(i, j): return i * NCIRC + (j % NCIRC)
def dist(a, b): return math.dist(verts[a], verts[b])

# --- walls: (cell1, cell2, v1, v2, growFlag, bendFlag) -----------------------
walls = []
for i in range(NRINGS):                      # circumferential
    for j in range(NCIRC):
        below = qid(i - 1, j) if i > 0 else CAP_B
        above = qid(i, j) if i < NAX else CAP_A
        walls.append((below, above, vid(i, j), vid(i, j + 1), 0, 0))
n_circ = len(walls)
for i in range(NAX):                         # axial
    for j in range(NCIRC):
        walls.append((qid(i, j - 1), qid(i, j), vid(i, j), vid(i + 1, j), 1, 1))

# Mature length: the longest axial wall in the hook (outer flank) times the
# measured 1.07x. Every axial wall shares it.
hook_axial = [dist(vid(i, j), vid(i + 1, j))
              for i in range(N_BASAL, N_BASAL + N_HOOK) for j in range(NCIRC)]
L_MATURE = max(hook_axial) * MATURE_FACTOR

# --- cells -------------------------------------------------------------------
cells = []
for i in range(NAX):
    if N_BASAL <= i < N_BASAL + N_HOOK:
        u = (i - N_BASAL + 0.5) / N_HOOK
        hook_ax = math.sin(math.pi * u) ** 0.3   # broad plateau over the hook
        phi_deg = u * BEND_DEG
    else:
        hook_ax, phi_deg = 0.0, 0.0
    for j in range(NCIRC):
        th = 2.0 * math.pi * (j + 0.5) / NCIRC
        inner_w = max(0.0, math.cos(th))         # 1 inner, 0 outer
        auxin = AUXIN_BASE + (AUXIN_PEAK - AUXIN_BASE) * hook_ax * inner_w ** 1.5
        acid = ACID_DARK_OUTER - (ACID_DARK_OUTER - ACID_DARK_INNER) * hook_ax * inner_w
        cells.append((auxin, acid, 0, math.degrees(th)))
cells.append((AUXIN_BASE, ACID_DARK_OUTER, 1, 0.0))
cells.append((AUXIN_BASE, ACID_DARK_OUTER, 1, 0.0))

# --- write --------------------------------------------------------------------
with open(OUT, "w") as f:
    f.write(f"{NCELLS} {len(walls)} {len(verts)}\n")
    for w, (c1, c2, v1, v2, *_r) in enumerate(walls):
        f.write(f"{w} {c1} {c2} {v1} {v2}\n")
    f.write(f"\n{len(verts)} 3\n")
    for (x, y, z) in verts:
        f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
    f.write(f"\n{len(walls)} 1 3\n")
    for (c1, c2, v1, v2, grow, bendf) in walls:
        f.write(f"{dist(v1, v2):.6f} {grow} {bendf} {L_MATURE if grow else 0:.6f}\n")
    f.write(f"\n{NCELLS} 9\n")
    for (auxin, acid, ctype, th) in cells:
        f.write(f"0 0 0 0 0 {auxin:.4f} {acid:.4f} {ctype} {th:.2f}\n")

mid = N_BASAL + N_HOOK // 2
li = dist(vid(mid, 0), vid(mid + 1, 0))
lo = dist(vid(mid, NCIRC // 2), vid(mid + 1, NCIRC // 2))
print(f"wrote {OUT}: {NCELLS} cells ({NQUAD} quads + 2 caps), "
      f"{len(walls)} walls, {len(verts)} vertices")
print(f"  mid-hook axial cell length: inner {li:.1f} um, outer {lo:.1f} um")
print(f"  mature length L_mature = {L_MATURE:.1f} um "
      f"-> outer can extend {L_MATURE/lo:.2f}x, inner {L_MATURE/li:.2f}x")
print(f"  (measured: outer 1.06x, inner 2.15x at t=8 h)")
