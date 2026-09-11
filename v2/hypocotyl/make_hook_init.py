#!/usr/bin/env python3
"""Generate a tissue init file for a 2D median longitudinal section of an
Arabidopsis apical hook, with realistic geometry and anatomy.

Geometry and parameters follow Walia, Carter et al. (2024) Dev Cell 59:3245
(hypocotyl radius r = 75.8 um, hook curvature R = 2.88 r, closed hook angle
~160 deg) and standard Arabidopsis hypocotyl anatomy (epidermis, two cortex
layers, endodermis, central stele; e.g. Gendreau et al. 1997 Plant Physiol).

Layout: cells indexed cell(slice, file) = slice*NFILES + file.
  files (j):  0 = outer epidermis ... 9 = inner epidermis (facing hook center)
  slices (i): basal straight, hook arc (160 deg), apical straight

Wall variables: [restLength, force(save), stiffFlag, epiGrowFlag, intGrowFlag, bendFlag, transFlag]
  stiffFlag = 1 on the two outer epidermal surfaces (thick outer walls)
  epiGrowFlag = 1 on epidermal longitudinal walls (saturating growth)
  intGrowFlag = 1 on interior longitudinal walls (sustained growth)
Cell variables: [0,0,0,0, auxin, cellType, hookPhiDeg]
  auxin: high on the inner hook side (DR5 pattern), baseline elsewhere
  cellType: 0 epidermis, 1 cortex, 2 endodermis, 3 stele
"""
import math
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else "hook.init"

# --- dimensions (um) --------------------------------------------------------
R_HYP = 75.8                      # hypocotyl radius (fitted, Walia et al.)
LAMBDA_R = 2.88                   # hook curvature ratio R/r (fitted)
R_AXIS = LAMBDA_R * R_HYP         # hook axis radius of curvature (218 um)
BEND_DEG = 200.0                  # generated bend; relaxes to ~160 deg
                                  # (closed hook) under turgor equilibration
BASAL_LEN = 150.0                 # straight hypocotyl below the hook
APICAL_LEN = 60.0                 # straight segment toward cotyledons
N_BASAL, N_HOOK, N_APICAL = 8, 32, 3

# radial cell files, outer surface -> inner surface
FILE_WIDTHS = [10.0, 22.0, 22.0, 12.0, 9.8, 9.8, 12.0, 22.0, 22.0, 10.0]
FILE_TYPES = [0, 1, 1, 2, 3, 3, 2, 1, 1, 0]
NFILES = len(FILE_WIDTHS)
assert abs(sum(FILE_WIDTHS) - 2 * R_HYP) < 1e-9

# Auxin: substantial levels transit the whole dark hypocotyl (represses
# growth everywhere in darkness -> hook maintenance); the inner hook side
# carries the DR5 maximum (delays inner-epidermis release after light).
AUXIN_BASE = 0.6
AUXIN_MAX = 1.0
AUXIN_FILE_WEIGHT = {9: 1.0, 8: 0.75, 7: 0.5}

# --- axis path ---------------------------------------------------------------
bend = math.radians(BEND_DEG)
lines = []  # (point, normal-toward-center) per slice boundary
for k in range(N_BASAL + 1):
    y = -BASAL_LEN + BASAL_LEN * k / N_BASAL
    lines.append(((0.0, y), (1.0, 0.0)))
for k in range(1, N_HOOK + 1):
    phi = bend * k / N_HOOK
    p = (R_AXIS - R_AXIS * math.cos(phi), R_AXIS * math.sin(phi))
    n = (math.cos(phi), -math.sin(phi))
    lines.append((p, n))
end_p, end_n = lines[-1]
tangent = (math.sin(bend), math.cos(bend))  # tangent at end of arc
for k in range(1, N_APICAL + 1):
    s = APICAL_LEN * k / N_APICAL
    lines.append(((end_p[0] + tangent[0] * s, end_p[1] + tangent[1] * s), end_n))

NSLICES = len(lines) - 1
NCELLS = NSLICES * NFILES

# file boundary offsets from the axis (outer -R_HYP ... inner +R_HYP)
cum = [-R_HYP]
for w in FILE_WIDTHS:
    cum.append(cum[-1] + w)

# vertices: v(line l, boundary fb)
NV_PER_LINE = NFILES + 1
verts = []
for (p, n) in lines:
    for fb in range(NV_PER_LINE):
        verts.append((p[0] + cum[fb] * n[0], p[1] + cum[fb] * n[1]))

def vid(l, fb):
    return l * NV_PER_LINE + fb

def cid(i, j):
    return i * NFILES + j

def dist(a, b):
    return math.hypot(verts[a][0] - verts[b][0], verts[a][1] - verts[b][1])

# --- walls -------------------------------------------------------------------
# (cell1, cell2, v1, v2, stiffFlag, epiGrowFlag, intGrowFlag)
walls = []
# transverse (radial) walls at each slice boundary line l, file j
for l in range(NSLICES + 1):
    for j in range(NFILES):
        below = cid(l - 1, j) if l > 0 else -1
        above = cid(l, j) if l < NSLICES else -1
        walls.append((above, below, vid(l, j), vid(l, j + 1), 0, 0, 0, 0, 1))
# longitudinal walls at file boundary fb, slice i
for fb in range(NV_PER_LINE):
    # thick, stiff outer wall on the OUTER hook side only; the inner
    # exterior wall thins/softens during opening (Walia et al. Fig S2)
    stiff = 1 if fb == 0 else 0
    # outer epidermis (fb 0,1): saturated, no growth (extends ~1.07x in
    # experiments); inner epidermis (fb 9,10): stress-driven growth; interior
    # (fb 2..8): active elongation engine
    if fb in (0, 1):
        epi, inter = 0, 0
    elif fb in (NFILES - 1, NFILES):
        epi, inter = 1, 0
    else:
        epi, inter = 0, 1
    for i in range(NSLICES):
        left = cid(i, fb - 1) if fb > 0 else -1
        right = cid(i, fb) if fb < NFILES else -1
        walls.append((left, right, vid(i, fb), vid(i + 1, fb), stiff, epi,
                      inter, 1, 0))

# --- cell variables -----------------------------------------------------------
cells = []
for i in range(NSLICES):
    if N_BASAL <= i < N_BASAL + N_HOOK:
        phi_deg = (i - N_BASAL + 0.5) / N_HOOK * BEND_DEG
        # wide plateau: the DR5 domain covers the whole inner hook side
        hook_profile = math.sin(math.pi * (i - N_BASAL + 0.5) / N_HOOK) ** 0.25
    else:
        phi_deg = 0.0
        hook_profile = 0.0
    for j in range(NFILES):
        auxin = AUXIN_BASE + (AUXIN_MAX - AUXIN_BASE) * hook_profile * \
            AUXIN_FILE_WEIGHT.get(j, 0.0)
        cells.append((auxin, FILE_TYPES[j], phi_deg))

# --- write init ----------------------------------------------------------------
with open(OUT, "w") as f:
    f.write(f"{NCELLS} {len(walls)} {len(verts)}\n")
    for w, (c1, c2, v1, v2, *_) in enumerate(walls):
        f.write(f"{w} {c1} {c2} {v1} {v2}\n")
    f.write(f"\n{len(verts)} 2\n")
    for (x, y) in verts:
        f.write(f"{x:.6f} {y:.6f}\n")
    f.write(f"\n{len(walls)} 1 6\n")
    for (c1, c2, v1, v2, stiff, epi, inter, bend, trans) in walls:
        f.write(f"{dist(v1, v2):.6f} 0.0 {stiff} {epi} {inter} {bend} "
                f"{trans}\n")
    f.write(f"\n{NCELLS} 7\n")
    for (auxin, ctype, phi) in cells:
        f.write(f"0 0 0 0 {auxin:.4f} {ctype} {phi:.2f}\n")

inner_hook_len = dist(vid(N_BASAL, NFILES - 1), vid(N_BASAL + 1, NFILES - 1))
outer_hook_len = dist(vid(N_BASAL, 0), vid(N_BASAL + 1, 0))
print(f"wrote {OUT}: {NCELLS} cells, {len(walls)} walls, {len(verts)} vertices")
print(f"hook: {N_HOOK} slices over {BEND_DEG} deg; inner epidermal cell "
      f"~{inner_hook_len:.1f} um, outer ~{outer_hook_len:.1f} um")
