# Apical hook opening — cell-resolved 2D model

A mechanistic, cell-resolved implementation of apical hook opening in the
tissue v2 simulator, based on:

> Walia A.†, **Carter R.**†, Wightman R., Meyerowitz E.M., Jönsson H.,
> Jones A.M. (2024) *Differential growth is an emergent property of
> mechanochemical feedback mechanisms in curved plant organs.*
> Developmental Cell 59:3245-3258. doi:10.1016/j.devcel.2024.09.021
> (PDF in this directory; model code: zenodo.org/records/13379829)

The paper models the hook as a toroidal pressure vessel with two ODEs (inner
and outer side length). This directory re-embodies its mechanism in a
**2D vertex model of the median longitudinal section** with realistic
geometry and anatomy, so that differential growth, tissue stresses and organ
shape emerge from cell-level rules.

## Geometry (make_hook_init.py)

- Hypocotyl radius r = 75.8 um and hook curvature R = 2.88 r: the paper's
  fitted values; closed hook angle ~156 deg after turgor equilibration
  (experimental t0: 158.7 deg).
- Anatomy across the diameter (10 cell files): epidermis (10 um), two cortex
  layers (22 + 22 um), endodermis (12 um), stele (2 x 9.8 um), mirrored
  (cf. Gendreau et al. 1997 for hypocotyl anatomy).
- 43 axial slices (basal straight, 32 hook slices, apical straight);
  430 cells, 913 walls. Inner hook epidermal cells ~13 um, outer ~26 um at
  start (geometry-given, matching the measured inner:outer asymmetry).
- Auxin (cell variable 4): baseline 0.6 everywhere in darkness (transiting
  auxin represses growth = maintenance phase); DR5-like maximum (1.0) on the
  inner hook side with a wide plateau. Illumination at t=0 starts first-order
  decay (k = 0.9/h).

## Mechanism (hook.model), mapped from the paper

| paper ingredient | implementation |
|---|---|
| uniform turgor loading | `Pressure2D::AreaPotential` (P on every cell) |
| thick outer epidermal wall = growth constraint | outer-surface walls 33x stiffer (`WallMechanics::Spring` K_force2, stiffFlag); inner exterior wall soft (it thins during opening, Fig S2) |
| growth saturation at maximal cell length | `WallGrowth::StrainSaturationHill` L_max: outer epidermis frozen (extends 1.07x in experiments), inner epidermis capped at 40 um (3x; Fig 1 fold change) |
| subepidermal longitudinal force sL(t) | active auxin-gated elongation of interior files (`strain_flag=0` mode), gated to stop against compression |
| auxin gates growth (YUC6-OX blocks opening) | repressing Hill (K=0.3, n=6) on every growth term |
| light-triggered auxin depletion | `Degradation::One` on the auxin variable (t=0 = illumination) |
| tissue stresses (qua2 cracks: epidermal tension) | emerges: outer wall tension >> interior; interior compressed |

Numerical stabilizers (new v2 reactions, both physically motivated):
`CenterTriangulation::EdgeSpring` + `EdgeRelaxation` make cell interiors
viscoelastic (instantaneously stiff, fluid over ~10 min — suppresses the
quad-mesh shear mode), and `WallMechanics::BendingChain` gives wall files
plate-like bending stiffness (suppresses wrinkling without resisting
organ-scale shear).

## Running

```sh
python3 make_hook_init.py hook.init   # generate geometry
python3 pipeline.py                   # equilibrate (dark), then light + dark runs
python3 analyze_hook.py run_light.out # hook angle + arc lengths vs time
# animation (VTK + ParaView):
../build/simulator hook.model hook_eq.init solver_anim.rk5 -vtk_output vtk
open -a ParaView-5.10.0 hook_opening.pvsm
```

The pipeline first relaxes the tissue to its pressurized dark steady state
(growth frozen), the biological maintenance phase, then runs illumination
(hook.model) and dark control (hook_dark.model) from that state.

## Results (this parameterization)

| t (h) | sim light | sim dark | experiment (light) |
|---|---|---|---|
| 0 | 156 | 156 | 159 |
| 2 | 140 | 146 | 133 |
| 4 | 119 | 139 | 82 |
| 6 | 98 | 132 | 53 |
| 8 | 80 | 126 | 44 |
| 10 | 66 | 122 | 35 |

- Differential growth emerges: at t=10 the inner epidermal arc has extended
  1.44x, the outer 1.05x (experiment: inner 2-4x at full opening, outer
  ~1.07x). Opening completes (<40 deg) by t~16 h.
- Dark control: hook maintained, slow drift (matches the slow dark opening
  reported in the paper's Fig 5 scenario).
- Tissue-stress state matches qua2 observations: epidermis under tension
  (outer wall carrying ~30x interior levels), interior under compression.

## Known limitation

Opening runs ~1.5-2x slower than experiment in mid-phase. A 2D median
section cannot represent the hoop-stress asymmetry of the 3D toroid
(higher circumferential stress on the inner side, Walia et al. Fig 3B),
whose in-plane consequence is an artificial net closing moment carried by
the outer wall tension. This is precisely why the paper's own quantitative
model uses the analytic toroidal formulas; a 3D shell version (legacy
`hypocotyl3D` TRBS reactions, not yet ported) would remove the discrepancy.
The `Pressure2D::CapForce` reaction documents why the naive 2D fix fails
for closed hooks.


## 3D TRBS shell model (hook3d)

The definitive version: the hypocotyl epidermis as a closed, turgid 3D
surface (720-cell tube + end caps at full resolution; production mesh 350
cells), with proper 2D wall elasticity via triangular biquadratic springs
(`VertexFromTRBScenterTriangulation`, exact port of the legacy TRBS force
kernel) and shell-normal turgor (`Pressure3D::CenterTriangulation`). Hoop
stress, the pressure cap force and the toroidal stress distribution are real
physics here, so the hook opens the way the organ does: the inner side
expands longitudinally and the arms follow.

Mechanism (hook3d.model): uniform turgor loads the shell; axial walls and
internal edges yield above a strain threshold (Lockhart), saturating at a
maximal length (outer hook cells are born near saturation - the
geometry-given differential); auxin (high on the inner side in darkness)
represses yielding; illumination at t=0 depletes it.

Two numerical lessons encoded here:
1. The overdamped per-vertex drag of the framework suppresses organ-scale
   rotations (the arm swing) by orders of magnitude - the reason the paper's
   own model used quasi-static Newton-Raphson. Verified by a mobility test
   (the grown state unrolls at 7 deg/h at 50x force scaling, 0.14 deg/h at
   baseline). Fixed by scaling Y and P together (identical strains, faster
   mechanics) and keeping the arms short.
2. TRBS needs guarded acos/Heron evaluations so bad adaptive trial steps are
   rejected by error control instead of NaN-poisoning the state.

Results (`pipeline3d.py`; dark-equilibrated start at 161 deg):

| t (h) | sim light | sim dark | experiment (light) |
|---|---|---|---|
| 0 | 161 | 161 | 159 |
| 2 | 155 | 160 | 133 |
| 5 | 104 | 161 | 63 |
| 8 | 63 | 161 | 44 |
| 10 | 43 | 162 | 35 |

Inner hook arc extends 2.08x (exp 2-4x), outer 1.13x (exp 1.07x); the dark
control is fully maintained. The ~1.5 h mid-phase lag tracks the auxin
clearance time constant (k_d, k_growth are the tuning knobs).
