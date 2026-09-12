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


## Performance

Measured on the 350-cell hook shell (Apple M1), full pipeline = 2.5 h
equilibration + 10 h light + 10 h dark (the latter two in parallel):

| configuration | full pipeline | per simulated hour |
|---|---|---|
| initial anisotropic version | 75 min | 390 s |
| + stress-refresh interval | 56 min | 292 s |
| + Hill-factor caching | 47 min | 244 s |
| + solver tolerance 1e-4 -> 1e-3 | **94 s** | **7 s** |

**48x faster, with a bit-for-bit unchanged trajectory.** The three changes:

1. **Stress-state refresh interval** (`VertexFromTRBScenterTriangulation`
   parameter 3). Profiling showed 37% of runtime in the per-step stress
   tensor + Jacobi pass. CMT reorientation is an hours-scale process, so the
   material is refreshed every 0.02 h instead of every step (-25%).
2. **Eq. 3 Hill factor cached per cell** instead of recomputed per wall-cell
   pair, hoisting two `pow()` calls out of the inner loop (-16%).
3. **Solver tolerance.** The dominant cost: `eps=1e-4` forced tiny steps to
   resolve elastic transients that the growth-driven trajectory does not
   depend on — the system is strongly attracted to force balance, so the
   slow manifold is insensitive to them. Validated by comparing 1e-4, 1e-3
   and 1e-2 over the full 10 h run: identical hook angles (162.7 -> 116.1
   deg) and inner fold (1.759x) in all three. 1e-3 is the shipped default
   (1e-2 halves runtime again if needed).

What does **not** help at this scale:

- **Threading**: 122 s (1 thread) vs 123 s (8 threads). With 350 cells and
  708 walls each parallel region is a few microseconds of work, comparable
  to the synchronization cost. The thread pool's grain thresholds correctly
  keep these loops serial; threading pays from ~10k cells upward (the 40k
  cell benchmark in ../README.md). Run independent conditions (light/dark,
  parameter sweeps) as separate processes instead - that is what
  `pipeline3d.py` does.
- **Larger max step** (`h_max` 0.02 -> 0.1): no change, error control
  already chooses smaller steps.

Remaining opportunity: a **quasi-static solver** (Newton or FIRE
minimization of the elastic energy between growth steps, as the paper's own
Python model used). Explicit integration costs steps proportional to
Y/growth-rate, which is exactly the ratio that has to be large for the
quasi-static approximation to hold - so the cost is structural, not a
constant factor. That would also remove the need to scale Y and P together
for mobility. Estimated further 5-20x, and it would let the mesh be refined
without a step-size penalty.

## Visualising the material anisotropy

`pvhook3d_aniso.py` builds `hook3d_anisotropy.pvsm`:

```sh
../build/simulator hook3d.model hook3d_eq.init solver3d_vtk.rk5 \
    -centerTri_init -vtk_output vtk3d
/Applications/ParaView-5.10.0.app/Contents/bin/pvpython pvhook3d_aniso.py 10
open -a ParaView-5.10.0 hook3d_anisotropy.pvsm
```

The shell is coloured by stress anisotropy a = 1 - s2/s1, and a white line
glyph at each cell centre lies along that cell's maximal principal stress
direction - the axis the fibre stiffness is redistributed onto by Eq. 3,
i.e. the model's prediction of the CMT/cellulose orientation. Glyph length
scales with a (isotropic cells show a dot, strongly anisotropic cells a
bar).

In the closed hook the bars wrap around the tube (circumferential, matching
dark CMT arrays). During opening the outer-flank bars rotate to run along
the tube while the inner side stays circumferential and its anisotropy
*rises* - the inner/outer switch asymmetry of Figs 3E-F and 4E.
