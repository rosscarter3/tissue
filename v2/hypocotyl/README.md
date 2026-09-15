# Apical hook opening — a cell-resolved 3D model

A biologically parameterised model of *Arabidopsis* apical hook opening,
built on the Tissue v2 solver. The hypocotyl epidermis is represented as a
closed, turgid, centre-triangulated 3D shell whose walls have direction-
dependent stiffness set by the local stress state — the mechanism proposed in

> A. Walia, R. Carter, *et al.* (2024) "Dynamic cell wall anisotropy governs
> differential growth during apical hook development in *Arabidopsis*."
> *Developmental Cell* **59**:3245–3259.

Given a dark-grown hook and a light stimulus at t = 0, the model reproduces
the measured opening kinetics (RMSE 11.9° over 0-10 h), the inner/outer
tissue-length fold changes, the inner/outer microtubule reorientation
asymmetry, and the auxin- and pH-dependent perturbations. It does **not**
reproduce the microtubule and cellulose perturbations — the fibre stiffness
turns out to carry too little load for wall anisotropy to drive anything.
That discrepancy, and what to do about it, is documented below.


## Quick start

```sh
python3 make_init.py                 # generate the geometry -> hook.init
python3 pipeline.py                  # equilibrate in the dark, then run light + dark
python3 analyze.py run_light.out     # observables vs experiment
python3 run_perturbations.py         # oryzalin / isoxaben / YUC6-OX / low light
python3 plot_compare.py 12           # all six comparison figures -> figures/
```

ParaView animation with the predicted CMT axes:

```sh
/Applications/ParaView-5.10.0.app/Contents/bin/pvpython paraview_state.py 12
open -a ParaView-5.10.0 hook.pvsm
```


## Files

| file | role |
|---|---|
| `make_init.py` | geometry + initial auxin/pH/wall fields → `hook.init` |
| `hook.model` | the ten reactions (mechanics, material, growth, chemistry) |
| `solver.rk5`, `solver_eq.rk5`, `solver_vtk.rk5` | production / equilibration / animation settings |
| `pipeline.py` | dark equilibration → parallel light + dark + VTK runs |
| `run_perturbations.py` | one parameter change per treatment, run in parallel |
| `analyze.py` | hook angle, arc fold change, MT angle, auxin, pH |
| `expdata.py` | the published measurements, parsed from the authors' package |
| `plot_compare.py` | six experiment-vs-simulation figures |
| `paraview_state.py` | ParaView state + renders |
| `refdata/walia2024/` | the authors' `data.py` + `length_measurements.csv` |

Generated files (`*.init`, `*.out`, `vtk/`, `*.pvsm`, `equil.model`,
`hook_dark.model`, `pt_*.model`) are not tracked.


## Geometry (`make_init.py`)

Anatomy from the paper's fitted toroid and from hypocotyl histology:

- hypocotyl radius r = 75.8 µm, hook curvature R = 2.88 r (both fitted in the
  paper), generated bend 172° relaxing to ~159° under turgor
- 16 epidermal cell files around the circumference (*Arabidopsis* has ~16–24)
- 3 basal + 26 hook + 3 apical rings, with a gentle taper over the apical
  segment where the cotyledons attach
- 514 cells, 1040 walls, 528 vertices

The inner/outer length difference is not imposed — it is geometry. On the
curved segment the inner surface sits at radius R − r and the outer at R + r,
so inner cells are born short and outer cells long.

The **mature length is shared**: every hook epidermal cell has the same
`Lmax`, set to 1.05× the longest (outer-flank) cell. That one assumption,
with no free inner/outer parameter, predicts both measured fold changes:

| | model | measured (t = 8 h) |
|---|---|---|
| outer flank | 1.05× | 1.055× |
| inner flank | 2.17× | 2.149× |

Chemistry initialised from Figs 2 and S3: a DR5 auxin maximum on the inner
flank (1.00 vs 0.45 baseline), and apoplastic pH with the inner flank more
alkaline than the outer in darkness (acidification 0.05 vs 0.25).


## Mechanism (`hook.model`)

| reaction | role |
|---|---|
| `CenterTriangulation::Initiate` | centre vertex per cell |
| `VertexFromTRBScenterTriangulation` | wall elasticity (triangular biquadratic springs), writes the cell stress tensor |
| `WallMechanics::FiberSpring` | CMT-guided cellulose stiffness, redistributed by Eq. 3 |
| `Pressure3D::CenterTriangulation` | turgor normal to the shell |
| `WallGrowth::AcidGrowth` | acid-growth wall yielding |
| `CenterTriangulation::WallGrowth::StrainSaturationHill` | internal edges follow the surface |
| `Degradation::One` ×2, `Creation::Zero` | auxin depletion, apoplast acidification |
| `VertexNoUpdateFromIndex` | clamps the apical end-cap ring |

**Dynamic anisotropic material.** The wall's fibre stiffness Y_f is split
between the two principal stress directions by Eq. 3 of the paper,

    g(a) = a^n / ((1−a)^n k^n + a^n),   k = 0.7, n = 1.2

where a = 1 − σ₂/σ₁ is the local stress anisotropy. At a = 0 both directions
get Y_f/2; at a = 1 all of it lies along the maximal stress direction. The
matrix modulus Y_m is isotropic; the two are calibrated to the paper's
Y_m : Y_f = 75 : 100.

**Acid growth.** Illumination releases the auxin gate and acidifies the
apoplast; walls then yield above a Lockhart threshold and saturate at the
mature length:

    dL/dt = k (ε − ε_th)₊ · L · (1 − L/Lmax)₊ · G_auxin · G_acid
    G_auxin = K_a^n / (K_a^n + auxin^n)      auxin represses extension
    G_acid  = acid^m / (K_c^m + acid^m)      acidification permits it

Auxin decays with τ = 2 h after light (so the hook barely moves for the first
~2 h, as measured); acidification completes within 30 min at both flanks.


## Results

Dark-equilibrated start, light at t = 0, angles in degrees
(`python3 analyze.py run_light.out`):

| t (h) | sim angle | exp angle | sim inner | exp inner | sim outer | exp outer | MT inner | MT outer |
|---|---|---|---|---|---|---|---|---|
| 0 | 159.2 | 158.7 | 1.00 | 1.00 | 1.00 | 1.00 | 90 | 90 |
| 1 | 159.9 | 155.6 | 1.42 | – | 1.13 | – | 89 | 3 |
| 2 | 121.5 | 133.4 | 1.71 | 1.36 | 1.14 | 1.07 | 89 | 9 |
| 3 | 97.0 | 102.9 | 1.83 | – | 1.13 | – | 90 | 22 |
| 4 | 79.2 | 81.5 | 1.89 | 1.92 | 1.13 | 1.04 | 90 | 33 |
| 5 | 64.2 | 63.2 | 1.95 | – | 1.12 | – | 90 | 41 |
| 6 | 50.5 | 53.3 | 2.00 | 2.04 | 1.12 | 1.06 | 90 | 49 |
| 8 | 26.9 | 44.3 | 2.10 | 2.15 | 1.11 | 1.05 | 89 | 64 |
| 10 | 12.5 | 34.5 | 2.16 | – | 1.11 | – | 89 | 80 |

**Hook-angle RMSE 11.9°** over 0–10 h, tracking within ~4° through the first
6 h. The model over-opens late: real hooks stall near 35–45° once the
cotyledons separate, which this model has no representation of (the paper's
own `a2` replicate series bottoms out at 22.7°, so part of the late spread is
experimental). The dark control never opens.

**Microtubule reorientation.** MT angle is read out as the predicted CMT axis,
i.e. the maximal principal stress direction (90° = circumferential/transverse,
0° = longitudinal). The inner flank holds 89–90° circumferential throughout
while the outer flank switches to longitudinal within the first hour and
rotates back as the organ straightens — the inner/outer switch asymmetry of
Figs 3E–F and 4E, emergent rather than imposed.

**Perturbations** (`run_perturbations.py`), each a single parameter change.
Two of the five reproduce the experiment and three do not — see the
limitation below:

Hook angle in degrees; WT light is 79.2° at 4 h and 12.5° at 10 h.

| treatment | what it changes | 4 h | 10 h | vs experiment |
|---|---|---|---|---|
| YUC6-OX | auxin gate never released (k_d = 0) | 121.0 | 53.8 | ✓ opening blocked (exp plateaus ~113°) |
| low light | slower auxin depletion, weaker acidification | 121.2 | 29.2 | ✓ opening slowed |
| oryzalin | MTs depolymerised → isotropic fibre (K_hill → 50) | 79.2 | 12.5 | ✗ **identical to WT**; exp is blocked |
| isoxaben 100 nM | less cellulose (Y_f 1350 → 470) | 78.7 | 11.8 | ✗ no effect (marginally *faster*); exp is blocked |
| isoxaben 600 nM | Y_f → 200 | 78.6 | 11.5 | ✗ no effect; exp is blocked |

The treatments acting on the *chemical gates* work; the treatments acting on
the *wall fibre* do nothing at all.


## The main open problem: the fibre carries almost no load

The wall's fibre stiffness is mechanically negligible in this
parameterisation. Setting `Y_fiber` to **zero** — removing the CMT-guided
cellulose entirely — costs at most 3.4° of hook angle over a 4 h run in which
the organ opens by 80°:

| t (h) | 0 | 1 | 2 | 3 | 3.75 |
|---|---|---|---|---|---|
| Y_f = 1350 (WT) | 159.2 | 159.9 | 121.5 | 97.0 | 83.3 |
| Y_f = 0 | 159.2 | 157.8 | 118.1 | 93.7 | 79.9 |

So the fibre accounts for roughly 4% of the effect, and in the *opposite*
direction to the experiment: removing it makes the hook open slightly faster
(less stiffness, more strain, more yielding), whereas isoxaben blocks opening.

The matrix modulus (Y_m = 20000) dominates the fibre term, which enters
`WallMechanics::FiberSpring` as `Y_f · orient · area/(2d)` with
`orient ∈ [0.5(1−g), 0.5(1+g)]`. Even full anisotropy therefore redistributes
a small fraction of the total wall stiffness, so Eq. 3 has no mechanical
consequence.

What this means for the results above:

- The **stress anisotropy and CMT reorientation are a faithful readout** —
  the inner/outer switch asymmetry in `fig4` and `fig5` is a real prediction
  of the stress solver, and it matches Figs 3E–F and 4E.
- But anisotropy is **not a driver** of opening in this model. The
  differential growth comes entirely from the geometric `Lmax` differential
  and the auxin/pH gates. That is why oryzalin and isoxaben do nothing, and
  it is a weaker claim than the paper's.

The fix is a recalibration of Y_m : Y_f toward the paper's 75 : 100, i.e.
making the fibre carry most of the load rather than ~1/15 of it. A probe run
at `Y_f = 27000` (20x) does separate oryzalin from its own control — 162.9°
vs 161.0° at 1 h, against 0.2° of separation at `Y_f = 1350` — so the
direction is right. That is an early-time probe only (the stiffer wall makes
the ODE about 10x more expensive per unit simulated time, so the run had not
reached the phase where the flanks diverge); whether 20x is *enough* is still
open.

This is therefore a re-tuning exercise rather than a one-line parameter
change: Y_m, Y_f, turgor, the growth rate and the mobility scaling (see
numerical lesson 1) all trade against each other, and the current fit
(RMSE 11.9°) was tuned against the existing balance. It is the single most
valuable next step for this model.


## Figures

`plot_compare.py` writes six figures to `figures/`, each overlaying the
published measurements (from `refdata/walia2024/`) on the simulation:

1. `fig1_opening` — hook angle vs time, light and dark
2. `fig2_lengths` — inner/middle/outer cell length fold change
3. `fig3_chemistry` — auxin depletion and apoplastic acidification
4. `fig4_anisotropy` — stress anisotropy and MT angle, inner vs outer
5. `fig5_cmt_polar` — CMT orientation distribution around the circumference
6. `fig6_perturbations` — the five treatments above

`paraview_state.py` additionally renders the shell coloured by auxin with a
white line glyph at each cell centre along the predicted CMT axis, its length
scaled by the stress anisotropy. In the closed hook the bars wrap around the
tube; during opening the outer-flank bars rotate to run along it while the
inner side stays circumferential and its anisotropy *rises*.


## Why 3D

The first version of this model was a 2D median section. It does not work,
and the reason is instructive: a 2D sheet cannot carry hoop stress, so the
toroid's circumferential asymmetry (higher hoop stress on the inner side,
Fig. 3B) is absent and the outer wall's tension becomes an artificial net
*closing* moment. The hook hinged at the shoulders instead of unrolling. This
is why the paper's own quantitative model uses the analytic toroidal
formulas. The 3D shell makes hoop stress, the pressure cap force and the
toroidal stress distribution real physics, and the organ then opens the way
it does in life: the inner side extends longitudinally and the arms follow.


## Numerical lessons

1. **Overdamped drag suppresses organ-scale rotation.** The framework's
   per-vertex drag slows the arm swing by orders of magnitude — the reason
   the paper's model used quasi-static Newton–Raphson. Verified by a mobility
   test: the grown state unrolls at 7°/h at 50× force scaling and 0.14°/h at
   baseline. Fixed by scaling Y and P together (identical strains, faster
   mechanics).
2. **TRBS needs guarded `acos` and Heron evaluations** so that bad adaptive
   trial steps are rejected by error control rather than NaN-poisoning the
   state.
3. **Material feedback belongs in `update()`, not `derivs()`.** Recomputing
   the stress-driven anisotropy inside the derivative evaluation makes the
   right-hand side non-smooth and collapses the adaptive step size:
   34 min equilibration became 2 min when it moved to once per accepted step.
4. **Model-variant editing must fail loudly.** The dark-control variant was
   originally produced by string substitution on the model text; rewording a
   comment silently stopped the match, so "dark" runs were quietly depleting
   auxin. All tunable parameters now carry `@TAG` markers and the variant
   functions regex on those.


## Performance

The optimisation history below was measured during development on the
350-cell shell (Apple M1), full pipeline = equilibration + 10 h light + 10 h
dark, the latter two in parallel. The current 514-cell mesh is ~1.5x larger;
timings for it are in the table after.

| configuration | full pipeline | per simulated hour |
|---|---|---|
| initial anisotropic version | 75 min | 390 s |
| + stress-refresh interval | 56 min | 292 s |
| + Hill-factor caching | 47 min | 244 s |
| + solver tolerance 1e-4 → 1e-3 | **94 s** | **7 s** |

**48× faster with a bit-for-bit unchanged trajectory.** The three changes:

1. **Stress-state refresh interval** (`VertexFromTRBScenterTriangulation`
   parameter 3). Profiling put 37% of runtime in the per-step stress tensor +
   Jacobi pass. CMT reorientation is an hours-scale process, so the material
   is refreshed every 0.02 h rather than every step (−25%).
2. **Eq. 3 Hill factor cached per cell** instead of recomputed per wall-cell
   pair, hoisting two `pow()` calls out of the inner loop (−16%).
3. **Solver tolerance.** The dominant cost: `eps=1e-4` forced tiny steps to
   resolve elastic transients the growth-driven trajectory does not depend
   on. Validated by comparing 1e-4, 1e-3 and 1e-2 over a full run — identical
   hook angles and inner fold in all three.

What does **not** help at this scale:

- **Threading**: 122 s (1 thread) vs 123 s (8 threads). With ~500 cells each
  parallel region is a few microseconds, comparable to synchronisation cost.
  The thread pool's grain thresholds correctly keep these loops serial;
  threading pays from ~10k cells upward (see `../README.md`). Run independent
  conditions as separate processes instead — that is what `pipeline.py` and
  `run_perturbations.py` do.
- **Larger `h_max`** (0.02 → 0.1): no change, error control already chooses
  smaller steps.

Remaining opportunity: a **quasi-static solver** (Newton or FIRE minimisation
of the elastic energy between growth steps, as the paper's Python model
used). Explicit integration costs steps proportional to Y/growth-rate, which
is exactly the ratio that must be large for the quasi-static approximation to
hold — so the cost is structural, not a constant factor. That would also
remove the need to scale Y and P together for mobility, and would let the
mesh be refined without a step-size penalty.


## Known limitations

- **Late over-opening.** Cotyledon separation and the maintenance phase are
  not represented, so the model keeps opening past the ~35° experimental
  plateau.
- **Single mature length.** Inner and outer walls share one `Lmax`, so the
  differential vanishes exactly as the hook straightens. This is what makes
  the fold-change prediction parameter-free, but it also gives an exponential
  approach to straight; separate per-flank saturation lengths would trade
  that prediction for better late kinetics.
- **Dark control drifts the wrong way.** Experimental dark hooks are
  maintained and then open slowly (160° → 100° over 28 h). The simulated dark
  control instead closes slightly further (159° → ~180°): with auxin
  undepleted the inner flank is gated off, but the outer flank still creeps,
  which tightens the hook. A dark maintenance phase is not modelled.
- **Cell length vs tissue length.** `fig2`'s left panel plots the paper's
  single-cell fold change (inner reaches ~4x); the simulation curve there is
  a tissue-level arc measure and belongs with the middle panel
  (`length_measurements.csv`, inner 2.15x), which it matches to 0.01x.
- **Epidermis only.** Inner tissues are represented by turgor alone.
