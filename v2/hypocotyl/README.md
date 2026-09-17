# Apical hook opening — a cell-resolved 3D model

A biologically parameterised model of *Arabidopsis* apical hook opening,
built on the Tissue v2 solver. The hypocotyl epidermis is represented as a
closed, turgid, centre-triangulated 3D shell whose walls have direction-
dependent stiffness set by the local stress state — the mechanism proposed in

> A. Walia, R. Carter, *et al.* (2024) "Dynamic cell wall anisotropy governs
> differential growth during apical hook development in *Arabidopsis*."
> *Developmental Cell* **59**:3245–3259.

Given a dark-grown hook and a light stimulus at t = 0, the model reproduces
the measured opening kinetics (RMSE 11.9° over 0-10 h), the inner-flank
tissue-length fold change (2.16x against 2.149x measured; the outer flank
overshoots, 1.11x against 1.055x), the inner/outer microtubule reorientation
asymmetry (**but see the recalibration section — that one turns out to be an
artefact of the solver's drag lag, not a property of the growing
equilibrium**), and the direction of the auxin- and pH-dependent
perturbations
(though those act on the growth gate directly, so their direction is
guaranteed by construction and only their magnitude is a test — see below).
It does **not** reproduce the microtubule and cellulose perturbations — the
fibre stiffness turns out to carry too little load for wall anisotropy to
drive anything. That discrepancy, and what to do about it, is documented
below.


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
| `solver_quasistatic.rk5` | experimental: equilibrium-solving solver (see below) |
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
`Lmax`, set to 1.05× the longest (outer-flank) cell. That one assumption, with
no free inner/outer parameter, gives each flank a ceiling close to what is
measured:

| | ceiling set by the geometry | measured (t = 8 h) | reached in simulation |
|---|---|---|---|
| outer flank (**fitted**) | 1.05× | 1.055× | 1.11× |
| inner flank (**predicted**) | 2.17× | 2.149× | 2.10× (2.16× by t = 10 h) |

Which row is which matters. `MATURE_FACTOR` is chosen to put the outer ceiling
at the measured 1.055×, so that row is a fit of one parameter to one number.
The inner ceiling is then whatever the geometry gives — 2.17× against 2.149×
measured — with nothing further to adjust. That is the actual content of the
"shared mature length" claim: one fit buys one independent prediction.

The first column is a property of the initial file, not a result: it is what
the generated rest lengths allow, printed by `make_init.py`. Only the third
column is a simulation outcome, and for the outer flank it sits ~6% above the
ceiling, at 1.11× against 1.055× measured. Two reasons, neither of which is
walls growing past `Lmax` (they cannot — `WallGrowth::AcidGrowth` stops
outright once `1 - L/Lmax` reaches zero):

- the observable is a **deformed** arc length while `Lmax` bounds the **rest**
  length, and the mean wall elastic strain rises from 3.0% at t = 0 to 6.0% by
  t = 8 h as turgor is expressed against softening walls — about half the gap;
- `Lmax` is set from the *longest* outer cell, so shorter cells elsewhere in
  the sector have headroom above 1.05× of their own initial length.

The inner flank tracks its target closely either way. Quoted as three separate
columns because the geometric construction and the simulated outcome are
different claims, and conflating them would make the agreement look better
than it is.

Chemistry initialised from Figs 2 and S3: a DR5 auxin maximum on the inner
flank (1.00 vs 0.45 baseline), and apoplastic pH with the inner flank more
alkaline than the outer in darkness (acidification 0.05 vs 0.25).


## Parameter provenance

How many knobs bought the agreement, and which. Anything marked *fitted* was
adjusted against a model output and cannot also be evidence for it.

| parameter | source |
|---|---|
| r = 75.8 µm, R/r = 2.88 | **measured** — the paper's fitted toroid |
| 16 epidermal files | **measured** — *Arabidopsis* anatomy (~16–24) |
| k = 0.7, n = 1.2 (Eq. 3) | **measured** — the paper's values verbatim |
| Y_m : Y_f = 75 : 100 | **measured** — the paper's stiffness ratio |
| auxin 0.45 / 1.00, pH 0.25 / 0.05 | **measured** — Figs 2 and S3 |
| k_d = 0.5 (τ = 2 h) | **measured** — the hook barely moves for the first ~2 h |
| acid drive/relax = 2.0 (τ = 0.5 h) | **measured** — both flanks acidify within 30 min |
| **k_growth = 80** | **fitted** — set so opening runs on the measured timescale |
| **MATURE_FACTOR = 1.05** | **fitted** — set to the measured outer fold change |
| **Y, Y_f, P absolute scale (×25)** | **fitted in effect** — see below |
| Poisson 0.3, yield strain 0.005 | chosen, standard values |
| K_auxin 0.55, n 6; K_acid 0.45, n 3 | chosen by hand; no independent measurement, and not demonstrably held fixed during development |
| CT internal-edge L_max = 30 | chosen |
| generated bend 172°, apex taper 0.86 | chosen so the relaxed shape starts near the measured 159° |
| stress refresh 0.02 h | numerical; measured to move the hook angle ≤2×10⁻³° vs refreshing every step |
| solver tolerance 1e-2 | numerical; measured to agree with 1e-3 to 1.4×10⁻⁴° on this model |

The third fitted entry deserves its name. Scaling Y, Y_fiber and P together
leaves the *equilibrium* shape untouched, which is why it was treated as a free
numerical choice — but under `RK5Adaptive` the shell never reaches equilibrium,
and the scale sets how far it lags (see the quasi-static section). The lag sets
the opening rate. So the ×25 is a second knob on the same observable
`k_growth` was fitted to, and the two are not independent.

That is three effective degrees of freedom behind the opening curve, one of
which (MATURE_FACTOR) also buys the independent inner-flank prediction. The
neck of the argument is that nothing in the list was fitted to the microtubule
reorientation, the anisotropy magnitudes, or the perturbation responses.


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
6 h — but read that as fit quality, not as a prediction. Three effective
degrees of freedom sit behind this curve (see *Parameter provenance*):
`k_growth`, `MATURE_FACTOR`, and the ×25 scaling of Y and P, which is
shape-neutral at equilibrium but sets the drag lag that partly determines the
opening rate under `RK5Adaptive`. What is *not* fitted to this curve, and so
carries the evidential weight, is everything else: the inner-flank fold
change, the microtubule reorientation, the anisotropy magnitudes, and the
perturbation responses. (Not being fitted is necessary for a result to count
as evidence, not sufficient — the microtubule reorientation is unfitted and
still turned out to be a solver artefact.)

The model over-opens late: real hooks stall near 35–45° once the
cotyledons separate, which this model has no representation of (the paper's
own `a2` replicate series bottoms out at 22.7°, so part of the late spread is
experimental). The dark control never opens — but that is a consistency
check, not evidence: the dark variant zeroes auxin depletion and both
acidification terms, so every gate on growth is shut by construction and the
hook could not open. It confirms the variant machinery does what it claims
(the `@TAG` substitution failure below is exactly what it is there to catch),
and nothing about the biology.

**Microtubule reorientation.** MT angle is read out as the predicted CMT axis,
i.e. the maximal principal stress direction (90° = circumferential/transverse,
0° = longitudinal). The inner flank holds 89–90° circumferential throughout
while the outer flank switches to longitudinal within the first hour and
rotates back as the organ straightens — the inner/outer switch asymmetry of
Figs 3E–F and 4E, emergent rather than imposed. **This result does not
survive removal of the drag lag**: at force balance the outer flank stays
circumferential throughout, so the switch is a transient of the relaxation
rather than a property of the growing equilibrium. See *Recalibrating to the
paper's stiffness ratio* for the matched-hook-angle control.

**Perturbations** (`run_perturbations.py`), each a single parameter change.
Hook angle in degrees; WT light is 79.2° at 4 h and 12.5° at 10 h. The two
groups have to be read differently, because only one of them is a test:

| treatment | what it changes | 4 h | 10 h | separation from control | measured separation |
|---|---|---|---|---|---|
| YUC6-OX | auxin gate never released (k_d = 0) | 121.0 | 53.8 | 41.3° | 78.7° |
| low light | slower auxin depletion, weaker acidification | 121.2 | 29.2 | 16.7° | 11.6° |
| oryzalin | MTs depolymerised → isotropic fibre (K_hill → 50) | 79.2 | 12.5 | **0.0°** | blocked |
| isoxaben 100 nM | less cellulose (Y_f 1350 → 470) | 78.7 | 11.8 | −0.7° | blocked |
| isoxaben 600 nM | Y_f → 200 | 78.6 | 11.5 | −1.0° | blocked |

The first two act **directly on the growth gate**: with auxin never depleted,
`G_auxin` stays low and growth is suppressed, so a slowdown is guaranteed by
construction and its *direction* tests nothing. Only the magnitude is a
result, and there the model under-blocks — YUC6-OX separates from WT by 41°
where the measurement separates by 79°, about half. Low light is closer
(17° against 12°), though the model's absolute angles run low throughout
because it over-opens late.

The last three act on the **wall fibre** and are genuine, failed predictions:
oryzalin is bit-identical to wild type, and isoxaben moves the angle by ~1° in
the *wrong* direction. Those are the ones that matter, and they fail for the
reason set out next.


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

**This is not an artefact of the drag lag.** Repeating the control under
`QuasiStatic`, where the shell is at force balance every step and there is no
lag by construction, gives the same answer — 4.0° at peak against 3.4° under
RK5Adaptive, about 4.5% of the opening either way:

| t (h) | 0.75 | 1.00 | 1.50 |
|---|---|---|---|
| Y_f = 1350 | 113.3 | 69.3 | 15.6 |
| Y_f = 0 | 110.2 | 65.3 | 13.6 |

The explanation is simply how little of the force balance the fibre carries:
the fibre carries **2.67%** of it. That figure is measured, not inferred from
the parameters: relaxing the same geometry with and without the fibre and
comparing mean wall strain gives `K_f/K_m = eps(Y_f=0)/eps(Y_f) - 1`, which
needs no assumption that the two force laws share units. (Dividing the raw
parameters, Y_f/(Y_m+Y_f) = 6.3%, overstates it 2.4-fold — they are not
comparable, since the fibre term carries an extra `area/(2d)` width factor.)
A strain-gated growth law then turns a 2.7% stiffness change into a
few-percent growth change. A four-point dose-response under `QuasiStatic` confirms the gain is
linear in the fibre's share, measured at peak sensitivity (t = 1 h):

| Y_f | share of stiffness | effect (deg) | linear prediction | ratio |
|---|---|---|---|---|
| 1350 | 2.67% | 0.000 | – | – |
| 675 | 1.35% | 1.965 | 2.007 | 0.98 |
| 337.5 | 0.68% | 2.977 | 3.010 | 0.99 |
| 0 | 0% | 4.013 | 4.013 | 1.00 |

within 2% of proportional throughout. (These are re-runs with the relaxation
cap raised to 100000. The first version of this table came from runs that
silently hit a 3000 cap on some growth steps, i.e. were not fully relaxed;
the converged numbers are *cleaner*, moving the ratios from 0.90/0.97 to
0.98/0.99, but the earlier ones should not have been quoted. Always check the
solver's cap warning.)

**The quasi-static answer converges in the growth step.** That is not
automatic — a model whose growth law has no braking term can fail to have a
quasi-static limit at all — so it is worth demonstrating rather than assuming
(0-2 h, relaxation cap 100000, zero cap hits):

| h_growth | angle at 2 h | difference | inner fold |
|---|---|---|---|
| 0.02 | 7.215 | – | 2.1975 |
| 0.01 | 6.286 | −0.929 | 2.2044 |
| 0.005 | 5.917 | −0.369 | 2.2072 |
| 0.0025 | 5.752 | −0.165 | 2.2085 |

Successive differences shrink by a factor of ~0.42, extrapolating to ~5.6°
and an inner fold of ~2.209. The production runs use `RK5Adaptive` regardless;
this establishes that the `QuasiStatic` comparisons above rest on a limit that
exists. (An earlier hint of strong
non-linearity — isoxaben at Y_f = 470 appearing to cost only 0.5° — was an
artefact of comparing across different solvers and timepoints: the observable
saturates late, when the hook is nearly open, so *when* the effect is measured
matters more than it looks. Peak sensitivity is around t = 1 h.)

Note also that reducing Y_f cuts both the magnitude *and* the anisotropy of
the fibre term, and only the magnitude channel matters here: oryzalin
(`K_hill` -> 50, which makes the fibre isotropic at unchanged magnitude)
changes nothing at all, while Y_f -> 0 changes 4.0°. (Tested because a parallel model — pavement-cell lobing, whose
boundary is kinematically clamped — finds *exactly* zero for its analogous
stiffening channel, and predicted from that the effect should scale with the
softest mode's drag time. It does not: removing the lag entirely leaves the
effect intact. Note that scaling *one* of several competing force terms does
change the equilibrium; only scaling them all together leaves it invariant,
which is why the Y/P mobility rescaling is shape-neutral but zeroing Y_f is
not.)

The matrix modulus (Y_m = 20000) dominates the fibre term, which enters
`WallMechanics::FiberSpring` as `Y_f · orient · area/(2d)` with
`orient ∈ [0.5(1−g), 0.5(1+g)]`. Even full anisotropy therefore redistributes
a small fraction of the total wall stiffness, so Eq. 3 has no mechanical
consequence.

What this means for the results above:

- The **stress anisotropy is a faithful readout** of whatever stress state
  the solver is in — but the CMT switch it produces in `fig4` and `fig5` is a
  readout of the *lagged* state, and disappears at force balance. It is not
  the independent confirmation of Figs 3E–F it was presented as.
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


## Recalibrating to the paper's stiffness ratio

The section above says the fix is to move Y_m : Y_f toward the paper's
75 : 100 and refit. That was done. It works for one of the two failed
perturbation classes, cannot work for the other, and costs a result that had
been reported as a success. All three outcomes are worth having.

**The recipe.** Measuring the fibre share directly (above) gives
K_f/K_m = 0.0274 at Y_f = 1350, and the response is linear in Y_f, so the
paper's 57% share needs **Y_f = 65700**. Raising it that far stiffens the
shell, which would change the growth drive as well as its distribution, so
turgor is raised to **P = 36.2** to put mean wall strain back on its original
value (0.0304 against 0.0306, within 0.7%). Only the matrix:fibre *split* then
differs from the shipped model. `k_growth` refits to **13** under `QuasiStatic`
(scanned over 8/13/20 against the measured curve; at t = 1 h the model gives
155.8 deg against 155.6 measured, and at t = 1.5 h, 143.7 against ~144.5).

The quasi-static solver is not optional here. The fibre loads the *stiffest*
mode disproportionately: lambda_max rises from 1.5e5 to 1.4e6, a factor of 9,
where mean stiffness rises only 2.2x. Explicit integration costs h ~ 1/lambda
and so becomes 9x more expensive, while FIRE costs dt ~ 1/sqrt(lambda), 3x.

**What it fixes: oryzalin.** Depolymerising microtubules is modelled by
K_hill -> 50, which makes the fibre isotropic at unchanged magnitude. With the
fibre carrying 57% rather than 2.7%, that redistribution finally has something
to redistribute. Opening measured from each run's own frame 0 (QuasiStatic
relaxes before its first print, so the variants' starting shapes differ and
absolute angles are not comparable):

| t (h) | WT | oryzalin | isoxaben 100 nM | isoxaben 600 nM |
|---|---|---|---|---|
| 0.75 | 2.1 | 1.2 | 4.5 | 7.4 |
| 1.50 | 17.1 | 13.8 | 38.4 | 52.1 |
| 3.00 | 61.0 | 51.8 | 111.3 | 129.1 |
| 4.00 | 81.7 | 69.3 | 132.1 | 144.2 |

(degrees opened; all four runs converged, zero relaxation-cap hits). Oryzalin
now opens **15% slower** than wild type, monotone from t = 0.75 h onward,
against exactly 0.0 difference before recalibration. The direction matches
experiment.

**What it cannot fix: isoxaben.** Reducing cellulose is modelled by lowering
Y_f, and in a Lockhart model a softer wall yields *faster* - so opening
accelerates, by 62% and 77% at t = 4 h (and by more earlier), where the experiment shows it
blocked. Recalibration makes this worse rather than better, and no choice of
stiffness can reverse it: the sign is forced by the growth law. The real
mechanism is that cellulose-synthesis inhibition starves new wall deposition,
which is a limit on *synthesis*, not on *compliance*, and this model has no
wall-synthesis term at all. Reproducing isoxaben needs a new term, not a new
parameter.

**What it costs: the microtubule reorientation was a lag artefact.** The
outer-flank switch from circumferential to longitudinal, reported above as
emergent and matching Figs 3E-F, does not survive removal of the drag lag.
Compared at *matched hook angle* rather than matched time, so the three runs
are compared at the same shape:

| outer-flank MT angle | hook 150 | hook 140 | hook 130 | hook 120 |
|---|---|---|---|---|
| RK5, Y_f = 1350 (shipped) | 3 | 5 | 6 | 13 |
| QuasiStatic, Y_f = 1350 | 87 | 88 | 88 | 88 |
| QuasiStatic, Y_f = 65700 | 57 | 73 | 77 | 81 |

(90 = circumferential, 0 = longitudinal.) At force balance the outer flank
stays circumferential at every stage of opening. The switch appears only when
the shell lags: the arms swing behind where the rest lengths put them, which
loads the outer flank longitudinally as a *transient of the relaxation*, not
as a property of the growing equilibrium. Recalibrating recovers a little of
the rotation (57 deg at hook 150) but nothing resembling a switch.

This matters beyond the model. A thin-walled pressure vessel has hoop stress
twice longitudinal, so maximal principal stress is circumferential unless
differential growth overturns it - and here it does not. If cortical
microtubules follow maximal principal stress, a quasi-static model does not
produce the measured outer-flank switch. Either the walls carry real
viscoelastic lag on the growth timescale, or microtubules follow something
other than the instantaneous stress direction - strain rate would be the
obvious candidate, and is not implemented here.

**Status: measured, not adopted.** The recalibrated configuration is not the
shipped model. It trades a working microtubule prediction for a working
oryzalin prediction, worsens isoxaben, and has been validated only to t = 4 h -
the fold changes, the dark control and the full 12 h trajectory have not been
re-run against it. The recipe above reproduces it in full
(Y_f = 65700, P = 36.2, k_growth = 13, `QuasiStatic`).


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

**48× faster, with the trajectory unchanged to 2×10⁻³ degrees.** Not
bit-for-bit, which is what this used to claim. Two of the three changes below
are exact refactors; the stress-refresh interval is not, and holding the
material constant for ~1000 solver steps instead of updating it every step
moves the hook angle by up to 2×10⁻³° over 0.05 h (measured against
`interval = 0`). Far below anything that matters here, but the justification
for the interval — that CMT reorientation is an hours-scale process — is an
argument about the *biology* being insensitive, not about the arithmetic being
identical, and it was quoted as the latter. The three changes:

1. **Stress-state refresh interval** (`VertexFromTRBScenterTriangulation`
   parameter 3). Profiling put 37% of runtime in the per-step stress tensor +
   Jacobi pass. CMT reorientation is an hours-scale process, so the material
   is refreshed every 0.02 h rather than every step (−25%).
2. **Eq. 3 Hill factor cached per cell** instead of recomputed per wall-cell
   pair, hoisting two `pow()` calls out of the inner loop (−16%).
3. **Solver tolerance.** The dominant cost: `eps=1e-4` forced tiny steps to
   resolve elastic transients the growth-driven trajectory does not depend
   on. Validated by comparing 1e-4, 1e-3 and 1e-2 over a full run of the
   350-cell model, and re-checked on the current 514-cell one: 1e-2 and 1e-3
   agree to 1.4×10⁻⁴ degrees, i.e. seven significant figures, not "identical".

What does **not** help at this scale:

- **Threading**: 122 s (1 thread) vs 123 s (8 threads). With ~500 cells each
  parallel region is a few microseconds, comparable to synchronisation cost.
  The thread pool's grain thresholds correctly keep these loops serial;
  threading pays from ~10k cells upward (see `../README.md`). Run independent
  conditions as separate processes instead — that is what `pipeline.py` and
  `run_perturbations.py` do.
- **Larger `h_max`** (0.02 → 0.1): no change, error control already chooses
  smaller steps.

### The quasi-static solver, and what it revealed

`QuasiStatic` (v2 solver, `solver_quasistatic.rk5`) does what this section
used to propose: it steps growth and then solves mechanical equilibrium with
FIRE relaxation instead of integrating toward it. Explicit integration is
bounded by the *stiffest* elastic mode (mean `h` here is 2e-5 h) while the
process being modelled runs on the *softest* one, a ratio near 1e5.

It works, and it is about 2x faster than `RK5Adaptive` on this model — but it
does **not** reproduce the trajectory, and that is the interesting part:

| t (h) | 0 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|
| RK5Adaptive | 159.2 | 170.5 | 159.9 | 139.1 | 121.5 |
| QuasiStatic | 159.0 | 149.0 | 69.3 | 15.6 | 6.3 |

Halving the growth step and tightening the force tolerance tenfold moves the
QuasiStatic column by under 3°, so this is not discretisation error. The two
solvers are answering different questions: with per-vertex drag, the shell
lags far behind force balance, and **the opening kinetics this model was
tuned to reproduce are substantially that lag rather than growth kinetics.**
Under true equilibrium the same `k_growth` opens the hook roughly six times
too fast.

That matters beyond performance. It is a second symptom of the mechanics
being mis-scaled (the first being the fibre load share above), and it points
the same way: `k_growth` needs recalibrating against equilibrium mechanics,
after which wall stiffness would actually influence the opening rate — which
is the precondition for the oryzalin and isoxaben predictions to work at all.
Until that recalibration is done, `solver.rk5` (RK5Adaptive) remains the one
the published numbers here come from.

Note that uniformly scaling Y, Y_fiber and P buys nothing under QuasiStatic
(62 s scaled vs 74 s at biological stiffness for the same 2 h): FIRE converges
in O(sqrt(condition number)) iterations and uniform scaling leaves the
condition number unchanged. What it buys is that the scaling is no longer
*needed*, since equilibrium does not depend on it.


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
