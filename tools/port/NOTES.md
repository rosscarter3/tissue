# Reaction port: method and validated batches

Hand-maintained companion to the generated `STATUS.md`.

## Method

Ported reactions are checked against the behaviour they replace, not merely
compiled:

```sh
python3 tools/port/compare.py MODEL INIT SOLVER [--tol=1e-9]
```

drives both `bin/simulator` (legacy) and `build/simulator` over the same input
and diffs the numbers. It forces print flag 0 so models that natively emit VTK
or gnuplot still produce something comparable - without that a tutorial can
compare zero values and trivially pass.

Two classes of expected mismatch, both in the README under *Deliberate fixes
over legacy*: models with **wall** dynamics under `RK5Adaptive`, which legacy
mis-integrates, and anything seeded differently. Compare those under `Euler`
instead. Everything else should match to machine precision.

Fixtures live in `tests/port/`.

## Validated batches

### `legacy/degradation.cc` - 6 classes, 2026-09-17

`Degradation::Hill`, `HillN`, `TwoGeometric`, `OneWall`, `OneBoundary`,
`OneFromList`, each with its no-colon legacy alias.

Checked one at a time on the two-cell tutorial init. Five match to 0.000e+00
over 904 values under `RK5Adaptive`. `OneWall` is the exception and it is not a
porting error: it is the only one of the six acting on wall variables, so it
runs into the legacy RK5 defect. Under `Euler` it matches to 0.000e+00 as well,
which means this batch also reproduces that documented bug independently.

### `legacy/creation.cc` - 7 classes, 2026-09-17

`Creation::SpatialCylinder`, `SpatialRing`, `SpatialCoordinate`,
`SpatialPlane`, `FromList`, `OneGeometric`, `Sinus`, each with its legacy
alias (`Sinus` registers as `creationSinus`, lower case, in legacy).

All seven match to 0.000e+00 over 904 values under `RK5Adaptive`. None of them
touch wall variables, so none hits the legacy RK5 defect.

One legacy oddity reproduced deliberately: `Creation::FromList` divides by
`T.cell(k).calculateVolume()` where `k` is the loop counter, not the listed
cell `variableIndex(1,k)` it adds to. That looks like a defect, but models
written against legacy depend on it and the two coincide whenever the list is
`0,1,2,...`. Flagged in the source; worth revisiting if anyone relies on a
sparse list with the number flag set.

`Creation::Sinus` keeps legacy's `6.28` rather than `2*pi` - the period that
implies is what existing models were tuned against.

### `legacy/massAction.cc` - 6 classes, 2026-09-17

`MassAction::OneToTwo`, `OneToTwoWall`, `TwoToOneWall`, `GeneralWall`,
`HillSimple`, `GeneralEnzymatic`. Legacy registers no no-colon aliases for
these.

All six match to 0.000e+00 over 904 values - the three cell reactions under
`RK5Adaptive`, the three wall reactions under `Euler` for the reason above.

The wall variants treat wall variables as *paired*: each listed index `k` also
drives `k+1` with its own rate, which is how legacy represents the two sides of
a wall. `GeneralWall` and `GeneralEnzymatic` keep legacy's `rate > 0` gate, so
a zero or negative rate leaves the compartment untouched rather than running
the reaction backwards.

### `legacy/growth.cc` - 5 classes, 2026-09-17

`WallGrowth::StressSpatial`, `StressSpatialSingle`, `StressConcentrationHill`,
`ConstantStressEpidermalAsymmetric`, `Force`, with their no-colon aliases
(`WallGrowth::Force` has none in legacy).

All five match to 0.000e+00 over 904 values, compared under `Euler` since they
write wall lengths.

The stress term - summed wall variables, or the extension (d-L)/L when the
strain flag is set - was identical in four of the five and is now one helper.
`StressConcentrationHill` lets both neighbouring cells contribute to its Hill
factor, so the factor runs to 2 on an interior wall and 1 at the boundary;
that asymmetry is legacy behaviour and is preserved.

Still outstanding in this file: the six `CenterTriangulation` variants and
`Hypocotyl3D::StrainTRBS`, which need the TRBS strain machinery.

### Six small files, finished off - 12 classes, 2026-09-17

`legacy/pressure2D.cc` (3), `legacy/bending.cc` (4), `legacy/sisterVertex.cc`
(2), `legacy/turgorGrowth.cc` (1), `legacy/cellTime.cc` (1),
`legacy/dilution.cc` (1). Each of these six legacy files is now empty of
outstanding classes.

Seven of the twelve match legacy to 0.000e+00: `CellTimeDerivative`,
`Dilution::FromVertexDerivs`, `TurgorGrowth::WaterVolume`
(`WaterVolumeFromTurgor`, both flag settings),
`Pressure2D::{AreaPotentialTri,AreaPotentialTriSpatialThreshold,AreaPotentialTargetArea}`
and `Bending::{AngleInitiate,AngleRelax}` - the last two over 123376 values on
the 267-cell `meristem.init` as well as the two-cell fixture.

Two of these needed a model with more than one interior cell to be a real
test: `AreaPotentialTri`'s internal-cells-only flag and
`AreaPotentialTriSpatialThreshold`'s apex gate are both no-ops on a two-cell
mesh where every cell touches the background. Both were rerun on
`meristem.init` (267 cells, 573 vertices) and match exactly there.

Some reactions need a non-trivial setup or the comparison passes vacuously.
`Bending::Angle` on a square mesh has theta == theta_0 at every corner and so
produces no force at all - it "matched" legacy before it was tested against a
rest angle that differs from the current one. `Pressure2D::AreaPotentialTargetArea`
likewise needs the target area to be driven away from the actual area, so it
was paired with `Creation::Zero` (inflating), `Degradation::One` (contracting)
and the no-contraction flag in turn.

Three diverge from legacy on purpose - README items 8, 9 and 10, each a
wrong-index bug in the legacy source:

- `Bending::NeighborCenter` wrote force to `vertexDerivs[k]`, the cell-local
  wall counter used as a global vertex index.
- `Bending::Angle` added the central-vertex term to `jm` twice instead of once
  to `j`, so the three per-turn contributions did not sum to zero.
- `SisterVertex::SpringCellConc` tested `cell1` twice in its activity gate.

The first two are validated instead against an independent Python
reimplementation of the force law, `tests/port/bending_refcheck.py`: it applies
one explicit Euler step and checks the simulator lands in the same place
(0.000e+00 for both), and separately asserts that `Bending::Angle`'s three
contributions sum to zero, which is the property the legacy bug destroys.
`SisterVertex::SpringCellConc` is compared against legacy directly and matches
exactly on every input the gate bug does not reach (including pair breaking
under `BreakLength`); reversing the order of a sister pair whose two cells
differ is what exposes it.

Two legacy warts in `bending.cc` are kept rather than fixed, and documented in
the source: pi written as `3.14159`, and `Bending::Angle`'s clamp of
`cos(theta)` to +/-0.999, which leaves a straight chain with a residual
0.045 rad of bending strain. Both are load-bearing - the rest angle written by
`AngleInitiate` is compared against the angle computed by `Angle`, so the two
must use the same constant, and init files carry rest angles saved from legacy
runs.

`SisterVertex::InitiateFromFile` reads a file named literally `sister` from the
working directory. That is hard-coded in legacy and kept, so existing model
directories still work.

### `legacy/transport.cc` - 9 classes, 2026-09-17

`Diffusion::{MembraneSimple,SimpleOne,ConductiveSimple,2D}`,
`ActiveTransportCellEfflux`, `ActiveTransportCellEffluxMM`,
`DiffusionActiveTransportCell`, `ActiveTransportWall`,
`InfluxActiveTransportCell`, with their aliases. `legacy/transport.cc` is now
empty of outstanding classes.

Every init shipped with the tutorials leaves **wall variables at zero**, and
these reactions are driven by paired wall variables (PIN, AUX/LAX, wall
auxin). Comparing on a stock init therefore transports nothing and passes
vacuously. `tests/port/seed_init.py` rewrites an init with deterministic
non-zero values in both tables; the fixtures used here are seeded copies of
`twoSquare.init` and the 267-cell `meristem.init`.

All nine match legacy to 0.000e+00 - the cell-only reactions under
`RK5Adaptive`, the wall-writing ones under `Euler`.

Two of them take an optional third index level that records the flux in a
paired wall variable, and those two forms **do** differ from legacy. The
difference is entirely in that diagnostic field; the integrated trajectory is
bit-identical. It is README item 5 (redundant duplicate derivative evaluations
removed), and it shows up here because these two writes go to `wallData`
during `derivs` and so are not side-effect free:

- `InfluxActiveTransportCell` *accumulates* (`+=`) each evaluation. Legacy's
  total is larger than v2's by exactly (2n+1)/(n+1) for n Euler steps -
  measured 1.5010, 1.6680, 1.8016, 1.8907 at n = 1, 2, 4, 8 against predicted
  1.5, 1.6667, 1.8, 1.8889. Legacy evaluates derivs twice per step and once
  more before the first print; v2 evaluates once per step. The per-evaluation
  contribution is identical, which is what that exact ratio shows.
- `DiffusionActiveTransportCell` *assigns*, so only the last evaluation counts
  and the two agree exactly at every print point after the first. They differ
  at t=0 only, where legacy's extra pre-print evaluation has already
  overwritten the field and v2 still shows the initial state.

Neither field is read by the reaction or fed back on. Read it as a diagnostic
whose scale depends on the solver, not as a rate.

Two further legacy behaviours, deliberately treated differently from each
other:

`Diffusion::2D` reads its two cell centres from `Cell::positionFromVertex()`,
the no-argument overload, which returns the *cached* vertex positions rather
than the live `vertexData` in the same expression as the wall length and cell
area. The cache is only refreshed by `BaseSolver::setTissueVariables()`, i.e.
at print points. v2 reads all three from `vertexData` (README item 11). This
is not a guess: on a static mesh the two agree to
0.000e+00, and adding a pressure term that moves the vertices makes them
diverge by 2.5%.

`Diffusion::2D` also visits every interior wall once from each side with no
ordering guard, so each flux is applied twice and the effective diffusion
constant is 2*p_0. That one is **kept**, since halving it would silently
rescale every model written against it - the doubling is reproduced
explicitly, and the bit-exact match on a static mesh confirms it.

### `legacy/membraneCycling.cc` + `membraneCyclingAll.cc` - 14 classes, 2026-09-17

All eleven registered `MembraneCycling::` reactions and all three
`MembraneCyclingAll::` ones. Both files are now empty of outstanding classes.
(`MembraneCycling::CrossMembraneLinear` exists in the source but the legacy
factory has no entry for it, so no model can reach it; not ported.)

These share one shape - walk each cell's walls, move carrier on at k_on and
off at k_off across the membrane facing that cell, modulated by something -
and the port factors that into one `cycle()` helper taking the net on-rate per
membrane. The `MembraneCyclingAll::` variants are the same rules with the "the
wall must separate two real cells" test dropped.

Ten match legacy to 0.000e+00, compared under `Euler` since they all write
wall variables. They were compared on a seeded init
(`tests/port/seed_init.py`) for the reason given in the transport entry: the
carrier lives in paired wall variables, which every shipped init leaves at
zero.

The other four diverge on purpose - README item 12. Each has a typo in one of
its two mirror-image branches, so the result depends on which of a wall's
cells the init file happened to list first. That is not a property any of
these rules can sensibly have, which gives a test that does not need legacy at
all: `tests/port/membrane_swapcheck.py` rewrites the init with every wall's
cell1/cell2 flipped and its paired wall variables swapped to match, and
compares the cell variables.

    reaction                                        legacy      v2
    LocalWallFeedbackLinear                         1.001e-01   0.000e+00
    PINFeedbackLinear                               7.890e-02   0.000e+00
    PINFeedbackNonLinear                            4.572e-03   0.000e+00
    All::LocalWallFeedbackNonLinearInhibition       1.644e-02   0.000e+00

Run against the ten reactions that needed no fix, legacy is invariant to
0.000e+00 too, so the check is discriminating rather than trivially satisfied.

Two legacy behaviours kept: `CellFluxExocytosis` removes carrier from a
membrane in proportion to the square of the outward auxin flux it carries but
does nothing at all to an inward-facing one - an asymmetry, not an oversight
in transcription. And `LocalWallFeedbackNonLinearInhibition`'s two terms are
inhibited by different functions, `1/(x^n + k^n)` against
`1/(1 + x^n/k^n)`, which differ by a factor `k^n`; that rescales k_on against
k_off but leaves both monotonically inhibited, so it is preserved.

`TISSUE_REGISTER_REACTION` pastes the type name into an identifier, so
template instantiations are registered through a `using` alias. Several of
these reactions differ only in a boundary flag or a linear/non-linear switch,
and are one template each.

### TRBS: shared kernel, then `VertexFromTRBS` + ConcentrationHill, 2026-09-17

`legacy/mechanicalTRBS.cc` is 14.5k lines, and most of that is the same
~200-line triangle routine pasted seven times with a handful of lines changed
each. Porting the remaining variants one at a time would have imported that
duplication wholesale, so the kernel came out first:
`src/tissue/reactions/trbs_core.h` (318 lines) holds the element geometry,
the stiffness from the Lame pair, the node forces, the per-triangle Cauchy
stress and its rotation to the global frame, the Jacobi diagonalization, and
the principal-pair extraction.

The already-ported `VertexFromTRBScenterTriangulation` was rewritten onto it
first and re-checked: still 0.000e+00 against legacy over 946 values, so the
extraction is behaviour preserving to the bit. Its per-cell loop then came out
as `ctTrbsEvaluate()`, taking Young's modulus as a per-cell callback, which is
the only thing the center-triangulated variants actually vary.

Two new reactions on that base, both 0.000e+00 against legacy:

- `VertexFromTRBScenterTriangulationConcentrationHill` - the same material
  with `Y = Y_min + Y_max K^n/(K^n + c^n)`, i.e. a species that softens the
  wall as it accumulates. Legacy ships it as a 250-line copy of the base
  reaction with that one line different; here it is the callback and 40 lines
  of interface. Matches under both `Euler` and `RK5Adaptive`.
- `VertexFromTRBS` - the same elasticity with no center triangulation: the
  cell is one triangle and its three walls are the edges. Needed a triangular
  3D fixture (`tests/port/tri3D.init`), since every mesh shipped with the
  tutorials is quadrilateral. Exact under `Euler`; 3e-16 under `RK5Adaptive`,
  which is the last-bit cost of this port computing cotangents algebraically
  instead of through legacy's `acos`/`tan` round trip, and of the sorted Heron
  formula for the resting area.

Both accept an optional third index level to store the cell stress state,
matching the extension already on the base reaction; legacy has no equivalent.

One behaviour change while the principal-pair code was being extracted, and it
is not cosmetic. The rule picked the two *largest* eigenvalues as the in-plane
pair, on the argument that a membrane carries no load through its normal so
the normal eigenvalue is ~0. That holds in tension. In compression all three
in-plane values are negative, the ~0 normal sorts to the top, and the reaction
reported the shell normal as the principal stress direction with zero
anisotropy and zero sigma1 - wrong, and silently so. Ranking by |eigenvalue|
instead is identical in tension and correct in compression. Measured on a
center-triangulated sheet:

    tension      before dir=(0,1,0) a=0.675896 s1=0.189098
                 after  dir=(0,1,0) a=0.675896 s1=0.189098   (unchanged)
    compression  before dir=(0,0,1) a=0        s1=0
                 after  dir=(0,1,0) a=0.653949 s1=-0.210907

The anisotropy is now `1 - |s2|/|s1|` so it stays in [0,1] when both principal
values are negative, and sigma1 keeps its sign. Nothing calibrated in tension
moves.

Still outstanding in `mechanicalTRBS.cc`: the five anisotropic (MT) variants,
about 5700 live lines between them. They add a fibre direction and a
transversely isotropic material, plus a `parameter(4)` switch over seven
anisotropy modes and some ad-hoc branches keyed on hard-coded cell variable
indices. The kernel above already takes the coefficient pair the transversely
isotropic form needs (`lambdaT + 2 mioT`, `2 mioT`), so the material part is
in place.

### `legacy/growth.cc` center-triangulation stress rules - 2 classes, 2026-09-17

`WallGrowth::CenterTriangulation::Stress` and `StressConcentrationHill`, with
their aliases. Both are Lockhart yielding applied to the internal
centre-to-vertex edges rather than to walls: an edge grows once its stretch
`(d - L)/L` passes a threshold, the Hill variant raising the rate by
`k_hill c^n/(K^n + c^n)`. They share one `ctStretchYield()` helper. Both
0.000e+00 against legacy (2.2e-16 for the Hill variant under `RK5Adaptive`).

Both legacy versions carry a `stress_flag` selecting a stored wall stress
instead of the stretch, and neither implements it: `Stress` exits at
construction, and `StressConcentrationHill` prints an error per vertex on
every evaluation and then grows nothing at all. Both reject the flag at
construction here, so a model setting it is told once instead of silently
producing no growth while filling stderr.

`tests/port/trbs.sh` drives all of the above plus the TRBS elasticity, and
`tests/port/tri3D.init` is the triangular 3D fixture `VertexFromTRBS` needs.

### `legacy/mechanicalSpring.cc`: shared driver + 4 variants, 2026-09-17

Same treatment as the TRBS kernel, for the same reason: every reaction in
`mechanicalSpring.cc` is one loop with two things changed - which walls act,
and how stiff each one is. The coefficient idiom `K (1/L - 1/d)`, zeroed only
when both `d` and `L` are non-positive, scaled by `frac_adh` once stretched,
is identical in all nineteen. That loop is now `wallSpringLoop()`, taking an
`active(w)` gate and a `stiffness(w)` callback and keeping both execution
paths the existing `WallMechanics::Spring` had (serial in legacy's exact
accumulation order for small tissues, a two-pass gather for large ones).

`WallMechanics::Spring` was moved onto it first and re-checked against legacy
before anything was added. Four variants follow, all 0.000e+00 against legacy
on the two-cell fixture and on the 267-cell mesh:

- `WallMechanics::SpringEpidermal` - boundary walls only.
- `WallMechanics::SpringEpidermalCell` - every wall of a cell that touches the
  background, not just the boundary wall itself.
- `WallMechanics::SpringConcentrationHill` - stiffness from an inhibitory Hill
  function summed over the wall's two cells. A boundary wall gets one
  contribution rather than two, so its stiffness runs to `K_min + K_max`
  instead of `K_min + 2 K_max`; that asymmetry is legacy's and is preserved.
- `VertexFromWallBoundarySpring` - `Spring` restricted to boundary walls.
  Unlike `SpringEpidermal` it leaves the saved force on interior walls
  untouched instead of zeroing it, which is also legacy's behaviour.

Any form that saves the force into a wall variable differs from legacy **at
t=0 only**, and this is the third place the same cause has turned up: the save
is a write to `wallData` from inside `derivs`, so it is not side-effect free
and it sees legacy's extra pre-print derivative evaluation (README item 5).
Measured on `WallMechanics::Spring` with a save index, one Euler step:

    legacy  t=0  0           t=0.04  -0.000789479
    v2      t=0  0.808275    t=0.04  -0.000789479

Legacy's 0 at t=0 is the force it computed before the first print, which is
genuinely zero there because the mesh starts at rest; v2's value is the
untouched initial contents of that wall variable. Every later print agrees
exactly. Worth noting that this was *pre-existing* - it was not introduced by
the refactor, and was confirmed by rebuilding the previous commit and getting
the identical 8.083e-01 discrepancy - it had simply never been exercised,
because the earlier check of `WallMechanics::Spring` never passed a save
index.

Still outstanding in this file: the microtubule (MT) spring variants, which
need a per-cell fibre direction and the angle between it and each wall, plus
`VertexFromExternalSpring*`, `cellcellRepulsion` and `vertexFromSubstrate`.

### `legacy/mechanicalSpring.cc`: the two MT springs, 2026-09-17

`VertexFromWallSpringMT` and `VertexFromWallSpringMTConcentrationHill`, both
0.000e+00 against legacy on the two-cell fixture and the 267-cell mesh.

These are the fibre-reinforcement rule: a wall running along a cell's
microtubule direction is soft, one running across it is stiff, through
`K = K_min + K_max (2 - cos^2 t1 - cos^2 t2)`. They read the direction from
*cell variables* - components at a given index followed by a flag that must
exceed 0.5 - not from the legacy direction machinery, so nothing unported is
needed. `mtCosSq()` is shared between them.

Two legacy quirks kept, both in the ConcentrationHill variant. Its fallback
for a cell with no direction is 0.5 rather than the 1.0 the plain MT variant
uses. And the concentration factor is computed *inside the same guard* as the
direction, so a cell whose direction flag is unset contributes a Hill factor
of 0 - its concentration is ignored entirely rather than softening the wall.
That couples two things that look independent, and I first transcribed it as
independent, which showed up immediately as a 7e-3 mismatch. Unlike the
cell1/cell2 asymmetries in membraneCycling.cc there is no invariance the rule
must satisfy here, so it is reproduced rather than "fixed". The harness tests
both a direction index where both fixture cells are flagged and one where only
one is, since that is the path the coupling changes.

### `legacy/mechanicalTRBS.cc`: the transversely isotropic material, 2026-09-17

`VertexFromTRBScenterTriangulationMT` - the fibre-reinforced material itself,
about 1800 live lines in legacy (573 of `derivs`, 1231 of `update`), and the
one the microtubule feedback in the hook work runs on.

It decomposes cleanly onto the kernel extracted earlier. The isotropic part of
the element force is the ordinary TRBS force evaluated with the *transverse*
moduli, which `stiffnessFrom(el, lambdaT + 2 mioT, 2 mioT)` already gave; the
fibre adds a correction from the extra Lame pair and the fibre direction
pulled back into the element's rest frame with the cofactor of F. That
correction is now `anisotropicDeltaForce()` in `trbs_core.h`, and the element
frame it needs - shape vectors, deformation gradient, local-to-global
rotation - came out as `localFrameOf()`, which the isotropic Cauchy stress was
already rebuilding inline. `addCauchyStress` was moved onto it and re-checked
at 0.000e+00 first.

`derivs` matched legacy on the first run, and all eleven MF (material) flags
match to 0.000e+00, under both plane stress and plane strain, with and without
wall dynamics. `update` - the pass that computes each cell's area-averaged
true strain and true stress, diagonalizes both, and stores the anisotropy
measures, energies, Mises stress and cell normal that other reactions read -
matched on its first run too, then needed three corrections found by the
harness:

- Legacy's `update` has a *different* material chain from its `derivs`: no
  branch at all for MF flag 9, and no ad-hoc `cellData[40] == 100` switch in
  MF flag 0. Neither difference looks intentional; both change results, so
  both are reproduced (`moduliFor(..., inUpdate)`).
- Legacy's two stress passes differ: the direct one ranks the principal values
  by magnitude and records the Mises stress, the neighbour-averaged one ranks
  by signed value and records none.
- The two Jacobi loops also differ - the stress pass finds its pivot at the
  top of the loop and so performs one extra rotation below threshold, the
  strain pass finds it before and refreshes at the bottom. Both forms are in
  `trbs_core.h`. Neither of legacy's loops is bounded; both are capped at 50
  sweeps here, since a degenerate tensor otherwise hangs the run.

One deliberate divergence, README item 13: legacy's neighbour weighting is
dead code, guarded by `size_t > -1`. See the README for the measurements. The
diagnosis was confirmed both ways - emulating the broken guard reproduces
legacy bit-for-bit, and instrumenting the fixed loop shows it reaching one
neighbour per cell on the two-cell fixture and contributing exactly
`w * neighbour`.

Three further legacy behaviours, documented in the source rather than changed:

- The **MT update flag** (parameter 9) is documented in the constructor as
  offering four feedback modes - force to stress, to strain, to
  perpendicular-strain - and the constructor validates it in 0-4. Only mode 1
  (set the direction from the TETA angle) does anything. Nothing in this
  reaction rotates the direction toward stress or strain; the stored direction
  is only ever normalized.
- **MF flag 10's Hill constants** are read from parameters 9 and 10 - the MT
  update flag and the double-resting-length flag - although the 13-parameter
  form exists precisely to carry them at 11 and 12, where nothing reads them.
  The mode is unusable as written: the Hill exponent can only be 0 or 1 and
  the half-max only 0-4, and neither flag can then be set independently. The
  11-parameter form keeps legacy's reading so existing models reproduce
  exactly; the 13-parameter form uses the documented positions.
- The **sliver branch**: a near-degenerate element makes legacy *double* the
  derivatives already accumulated on its three nodes and drop its own force.
  That is not a meaningful remedy, but it is reachable, so it is reproduced
  with its warning rather than reinterpreted.

Two smaller changes. The degenerate-fibre fallback draws from the seeded
generator the rest of the rewrite shares rather than libc `rand()`, so runs
stay reproducible; it is a fallback for a fibre lying along the element normal
and the validation runs do not reach it. And MF flags 0 and 6-9 read cell
variables 40 and 13 by hard-coded number, which legacy would run off the end
of the row without saying so - `initiate()` now checks the row is wide enough
and says which mode needs it.

The MF flags whose `derivs` writes a cell variable (1 writes the transverse
modulus, 5 reads the anisotropic energy, and MT update flag 1 writes the
direction) differ from legacy at the **t=0 print only**, the same
derivs-side-effect timing as the transport and spring reactions. The
trajectories are identical: measured on MF flag 1, legacy prints 3.60489 and
3.43234 for the two cells at every print, v2 prints the untouched initial
values at t=0 and then the same 3.60489 / 3.43234. The harness compares those
configurations with `--block=vertex`, which `tools/port/compare.py` now
supports for exactly this purpose.

Still outstanding in this file: `VertexFromTRBSMT` (the same material without
center triangulation), `VertexFromTRLScenterTriangulationMT`,
`VertexFromTRBScenterTriangulationConcentrationHillMT`,
`VertexFromTRBScenterTriangulationMTOpt`, and
`Hypocotyl3D::VertexFromTRBScenterTriangulationMT`. All five are variants of
the reaction ported here and now have the whole kernel available.

### `VertexFromTRBScenterTriangulationConcentrationHillMT`, 2026-09-17

The fibre-reinforced material with both moduli set by inhibitory Hill
functions of one cell concentration. Legacy ships it as its own ~1000-line
copy, but it is *not* the same material law as
`VertexFromTRBScenterTriangulationMT` with a different modulus - four things
differ, and each of the first three showed up as a shrinking mismatch as they
were found:

- the shear modulus convention is `Y/(1+p)` rather than `Y/(2(1+p))`, so the
  stiffness coefficients read `(lambdaT + mioT, mioT)`;
- the trace entering both stress terms is the cotangent form
  `(sum Delta_i cot_i)/(4 A_rest)`, not `tr(E)`;
- the fibre is pulled back to the rest frame **barycentrically** - take the
  point one unit along the fibre from the current centroid, express it in
  barycentric coordinates of the current triangle, map those onto the rest
  triangle, subtract the rest centroid - rather than by applying the cofactor
  of F. The two agree for a rigid motion and differ under shear, but the
  consequential part is that this one keeps the *length* of the result and
  scales the anisotropic Lame pair by it, so an element whose deformation
  shortens the fibre direction also weakens its anisotropy;
- the anisotropic force is not a stress tensor pushed through the shape
  vectors at all. It differentiates the invariants I1, I4 and I5 directly with
  respect to the node positions.

All four are now separate pieces in `trbs_core.h`
(`hillVariantDeltaS`, `barycentricFibre`, `anisotropyInvariantsFrom`,
`invariantAnisotropicForce`), so the two materials share the geometry and
differ only where they actually differ.

Exact against legacy (0.000e+00 over 1744 values) in the two-level form. The
six-level form, which stores the maximal and second strain and stress
directions, writes those from `derivs`, so it differs at the **t=0 print
only** - the same derivs-side-effect timing as everywhere else. Measured: the
stored strain direction is `0 0 1 0` in both simulators at every print after
the first; at t=0 legacy shows the value its extra pre-print evaluation wrote
and v2 shows the untouched initial contents.

Two legacy slips are reproduced rather than corrected, because both change the
answer and neither has an invariance that forces the issue: the rest and
current edge arrays are cyclically rotated just before the invariant-derivative
block while `cotan` and `Delta` are not, and the `derIprim1` accumulation
multiplies by `position[m]` inside a loop over `i`, so the i-sum only ever
scales one position vector.

Still outstanding in `mechanicalTRBS.cc`: `VertexFromTRBSMT` (the same
material without center triangulation), `VertexFromTRLScenterTriangulationMT`,
`VertexFromTRBScenterTriangulationMTOpt` - which is not a material variant at
all but a simulated-annealing optimiser over the anisotropy direction, with
its own temperature and annealing-rate parameters - and
`Hypocotyl3D::VertexFromTRBScenterTriangulationMT`.

### `VertexFromTRBSMT`, 2026-09-17

The same transversely isotropic material with no center triangulation: the
cell is one triangle and its three walls are the edges, so everything the
center-triangulated version defers to an `update()` pass happens in `derivs()`
here - there is only ever one element per cell to average over.

Forces are exact against legacy (0.000e+00), for both material flags this
variant supports.

The stored diagnostics - strain and stress anisotropy, area ratio, the two
energies, the stress projected on the fibre, the stress tensor - are written
from `derivs`, and legacy samples them one derivative evaluation later than v2
does (README item 5 again). That is not a formula difference and it is
measurable as such: the gap is exactly linear in the step size, halving with
it, on the area ratio at the first print -

    h          legacy      v2          difference
    0.001      0.97536     0.975396    3.600e-05
    0.0005     0.975367    0.975386    1.900e-05
    0.00025    0.975371    0.97538     9.000e-06

- and the vertex trajectory is bit-exact throughout.

Three more legacy behaviours reproduced rather than corrected:

- The fibre is pulled back with `F^T a` here, where
  `VertexFromTRBScenterTriangulationMT` uses the cofactor of F and
  `...ConcentrationHillMT` uses a barycentric map. Three variants of the same
  material, three different pullbacks. Legacy also normalizes the result
  twice, which makes this variant's degenerate-direction fallback
  unreachable, so there is none in the port.
- MT update flag 1 writes another hard-coded four-cell setup, and another
  parameter collision: cells 0 and 2 take their angle from parameter 5, which
  is also the neighbour weight, and cells 1 and 3 from parameter 8. Cells
  beyond the fourth are left alone entirely.
- After the per-cell loop, **cell 0's iso and aniso energy slots are
  overwritten with the tissue totals**, destroying that cell's own values.
- Legacy stores the *maximal* strain value alongside the perpendicular
  direction, not the second one, where the center-triangulated variant stores
  the second.

Its neighbour-weighted stress pass carries the same dead
`size_t > -1` guard as the center-triangulated one (README item 13) and is
fixed the same way.

While porting it, legacy's four Jacobi loops became five, so they are now one
`jacobiEigen3(A, eig, style)` with five named styles in `trbs_core.h`. No two
of legacy's loops are the same - they differ in the rotation formula
(`0.5 atan` vs half-angle), in the angle used for a degenerate pair (pi/4 with
pi spelled 3.1415 in one place and 3.14159265 in another), in whether the
pivot is taken at the top of the loop or refreshed at the bottom, and in
whether the eigenvector columns are renormalized. The consolidation was
checked by re-running every TRBS harness at 0.000e+00 first.

### `VertexFromTRLScenterTriangulationMT`, 2026-09-17

Triangular *linear* elements on a center-triangulated cell: the same fibre
direction, the same stored state, but the mechanics are an orthotropic
St Venant-Kirchhoff triangle rather than biquadratic springs. There is no
spring term at all - the whole force is the stress pushed through the shape
vectors, so it reuses `pushDeltaS` and nothing else of the TRBS force path.

The stress is taken in the *fibre* frame: rotate the Green strain into it,
apply the orthotropic stiffness in Voigt form, rotate back. Two things about
that stiffness are worth noticing, and both are legacy's: the Poisson ratio in
the whole in-plane block is the *longitudinal* one, and the shear modulus is
built from the longitudinal Young's modulus alone - the transverse ratio
(parameter 3) only ever reaches the Lame constants, which this reaction
computes and then never uses.

Forces exact against legacy (0.000e+00) for every material flag it supports
(0, 1, 2, 5), under both plane assumptions and with neighbour weighting on.
Its stored diagnostics carry the same one-evaluation sampling offset as
`VertexFromTRBSMT`, linear in the step size (4.1e-5, 2.1e-5, 1.0e-5 at
h = 0.001, 0.0005, 0.00025), with the vertex trajectory bit-exact throughout.

Flag 2 treats `Y_matrix` as a *total* stiffness and sets the transverse
modulus to `Y_matrix - Y_fibre`, so it needs `Y_matrix > Y_fibre`; below that
the modulus goes negative and the `sqrt(Y_L Y_T)` in the stiffness matrix
returns NaN. That is a parameter-domain limit, not a port issue, but nothing
in legacy checks it.

### `Hypocotyl3D::VertexFromTRBScenterTriangulationMT`, 2026-09-17

Legacy keeps this in `hypocotyl3D.cc` as a full 2500-line copy of
`VertexFromTRBScenterTriangulationMT`. A token-level diff of both the derivs
and the update says it differs in exactly four places, so it is the same class
here with a `Hypocotyl` template parameter:

- it adds **MF flag -1**, a tissue-layer material keyed on cell variable 37:
  a fixed base modulus of 50 with the Poisson ratios pinned at 0.2, scaled per
  layer by the first four parameters - -1 epidermis (both directions), -2
  inner (hoop only, axial left alone), -3 anticlinal axial, -4 anticlinal
  transverse;
- it drops MF flags 9 and 10, and the ad-hoc cell-variable-40 switch inside
  flag 0;
- its derivs writes no cell variable at all, not even the transverse modulus
  that flag 1 reports in the base;
- under MF -1, parameters 2 and 3 are layer scalings rather than Poisson
  ratios, so legacy *skips its Poisson range check* for that flag. Reproduced:
  without it a perfectly good hypocotyl model is rejected, since those
  scalings are routinely above 0.5.

Everything else, the update pass and the neighbour averaging included, is the
same code.

Exact against legacy (0.000e+00 over 1744 values) for MF flags 0, 1, 2 and -1,
the last on a fixture carrying real layer codes rather than seeded noise -
without that, no layer branch fires and both cells get the same modulus, which
tests nothing. The neighbour-weighted cases diverge by README item 13, the
same dead `size_t > -1` guard as everywhere else.

A caution about reading this file: my first diff of the update missed its last
2400 tokens, because I bounded the function with a brace count that terminated
early on a comment. That made it look as though this variant had no
neighbour-averaging pass at all, and the port briefly omitted one. The
harness caught it - a neighbour-weighted run left the stress anisotropy at its
initial value while legacy wrote one.

One incidental find while building the harness for this variant: legacy's
unbounded Jacobi loops are reachable from ordinary parameters, not just
pathological ones. MF flag 2 treats `Y_matrix` as a total stiffness and sets
the transverse modulus to `Y_matrix - Y_fibre`, so any model with
`Y_fibre > Y_matrix` has a negative modulus, and the legacy binary then runs
indefinitely - measured still going after 25 s on the two-cell fixture, where
the whole run otherwise takes under a second. This port finishes immediately
because every Jacobi is capped at 50 sweeps (README item 14). The same trap
exists in `VertexFromTRLScenterTriangulationMT`, where a negative modulus also
makes `sqrt(Y_L Y_T)` return NaN. Nothing in legacy checks the sign.

### `VertexFromTRBScenterTriangulationMTOpt`, 2026-09-17

Named for a simulated-annealing search over the anisotropy directions that was
never finished. What is left is an **energy probe that applies no force**:
`derivs` writes nothing to any derivative array, `update` is a single
`std::cout` of four energies, the accept/reject step is commented out, and the
three annealing parameters (initial temperature and the two annealing rates)
are never read.

Its Monte Carlo proposal is inert too, which is the only reason it can be
compared against legacy at all. Vertex positions and cell centres are
perturbed by `posStep * (rand() - 0.5)` with `posStep` hard-coded to **zero**,
and the anisotropy direction *is* randomized into the shadow state but the
energy then reads the unperturbed direction straight out of `cellData`. So the
energy is always the energy of the current state and the reaction is
deterministic, despite re-seeding libc `rand()` from the wall clock on every
derivative evaluation. Two legacy runs confirm it: identical output, with and
without strain present.

That being so, the port implements the effective behaviour and describes the
dead scaffolding rather than transcribing it. The evidence that the two are
the same thing is the comparison: 0.000e+00 over 828 values for both material
flags, energies included, since those go to stdout and `compare.py` reads
them. And the property that matters for anyone using it - that adding it to a
model changes nothing - is checked directly: the tissue state is identical
with and without the reaction.

Two further legacy details reproduced: the pressure term is a tetrahedron
volume measured against a hard-coded apex at (0, 0, -70), and it is *assigned*
rather than accumulated inside the element loop, so only the last element of
each cell contributes to it.

`legacy/mechanicalTRBS.cc` is now empty of outstanding classes: 14.5k lines,
twelve reactions, all ported and all checked against legacy.

### `legacy/network.cc`: the AuxinModelSimple1 family, 2026-09-18

`network.cc` holds 41 auxin/PIN network models in 5700 lines - families of
near-identical reactions differing in one production term, one Hill factor, or
whether a flux is weighted by geometry. Nearly all share one shape: a per-cell
block of production and degradation ODEs, then a polarised transport block
that walks the cell's interior walls, builds a normalising sum over the
neighbours, and exports auxin in proportion to the PIN allocated to each
membrane. `forEachInteriorWall`, `neighbourSum` and `onBoundary` carry that
traversal; each reaction supplies its own chemistry.

First batch, all 0.000e+00 against legacy on the two-cell fixture and the
267-cell mesh:

- `AuxinModelSimple1` and `AuxinModel1` - the same model, the second weighting
  each flux by wall length over cell volume. One template, one flag.
- `AuxinModelSimple1Wall` - polarises on a *wall* variable instead of the
  neighbouring cell's. Legacy shifts the whole set up by the most negative
  value so a negative signal cannot invert the denominator, but applies that
  correction inside the accumulation loop, where it only sees the walls
  visited so far. Reproduced as written.
- `AuxinTransportCellCellNoGeometry` - transport only, with the polarisation
  signal a Hill function of the neighbour's X. Boundary walls contribute the
  baseline k1 to the normalising sum even though nothing moves across them.

One behaviour worth naming, since three of these share it: when the
normalising sum comes out at exactly zero, `AuxinModelSimple1` sets the
polarisation rate to **1**, not to 1/n. A cell with no signal anywhere around
it therefore exports at the full rate through every membrane rather than
spreading one unit between them.

The optional membrane-PIN store writes to `wallData` from `derivs`, so it
differs from legacy at the **t=0 print only** - the same derivs-side-effect
timing as the transport, spring and MT reactions. Measured on the interior
wall of the two-cell fixture: legacy 0.230822 / 0.280055 at t=0 where v2 still
shows the seeded values, and both 0.238324 / 0.296522 at every print after.

### `legacy/network.cc`: the SimpleROPModel family - 7 classes, 2026-09-18

Seven variants of one ROP/PIN polarisation model, all 0.000e+00 against
legacy. PIN cycles between the cell's cytoplasmic pool and each membrane, and
what makes a membrane accumulate PIN is a Hill function of the PIN on the
*opposite* face of the same wall - the positive feedback that polarises
neighbouring cells against each other.

They divide in two. Models 1-3 carry auxin in a wall compartment, secreting it
into the wall and diffusing it across to the opposite face; models 4-7 move
auxin straight from cell to cell and make PIN production saturate in the
cell's own auxin. Level 1 still takes two indices in 4-7, but the first (the
wall auxin) is no longer read.

Within each half the differences are one term each, so the per-membrane
chemistry is written out in every variant rather than hidden behind flags:

- 2 differs from 1 only in that PIN removal is no longer proportional to the
  wall auxin on that face.
- 3 differs from 2 only in writing the constant return term last rather than
  first - algebraically the same reaction, and legacy registers both names.
- 5 adds a self-limiting removal: PIN comes off a membrane faster the more it
  already carries.
- 6 also cycles PIN on *boundary* membranes, where there is no neighbour to
  transport to and no opposite face to feed back from; legacy reads the cell1
  face there whichever side the cell is on, which is reproduced.
- 7 gates everything on a wall marker: walls flagged 1 run the polarising
  model, walls flagged 0 leak auxin and PIN at fixed rates instead. Its leak
  branch has a slip that is harmless and reproduced - on the cell1 side it
  sends the leak to `cell1`, which is the cell itself, so the two
  contributions cancel and only the cell2 side of a flagged-0 wall actually
  moves anything.

Testing 7 needed a fixture carrying real 0/1 markers: with the seeded values
neither branch fires and the comparison passes without exercising the
reaction at all. Both branches are covered now, including an all-zero variant
that runs only the leak path.

### `legacy/network.cc`: the gradient models - 5 classes, 2026-09-18

Five reactions on the `SimpleROPModel7` skeleton - marked walls run the model,
unmarked ones leak - but with the membrane PIN feedback coming from an auxin
concentration rather than from the opposite face:

    dPIN_membrane/dt = p8 PIN_membrane - p9 PIN_cell N^n/(K^n + D^n)

Where N and D come from is the entire difference between the first three:
the cell's own auxin in both (up the internal gradient), K in the numerator so
the term is repressed instead (down the internal gradient), or the
*neighbour's* auxin in both (up the external gradient). One template, three
aliases. All 0.000e+00 against legacy.

The other two are `DownInternalGradientModel` with a change each, and both
carry legacy branching that is worth naming because it is not what the names
suggest:

- `...SingleCell` drops the auxin transport, so PIN cycles but nothing moves
  between cells. Its interior-wall test is absent from both cycling branches,
  so PIN cycles on *boundary* membranes too; and its leak is attached as an
  `else` to the cell2 branch alone, so an unmarked wall leaks when the cell is
  its cell1 and does nothing when the cell is its cell2.
- `...Geometric` scales the fluxes by wall length over cell volume, g_ij for
  the donor and g_ji for the receiver. The PIN taken *off* a membrane is not
  scaled while the PIN added to the cell is, so PIN is not conserved between
  the two. Its wall length comes from `Wall::length()`, legacy's cached copy
  of wall variable 0 refreshed only at print points; this reads that variable
  live, so the two agree exactly unless a model grows wall lengths - the same
  divergence as README item 11.

### `legacy/network.cc`: the AuxinROPModel family - 3 classes, 2026-09-18

`SimpleROPModel` with the ROP made explicit: a third species cycles between
the cell and each membrane with the opposite-face Hill feedback, and it is the
*ROP* on a membrane that pulls PIN out of the cell there. All three 0.000e+00
against legacy.

`AuxinROPModel2` adds two saturating terms - auxin secretion saturating in the
cell's auxin, and PIN removal a Hill of the membrane ROP rather than linear in
it - which shifts the parameter numbering by three.

**`AuxinROPModel3` is not symmetric in a wall's two faces, and this is worth
Ross's attention.** Its cell1 branch uses a new secretion law with an
AUX/LAX-dependent influx, `(p3 + p4 PIN) a - (p2 + p16 AUX) a_wall`, while its
cell2 branch still carries the plain `AuxinROPModel` expression with no AUX
term at all. So `parameter(16)` and the AUX cell variable are read only on the
cell1 side, and the model's behaviour depends on which cell an init file
happens to list first for each wall - an arbitrary artefact of the file.
Measured with `tests/port/membrane_swapcheck.py`: flipping every wall's
orientation moves the result by **2.524e-01**, where `AuxinROPModel` is
invariant to 0.000e+00.

Unlike the cell1/cell2 asymmetries fixed in `membraneCycling.cc`, this one is
not a typo that can be read off the code - the second branch is the older
expression, intact, so it looks like a half-finished edit and there is no way
to tell which form was intended. It is therefore reproduced as written and
flagged here rather than repaired. Both simulators show the same 25% swing,
which is the evidence that the port is faithful; whether the model should
behave that way is a question for whoever uses it.

### `Pressure3D::CenterTriangulation::Linear`, 2026-09-18

Ported out of order, ahead of the rest of `mechanical.cc`, because the
benchmark sweep showed it was the single reaction blocking **all 118 models**
of `publications/eng_et_al_2021` - every one of them was otherwise already
supported. It is also needed by `bozorg_etal_2016` and `3Dhypocotyl`.

Pressure on a center-triangulated shell that ramps in over a set time rather
than being applied at full strength from the start, which is how those models
inflate a tissue without kicking it. Per CT triangle the force is
`k_force * A * n_hat`, given to all three nodes - each gets the whole thing,
not a third of it. `areaFlag` selects `A`: 1/3 (no area weighting at all),
`Area/2`, or `Area/2` with only the z component applied through a second ramp.
A fourth parameter ramps from `p3` to `p0` instead of from zero and reports the
current pressure into a cell variable.

Exact against legacy for all three area flags. The four-parameter form differs
at the **t=0 print only**, since it writes that cell variable from `derivs`:
measured, legacy reports 0.1 (the ramp's starting pressure) at t=0 where this
build still shows the init value, and both report 0.12 at every print after.

One legacy quirk kept: under `areaFlag 2` the *first* ramp is never advanced -
`update` only steps `timeFactor1` for flags 0 and 1 - so the main force term
stays at zero for the whole run and only the z-only term does anything.

### `legacy/fiberModel.cc`, 2026-09-18

Ported next because the benchmark sweep put `FiberModel` at the top of the
blocking list: 26 of the 32 models that legacy runs and this build cannot are
blocked by it, all in `publications/bozorg_etal_2014`.

Only two of the file's six classes are ported, because only two are reachable:
`baseReaction.cc` registers `General` (also under the bare name `FiberModel`,
which is the spelling every model in the corpus uses) and `Deposition`, while
the entries for `Linear`, `LinearEquilibrium`, `Hill` and `HillEquilibrium`
are commented out upstream. Those four have never been constructible from a
model file, so there is nothing to validate them against and nothing using
them; porting them would be writing new code, not porting.

Both act entirely in `update()`, so a model containing only a FiberModel rule
moves no vertex and the comparison is of the cell-variable table alone.
Exact (0.000e+00) against legacy over 17 configurations: both gradual branches
with the velocity gate admitting one cell and then both, the direct branch,
all three initiation flags, Deposition with its gate open and closed, its
`init_flag 3` Hill branch, a cell below its stress floor, and both no-ops.

Two traps in building that fixture, both of which had made a case vacuous:

- The gradual branches also stop at `Y_M + Y_F`, and the seeded cell sat
  exactly at that ceiling, so it was excluded whatever the velocity gate said
  and the two thresholds gave identical output. Dropping the cell's `Y_L`
  below the ceiling makes the gate the only thing separating the cells.
- `Deposition` compares a cell's maximal stress against a window hard-coded in
  legacy as `[8, 14]`, and `seed_init.py` spreads its values over `(0.05,
  0.95)` - so every cell fell below the floor, every share was zero and the
  redistribution never ran. `tests/port/cellvar_init.py` (new) writes an init
  with cell variables chosen per cell, which is what puts one cell inside the
  window and one above or below it.

Legacy behaviour reproduced rather than repaired, in each case because the
published results came out of the code as written:

- `General`'s direct branch (`linear-hill_flag=2`) omits the matrix term `Y_M`
  that every other branch and the documented equation include, so it writes a
  fiber contribution where the other branches write a modulus. It also takes a
  second index at level 1 (`FiberL_index`, which the class documentation says
  receives the update) and never writes it - the result goes to
  `Young_L_index` like the rest.
- `Deposition` multiplies by `k_rate` but not by the step `h`, so its approach
  to the target depends on how often `update()` is called rather than on
  simulated time, and its stress window is not exposed as a parameter.
- `Deposition`'s `init_flag` 2 and 3 seed the global `rand()` from the clock,
  so those two flags were never reproducible run to run, in legacy either.
  Kept as-is rather than given a seeded generator; the comparison above uses
  flags 0, 1 and 3, and flag 3's randomness is in `initiate` only.

Next blocker for these models is `UpdateMTDirectionEquilibrium`.

### `legacy/directionReaction.cc`, 2026-09-18

All six classes, the whole file, because all six are registered and the file
is small. `UpdateMTDirectionEquilibrium` is what the benchmark wanted next
(the `bozorg_etal_2014` models reach it right after `FiberModel`); the other
five came along at no extra cost.

A direction occupies `dimension` consecutive cell variables from the given
index, so one index names two variables in 2D and three in 3D.

Exact (0.000e+00) against legacy over 19 configurations in 2D and 3D
(`tests/port/direction.sh`), covering both `Continous` forms, both `Update`
forms with each of their one-, two- and three-parameter gate combinations,
`ConcenHill`, `RotatingDirection`, and a no-op for each.

The fixture needed a second reaction to be meaningful at all. `UpdateMTDirection`
and `UpdateMTDirectionEquilibrium` both `initiate()` the direction *to* the
target, so a model containing only one of them has nothing left to follow and
is a no-op for the whole run - it "matched" legacy while doing nothing. Every
case for those two therefore runs `RotatingDirection` on the target index as
well, which keeps the target moving and doubles as that reaction's own test.
Checked, not assumed: opening the velocity gate, changing the stress
threshold, and going from two parameters to three each change v2's own output.

Legacy behaviour kept deliberately:

- `myMath::pi()` is the truncated literal `3.14159265`, and the angle folding
  in `ContinousMTDirection` is sensitive to it, so `kLegacyPi` is that literal
  rather than `M_PI` - as in `bending.cpp`, which needs a differently
  truncated one.
- `UpdateMTDirectionConcenHill`'s alignment test is hard-coded to three
  components while everything around it loops over `dimension`, so in a 2D
  model it also reads whichever cell variable follows each direction. Kept;
  changing it would change 2D results and no 2D model is known to use it.
- That same reaction writes the direction straight into `cellData` from
  `derivs()` with no `h`, so an adaptive solver applies it once per
  *derivative evaluation* rather than once per step. It is a direction update
  wearing a derivative's clothes; the gates in `UpdateMTDirectionEquilibrium`
  exist because that is the form you want when the target is noisy.
- `UpdateMTDirection` does not align the two axes before stepping, while
  `UpdateMTDirectionEquilibrium` does. A target pointing the other way
  therefore drags the direction through zero length in the first and takes the
  short way round in the second. This looks like an omission in the older
  class rather than a choice, but both are reachable from model files and the
  difference is visible in results, so both are reproduced as written.

Next blocker for the `bozorg_etal_2014` models is `CalculateAngleVectors`.

### `legacy/calculate.cc`, 2026-09-18

All five classes. This file was three of the benchmark's top blockers at once:
`CalculateAngleVectors` (45 models), `TemplateVolumeChange` (30, an alias for
`Calculate::TissueVolumeChange`), and `Calculate::VertexVelocity`, which is
what supplies the velocity that `UpdateMTDirectionEquilibrium` and
`FiberModel` gate on - the three batches only work together.

It also settles a question the benchmark raised. 18 of the models this build
could not run failed on `maxVelocity`, and legacy fails on them too: legacy
*removed* `Calculate::MaxVelocity`/`maxVelocity` and now exits telling the
author to use `Calculate::VertexVelocity`. Those models are stale for both
simulators, not a gap in this port. This build refuses the old names the same
way, and for the same reason - the quantity is a mean, never was a maximum.

Exact (0.000e+00, or 1e-15 float-ordering noise) against legacy over 12 of 14
configurations in 2D and 3D (`tests/port/calculate.sh`).

The two exceptions are both `TissueVolumeChange`, and the cause is measured
rather than assumed. Of the three quantities it stores, the two that are
functions of the state alone - the accumulated vertex displacement and the
summed vertex speed - match legacy exactly at every print. The third is stored
as a difference against its own previous value, so it reports the gap between
two consecutive *derivative evaluations*; legacy performs one more evaluation
before printing than this build does, so its gap spans a different pair.

Measured with Euler at four step sizes from 5e-4 to 4e-3, the offset is
exactly linear in h (rel/h constant at 2.28), which is a one-evaluation lag,
not an arithmetic difference. Under `RK5Adaptive` it is 1.6e-9 and shows only
at intermediate prints - the last stage there lands on the accepted end-of-step
state, so the two agree once the tissue settles. This is README item 5 seen
through one extra level of differencing, and the harness runs those two cases
at 1e-8 with the reason written next to them.

Other legacy behaviour kept:

- `Calculate::AngleVector` range-checks its axis parameter against 0, 1 and 2
  and then exits with status 0 and "The code should be modified for 3d" if it
  is given 1 or 2, so only axis 0 has ever run. This rejects 1 and 2 at
  construction instead, which names the offending rule instead of exiting
  silently part-way through a simulation.
- Its pi is the 7-digit `3.141592` while `AngleVectorXYplane`, ten lines away,
  uses `3.14159265`. Both are kept as they are, per class.
- `AngleVectorXYplane` is documented as storing `abs(cos(angle))` and stores
  `atan(z/|xy|)`, an angle in radians. The code is what models were written
  against.
- `AngleVectors` and `AngleVector` normalise their *input* vectors in place,
  so they rewrite the directions they are asked to measure.
- Neither `TissueVolumeChange` nor its alias measures a volume: the first
  quantity is a sum of vertex displacements.

With this batch every *reaction* in `bozorg_etal_2014/scripts/fig3/BC` is
recognised. Those models are now blocked on the direction *block*
(`StaticDirection`, `ParallellDirection`) - a separate subsystem, not a
reaction.

### Direction block: `StaticDirection` + `ParallellDirection`, 2026-09-18

Not a reaction but the third block of a model file, and after the three
batches above it was the only thing still stopping `bozorg_etal_2014` from
loading: 43 rules across that repository, and every one of them is this pair.

The subsystem is small. A model declares at most one direction as two rules -
how it changes during the run, and what happens to it when a cell divides -
and the solvers already called `initiateDirection`/`updateDirection` as
structural no-ops, so this adds the two rule types, a registry for each, the
parse, and the division hook.

Both rules are no-ops, and so are legacy's. The direction is stored in
`cellData`, which division already copies to the daughter, so "both daughters
keep the mother's direction" needs no code. `ParallellDirection` does have one
branch that does work - re-picking each cell's "directional wall" - but that
list is only ever populated by `WallDirection`, a different update rule that
is not ported, so the branch is unreachable. Porting `WallDirection` is what
would require writing it, and the comment in `direction.cpp` says so.

Validated by showing the block changes nothing, which is the claim being made:

- In **both** binaries, the `growth_and_division` tutorial with the block
  appended is byte-identical to the same model without it.
- Comparing legacy against this build on `random.model` gives exactly the same
  residual (1.938e-02, same first index) with the block and without it. That
  residual is the model's random division direction, not the direction block:
  the two simulators draw different random streams and so cut cells
  differently. Adding the block contributes exactly zero.
- An unported rule (`WallDirection`) is still refused by name at read time, so
  a model that needs real direction machinery fails loudly instead of running
  with a direction that silently never moves. This is the part worth keeping
  an eye on if more rules are added.

The division hook is called from inside `divideCell`, after the topology is
final but before the size-dependent variable split - legacy's call site. It
makes no difference to these two no-op rules, but it is where a real rule
would need to be to see the geometry legacy gives it.

### The `shortestPath` benchmark outlier was a legacy hang, 2026-09-19

Worth writing down because the first reading of it was wrong, and the wrong
reading was flattering.

The benchmark sweep left three `growthanddivision/modelFiles/vertexBased/
shortestPath*` models uncompared: this build finished them in 3.5-16s and
legacy exceeded the timeout. Re-run with a 900s budget, legacy still did not
finish any of the three. That looks like a speedup of at least 57x, and it is
not a speedup at all.

What it is: README item 3, legacy's unbounded wall-pair orientation loop in
`Division::ShortestPath2D`. Established in this order, and each step matters:

- **Same answers where legacy terminates.** Raising `V_threshold` from the
  shipped 1.2 to 40 leaves only a few cells above it. There the two agree
  exactly - 0.000e+00 over 35953 values on `meristem.init` - and the cell
  counts match step for step (268 at t=2, 278 at t=20). So the division search
  itself is a faithful port; nothing about *what* is computed differs.
- **Comparable speed where legacy terminates.** At `V_threshold` 40 and 20
  legacy runs at 1.3x and 1.1x of this build. It is not slow at this.
- **A cliff, not a slope.** At 10 and below legacy never returns. A slope
  would be a performance difference; a cliff is a hang.
- **Confirmed directly.** A legacy process left on the shipped threshold of
  1.2 ran for 16 hours without reaching its second print. `sample` put 100%
  of samples inside `getCandidates`, in the `do/while(flippedVectors)` loop.
  This build does the same run in 0.88s, reaching 7680 cells from 267.

The loop cannot terminate when the division point lies on a wall's supporting
line: both re-orientation tests are strict inequalities, so neither fires at
exactly zero, and the third test then swaps the two walls and repeats on the
same pair. v2 bounds it at 8 iterations and skips the pair.

Three things were wrong in how this was first reported and are now fixed:

- The throwaway classifier called these three "v2 gaps". They are the
  opposite: this build runs them and legacy does not.
- `summarise.py` lumped them in with genuine failures. It now buckets
  not-compared rows by *which side* failed, and names the hang case
  explicitly so no one reads a ratio into it. The same pass fixed a worse
  mislabel: 68 rows were being reported as "neither could run it" when
  `benchmark.py` had simply not bothered running legacy after this build
  failed, and recorded a sentinel rather than a real error.
- `benchmark.py` itself was never at fault - it only divides when both times
  exist. The ">57x" was prose, not arithmetic.

The general lesson for this corpus: a legacy timeout is never evidence of a
speedup until legacy has been shown to terminate on a smaller version of the
same problem.

### The `VertexNoUpdate*` clamps, 2026-09-19

Six of the family: `FromPosition`, `FromIndexHoldX/Y/Z`, `FromList` and
`Boundary`. `FromIndex` was already ported. `VertexNoUpdateBoundary` is the
one the benchmark wanted (3 models in `bozorg_etal_2016`); the rest are its
near neighbours and came cheaply.

**Not** ported, and named here so the gap is explicit: the four geometric
variants `BoundaryPtemplate`, `BoundaryPtemplateStatic`,
`BoundaryPtemplateStatic3D` and `Boundary3D`. Those do not clamp a coordinate
but project out the component of the derivative normal to the boundary edge,
using the cell's stored centre of mass, and the Static ones cache the normals
at `initiate`. They are a different and much larger job, and they block
nothing in the corpus.

Exact against legacy over 16 configurations in 2D and 3D
(`tests/port/constraints.sh`), with one annotated exception below.

The fixture is the whole story here, twice over.

First, the signature. `VertexNoUpdateFromList` takes **no** variable indices;
the port initially demanded one, and every model using it would have been
rejected. The harness caught it immediately because legacy refused the model
outright and produced no output at all - which is worth noting as the easiest
kind of mismatch to read, and the reason to compare against legacy rather
than to eyeball plausible-looking numbers.

Second, and worse: the first version of the harness was **entirely vacuous**.
A clamp only zeroes derivatives, so it needs something pushing the vertices,
and the natural pairing is a wall spring reading wall variable 0. But
variable 0 holds each wall's resting length, so on a seeded init the springs
are already at equilibrium, no vertex moves, and all sixteen cases "matched"
legacy while clamping a tissue that was going nowhere. Pointing the spring at
wall variable 1 - a seeded value unrelated to the geometry - puts the springs
out of equilibrium. Each clamped run is now checked against the unclamped one
and against a neighbouring configuration, so the axis, the direction and the
axis-subset all have to change the output.

That change immediately surfaced something the vacuous version could not:
`FromPosition` with direction +1 and threshold 1.5 differs from legacy in
exactly one printed value, one unit in the last of six significant digits
(0.131243 against 0.131242), in a wall variable rather than a vertex. It is a
print-rounding boundary and not a divergence: the same comparison is exact at
solver tolerances 1e-5, 1e-7 and 1e-9, and only the harness's 1e-6 lands on
it. That case runs at 2e-6 with the reason written beside it.

One place this build is deliberately louder than legacy: `HoldZ` on a 2D
tissue indexes past the end of the vertex row in legacy. This throws instead.

### `CenterTriangulation::VertexFromCellPressure` and `...Linear`, 2026-09-19

The centre-triangulated turgor pair from `mechanical.cc`;
`...PressureLinear` is the benchmark blocker (3 models in `bozorg_etal_2016`).
Unlike the 2D forms already ported, these treat a cell as a fan of triangles
around a stored centre and push each edge outwards, perpendicular to the edge
and within its triangle's plane. 3D only. The two share that geometry, so it
is one helper and the reactions differ only in how the pressure is computed:
a constant optionally scaled by a concentration, or a linear ramp from zero.

Exact (0.000e+00) against legacy over 11 configurations
(`tests/port/pressure_ct.sh`), covering both volume-normalisation settings,
concentration scaling, deflation, two ramp rates against the run length, both
legacy aliases, and the pair working against a wall spring.

Two legacy quirks kept:

- A concentration index of **0** means "no concentration", so cell variable 0
  cannot be used as the concentration. A model written against legacy that
  names index 0 means "constant", and this build reads it the same way.
- Both vertices of an edge receive the *whole* force rather than half each.

One thing dropped rather than reproduced: legacy's `update()` accumulates the
elapsed time in a function-level `static`, shared by every instance of the
class. It only ever feeds a branch that was disabled with `if (true)`, so it
changes nothing; the ramp here is a plain member.

Three fixture problems, all of the same family as the previous batches:

- The concentration case was **vacuous** on the obvious init. In
  `twoSquare3D.init` cell variable 3 is 1.0 in both cells, so scaling by it
  multiplies by one; the case passed while testing nothing. `ct3d.init` sets
  0.4 and 0.9 instead.
- Deflation is the ill-conditioned direction. The cell shrinks towards its own
  centre, so the centre-to-midpoint vector the force is built from shrinks
  towards zero and normalising it amplifies everything. At K = -0.8 the cell
  has collapsed by t = 2 and this build needs 176k accepted and 63k
  **rejected** steps to get there - the error controller failing repeatedly is
  the tell. The two agree exactly at 1e-12 up to t = 1 and then part, so this
  is the solver, not the arithmetic. The harness uses -0.1, which keeps the
  geometry valid for the whole run and still exercises the sign.
- The spring case needed a tighter solver for the same underlying reason. At
  the harness's 1e-6 error tolerance the answer is only good to about 1e-6, so
  once the two step sequences part company the results differ by roughly that
  much (measured: 4.7e-6). At 1e-9 they agree to 1e-12.

The general point, worth keeping: when a comparison fails, check whether the
solver tolerance is looser than the difference before suspecting the port.

### `Pressure3D::Triangular`, 2026-09-19

Pressure on a triangular face, pushing all three vertices along the face
normal. The last of the benchmark's per-reaction blockers (1 model in
`bozorg_etal_2014`). Exact against legacy over 8 configurations, appended to
`tests/port/pressure_ct.sh`; it needs triangular faces, so those cases use
`tri3D.init` rather than the square mesh the rest of that file uses.

Mostly dead code in legacy, reproduced as written because the coefficients
are what the published runs used:

- The six-parameter form (areaFlag 2 and 3) was meant to raise the pressure
  as the enclosed volume fell. The line accumulating that volume is commented
  out, so `totalVolume` is identically zero and the whole expression collapses
  to a constant, `(Vfactor - Pfactor) * k_force / (Vfactor - 1)`. `V0` is
  therefore unused.
- `Zplane` is unused in every form; its only use is commented out too.
- Each branch opens with a loop that computes and normalises every cell's
  normal and then does nothing with it. Its only surviving effect is
  validating that the cells are triangular, which this build does directly.

Both vertices - all three here - take the whole coefficient rather than a
share, as in the centre-triangulation pair above.

One place this build is stricter: legacy accepts areaFlag 2 or 3 with three
parameters and then reads `parameter(4)` and `parameter(5)` out of range.
This rejects that combination at construction instead.
