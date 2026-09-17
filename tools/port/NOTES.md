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
