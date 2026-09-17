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
