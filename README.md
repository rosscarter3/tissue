# Tissue

A ground-up rewrite of the [Tissue](https://gitlab.com/slcu/teamHJ/tissue)
vertex-based tissue simulator in modern C++20, with hardware acceleration and
full backwards compatibility with the legacy file formats and command line.

This is a fork. The original sources are preserved verbatim under `legacy/`
(and its README as `README-legacy.md`): they are the reference against which
the reaction port is checked, and the port is not finished — see
*Reaction coverage* below.

## Building

Requires CMake >= 3.20 and a C++20 compiler (tested with Apple clang 17).

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Note for conda users: make sure conda's `CXXFLAGS`/`LDFLAGS` are not exported
into the environment (they mix incompatible libc++ versions on macOS).

The simulator binary lands in `build/simulator` and is invoked exactly like
the legacy one:

```sh
build/simulator modelFile initFile solverParaFile [-init_output file] \
    [-init_output_format tissue] [-verbose 0|1] [-vtk_output dir] \
    [-centerTri_init]
```

## Backwards compatibility

Unchanged and verified against the legacy simulator:

- **Model files**: same grammar (`#` comments, whitespace tokens, reaction /
  compartment-change blocks), same reaction names including legacy aliases
  (`VertexFromWallSpring` = `WallMechanics::Spring`, `MoveVertexRadially` =
  `GrowthForce::Radial`, ...).
- **Init files**: same reader and writer (`-init_output` round-trips),
  including the `-centerTri_init` variant.
- **Solver files**: same four solvers (`RK5Adaptive`, `RK4`, `Euler`,
  `HeunIto`) with identical parameter layouts and print schedules.
- **Outputs**: VTK `.vtu`/`.pvd` files are byte-identical to legacy; init
  output and the text print formats (flags 0, 3, 4, 5, 77, 107) match.
- **Random numbers**: the legacy `ran3` generator was ported bit-exactly
  (including its default seed and quirks), so seeded reactions reproduce the
  same sequences.

Validation results on this machine (Apple M1):

| test | result |
|---|---|
| `mass_action` tutorial, RK5Adaptive, 15 time units | matches legacy to <= 6e-15 relative |
| `mesh/growth` model, Euler 10k steps | **bit-identical** to legacy |
| VTK output (`.pvd`, cell + wall `.vtu`) | **byte-identical** to legacy |
| meristem growth+division (705 divisions, 557 removals) | statistically consistent endpoints; topology checked every step |

## Hardware acceleration

- All state lives in flat, contiguous tables (`Matrix`) instead of nested
  `std::vector<std::vector<double>>` — cache friendly and autovectorizable
  (built with `-mcpu=native`).
- A dependency-free persistent thread pool (`TISSUE_NUM_THREADS` to override,
  `=1` to disable) parallelizes the derivative kernels, the Runge-Kutta stage
  combinations, and the error norms. Reductions are deterministic for a fixed
  thread count; small models run on the exact serial path.
- Force accumulation (springs) uses a two-pass gather at scale: per-wall
  coefficients, then per-vertex summation — race-free and deterministic for
  *any* thread count.

Measured on a generated 40 000-cell / 80 400-wall tissue (wall growth +
spring mechanics, Apple M1, all numbers re-measured 2026-09-14):

| job | legacy | v2 | speedup |
|---|---|---|---|
| Euler, 1000 fixed steps (identical work) | 18.4 s | 0.61 s | **30x** |
| RK5Adaptive t=0..5, eps=1e-5 | 104.7 s (2501 steps) | 0.24 s (13 steps) | **~430x** |
| RK5Adaptive t=0..50, eps=1e-8 | killed after 865 s, unfinished | 1.2-1.4 s | >600x |
| meristem tutorial (divisions, ~400 cells) | 4.7 s | 1.5 s | ~3x |

Two separable effects. **Raw throughput** is ~1.8x per derivative evaluation
(flat memory layout, fused kernels), and 30x per Euler step at 40k cells once
the redundant duplicate evaluation legacy performs each step is removed and the
per-step connectivity sweep is gated (below).
**Step count** is where the large factors come from: on any model with wall
dynamics the legacy RK5 wall-derivative bug (below) corrupts the embedded
error estimate, forcing 2501 steps where the correct scheme needs 13 — a
192x difference at the same tolerance. Models without wall dynamics take the
same steps as legacy and see only the throughput gain.

Threading adds nothing on this machine: these kernels are
memory-bandwidth-bound and repeated runs show 1-thread and 8-thread parity
even at 40k cells (3.9-4.4 s for the Euler-1000 benchmark, either way). What
matters much more is *not* threading small loops. The grain threshold — the
size below which a loop runs serially — was originally 1024 elements, which
is far below the real crossover: a parallel region costs a mutex, a
condition-variable broadcast and a join (tens of microseconds) against tens
of nanoseconds of work per element. On the 514-cell hook shell that put the
wall loops (1040 elements) and the cell table (8224) just over the line and
spent half the wall clock in the kernel:

| grain | real | user | sys |
|---|---|---|---|
| 1024 (was) | 58.2 s | 29.0 s | **29.1 s** |
| 8192 | 38.2 s | 24.5 s | 9.5 s |
| **65536 (now)** | **29.5 s** | 21.5 s | **1.0 s** |
| single-threaded | 29.9 s | 21.5 s | 1.2 s |

**1.5x to 2.6x on that model, with a bit-identical trajectory** (element-wise
loops cannot change results, and the one reduction — the scatter — is
separately guarded, see below). The spread is real and worth knowing: lock
contention gets much worse when the machine is otherwise busy, so the fix
helps least on an idle machine (2 h run: 262 s -> 170 s CPU, 1.5x) and most
when several conditions run in parallel, which is the normal case. Two
back-to-back A/B pairs of a 0.5 h run with four other jobs resident:

| grain | CPU (user + sys) | real |
|---|---|---|
| 1024 | 110.3 s (53.2 + 57.1) / 114.2 s (50.3 + 63.9) | 103.9 s / 85.3 s |
| 65536 | 41.8 s (40.5 + 1.3) / 40.1 s (39.2 + 0.8) | 48.2 s / 45.2 s |

i.e. 2.7x less CPU and 2.0x less wall clock, with system time down ~50x.
Independently confirmed on an unrelated workload (pavement-cell lobing, 37
cells / 2353 walls, 24 h simulated): ~270 s -> 149 s, 1.8x, with model
observables unchanged to four significant figures. Neutral at 40k cells.
Override with `TISSUE_GRAIN` to re-measure the crossover elsewhere.

The deterministic scatter (private per-partition copies of the target table,
reduced in partition order) carries a cost of O(threads x |target|) that does
not shrink as the loop is split, so it only pays when the loop is at least as
long as the target table is wide; it now checks that explicitly, and reuses
its scratch buffers instead of allocating a set per call. This is the real
reason wall-force threading shows little gain even at scale — a two-pass
gather over the vertex incidence lists has no such term and is the scalable
alternative.

Independent runs (conditions, parameter sweeps) parallelize perfectly as
separate processes.

**Transcendentals in the TRBS kernel.** Profiling the hook shell after the
grain fix put about a quarter of the force kernel in `tan` and `acos`. The
resting angles there enter only as cotangents, and as the sine and cosine of
one angle, while the law of cosines supplies the cosine directly - so the
round trip through an angle is avoidable:

    cot(acos(c)) = c / sqrt(1 - c^2),   cos(acos(c)) = c,   sin(acos(c)) = sqrt(1 - c^2)

all exact, and one `sqrt` in place of each transcendental. Back-to-back on a
0.5 h hook run: **37.2 s -> 26.2 s, 1.42x**, with the trajectory agreeing to
3.3e-6 degrees - round-off, as an algebraic identity should give. This is the
one place v2 departs from the legacy TRBS arithmetic; the departure is in
v2's favour, since the legacy form loses precision in the acos/tan round trip
near the clamp.

**The per-step connectivity sweep.** `checkConnectivity` is an O(cells + walls)
consistency check that ran after every accepted step. Topology can only change
when a division or removal rule fires, so for any model without compartment
changes - most mechanical models, the hook shell and the 40k-cell benchmark
included - it re-validated an unchanged mesh forever. It is now gated on
`checkCompartmentChange` reporting an actual change, so dividing models are
still checked on every step where they divide, which is where the bugs are.
Init files are still validated once at load.

On the 40k-cell benchmark that sweep *was* the benchmark: **3.72 s -> 0.61 s,
6.1x, bit-identical output.** On the 514-cell hook it is worth 3%.

## Deliberate fixes over legacy (documented divergences)

These legacy defects were **fixed**, not replicated. Results can therefore
differ from legacy for affected models — in v2's favor:

1. **RK5 Cash-Karp wall-stage bug** (`rungeKutta.cc:486`): legacy writes the
   stage-3 wall derivatives into the stage-2 slot (`ak2W`), clobbering stage 2
   and leaving stage 3 stale. Any model with wall dynamics (wall growth,
   dilution, ...) was mis-integrated by legacy RK5. v2 implements the correct
   scheme. (Cell and vertex variables were unaffected, which is why the
   mass-action comparison matches to machine precision.)
2. **RK4 vertex k3 bug** (`rungeKutta.cc:804`): legacy drops the k3 term from
   the vertex weighted sum. Fixed.
3. **`Division::ShortestPath2D` infinite loop**: when the division point lies
   exactly on a wall's supporting line, the legacy wall-pair orientation loop
   cycles forever — the shipped `shortestPath.model` tutorial reproducibly
   hangs the legacy simulator. v2 bounds the loop and skips the degenerate
   pair (with a warning).
4. **`WallGrowth::Strain` 4-parameter form**: legacy reads `parameter(4)` and
   `variableIndex(0,1)` out of range (undefined behavior). v2 skips the
   equilibrium gate in that form instead.
5. Redundant duplicate derivative evaluations per step in Euler/RK4 were
   removed (identical results for side-effect-free reactions).
6. File parsing validates sizes/headers and reports errors; the legacy binary
   was built with `-DNDEBUG` so all of its `assert` "validation" was compiled
   out and malformed files silently misparsed.
7. **`Pressure2D::AreaPotential` orientation**: legacy's force direction
   depended on the arbitrary orientation the sorting pass happened to pick for
   each cell, so `P > 0` deflated cells sorted clockwise. v2 orients the force
   by the signed area, so positive pressure always inflates.
8. **`Bending::NeighborCenter` wrong vertex** (`bending.cc:72`): legacy
   accumulates into `vertexDerivs[k]`, where `k` is the cell-local wall
   counter, not the global vertex index — so every cell pushed on vertices
   `0..numWall-1` of the whole tissue regardless of which vertices the angle
   belonged to. v2 writes to the intended vertex.
9. **`Bending::Angle` wrong vertex** (`bending.cc:158`): legacy adds the
   central-vertex term to `vertexDerivs[jm]` a second time instead of to
   `vertexDerivs[j]`. The turning vertex therefore felt no force, its
   predecessor felt two, and the three contributions did not sum to zero, so
   the reaction injected net momentum. v2 writes to `j`; the forces are then
   exactly `-dE/dx` for `E = (k/2)(theta - theta_0)^2`.
10. **`SisterVertex::SpringCellConc` gate** (`sisterVertex.cc:297`): legacy's
    activity test reads `cellData[cell1] > 0 || cellData[cell1] > 0`, checking
    the first cell twice, so a pair was inert whenever the first vertex's cell
    was zero however strongly the second expressed — making the force depend on
    the order the pair was listed in. v2 tests `cell2` in the second clause, as
    legacy's own comment says it intends.

11. **`Diffusion::2D` stale cell centres**: legacy reads the two cell centres
    from `Cell::positionFromVertex()`, the no-argument overload, which returns
    each vertex's *cached* position rather than the live `vertexData` used for
    the wall length and cell area in the same expression. The cache is only
    refreshed by `BaseSolver::setTissueVariables()`, i.e. at print points, so
    on a moving mesh the centre-to-centre distance lagged the rest of the
    geometry by up to a print interval. v2 reads all three from `vertexData`
    (measured: identical on a static mesh, 2.5% divergence once the vertices
    move).

12. **Four broken cell1/cell2 symmetries in `membraneCycling.cc` and
    `membraneCyclingAll.cc`**. Which of a wall's two cells an init file stores
    as `cell1` is arbitrary, and the paired wall variables are stored in that
    same order, so relabelling both together must leave the dynamics
    unchanged. Four reactions break that, each through a typo in one of the
    two mirror-image branches:
    - `MembraneCycling::LocalWallFeedbackLinear` and
      `MembraneCycling::PINFeedbackLinear` index the **wall** table with the
      **cell** index (`wallData[i]` for `wallData[j]`) in their off-rate term,
      reading an unrelated wall - or out of range in any tissue with more
      cells than walls.
    - `MembraneCycling::PINFeedbackNonLinear` raises the membrane load to
      `parameter(2)`, the Hill half-max, instead of `parameter(3)`, the Hill
      exponent, so one face of every wall ran a different power law from the
      other.
    - `MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition` reads the
      wall signal from the first face in one term of its second branch and
      from the second face in the other.

    All four are fixed. They are validated by the symmetry itself rather than
    against legacy: `tests/port/membrane_swapcheck.py` reruns a model with
    every wall's orientation flipped and checks the cell variables are
    unchanged. v2 is invariant to 0.000e+00 for all four; legacy shifts by
    0.5-10%. The same check passes for legacy on the ten reactions in those
    files that needed no fix, so it is discriminating rather than trivially
    satisfied.

13. **`VertexFromTRBScenterTriangulationMT` neighbour weighting never
    happens** (`mechanicalTRBS.cc:5568`): the guard on the neighbour list is
    `neighbor[nn] < numCells && neighbor[nn] > -1`, and `neighbor` is a
    `std::vector<size_t>`, so the second test compares an unsigned value
    against `SIZE_MAX` and is **always false**. The whole neighbour-averaging
    loop is dead. What `parameter(5)` actually does in legacy is scale each
    cell's own stress by `(1 - neighbourweight)` before the anisotropy is
    taken from it - so a model asking for spatial smoothing of the stress
    field the microtubules respond to got no smoothing and a silently weaker
    stress instead. Measured on a fixture whose two cells agree to under 1%,
    where a real average must barely move the answer:

        weight   legacy       v2 (fixed)
        0.0001   0.00112498   0.00112509
        0.2      0.000900072  0.00112501
        0.4      0.000675054  0.00112497
        0.8      0.000225018  0.00112503

    Legacy tracks `(1 - w)` exactly. The diagnosis is confirmed the other way
    too: re-emulating the always-false guard reproduces legacy bit-for-bit.

14. **Unbounded Jacobi loops** (`mechanicalTRBS.cc`, five places): none of
    legacy's eigenvalue loops has an iteration bound, so a tensor that does
    not converge hangs the simulation rather than failing. It is reachable
    from ordinary parameter choices: `VertexFromTRBScenterTriangulationMT`
    with MF flag 2 and `Y_fibre > Y_matrix` gives a negative transverse
    modulus, and the legacy binary then runs indefinitely (still going after
    25 s on a two-cell mesh where the whole simulation takes well under a
    second). Every loop here is capped at 50 sweeps, which is far above what
    a symmetric 3x3 needs, so the run completes. That is a divergence only
    where legacy would not have terminated at all.

Item 13 applies to `VertexFromTRBSMT`,
`VertexFromTRLScenterTriangulationMT` and
`Hypocotyl3D::VertexFromTRBScenterTriangulationMT` as well: all four carry the
same dead guard, and all four are fixed.

Items 8-10, 12 and 13 mean those reactions cannot be compared bit-for-bit against
legacy, and should not be. `Bending::Angle` and `Bending::NeighborCenter` are
instead checked against an independent reimplementation of the force law in
`tests/port/bending_refcheck.py`; `SisterVertex::SpringCellConc` is compared
against legacy on the inputs the gate bug does not reach, and deliberately
differs on the others.

One legacy wart is **kept** rather than fixed: `Bending::{Angle,AngleInitiate,
AngleRelax}` write pi as the literal `3.14159`, and `Bending::Angle` clamps
`cos(theta)` to +/-0.999 so that a straight chain is not a stress-free state
(the residual is 0.045 rad). Both are load-bearing for models calibrated
against legacy and for rest angles stored in existing init files, so they are
documented in `src/tissue/reactions/bending.cpp` rather than changed.

## What is ported so far

Reactions: see `tools/port/STATUS.md`, which is generated from the sources by
`tools/port/inventory.py` and lists every legacy reaction still outstanding,
grouped by the legacy file that implements it. `tools/port/NOTES.md` records
how each batch was validated. Every reaction the shipped tutorials use is
ported, so everything in `examples/` runs.

Compartment changes: `Division::VolumeRandomDirection`,
`Division::ShortestPath2D`, `RemovalOutsideRadius`.

Solvers: `RK5Adaptive`, `RK4`, `Euler`, `HeunIto`, and `QuasiStatic` (new in
v2, no legacy counterpart): growth stepping with mechanical equilibrium solved
by FIRE relaxation rather than integrated. Explicit integration of overdamped
mechanics is bounded by the *stiffest* elastic mode while the modelled process
runs on the *softest* one — a ratio near 1e5 on the hook shell, where mean
`h` is 2e-5 h. Solving equilibrium instead removes that bound and the need to
scale Y and P for mobility. It does not reproduce `RK5Adaptive` on a model
tuned against the lagged dynamics; see the header of `quasi_static.cpp`.

A reaction that prescribes *velocity* rather than force — imposed tissue-level
dilation (`GrowthForce::Radial`/`MoveVertexRadially`), a moving boundary —
cannot be relaxed: the term never vanishes, so a naive relaxer drives it
without bound instead of converging. Such reactions declare themselves with
`Reaction::prescribesVelocity()` and supply `velocityDerivs()`; `QuasiStatic`
then integrates that contribution over the growth step and holds it out of the
relaxation. Setting a derivative to zero (a clamp, `VertexNoUpdateFromIndex`)
is not prescribed velocity — zero is consistent with force balance and needs
no declaration.

Print flags: 0, 1, 2 (VTK), 3, 4, 5 (gnuplot), 77, 107.

This covers every model shipped in `examples/tutorials`. A model using an
unported reaction/rule fails fast with a clear message naming it. The legacy
catalog (~150 further reaction classes, direction machinery, PLY output,
special-purpose print flags) can be ported incrementally: each reaction is a
self-contained class registered with `TISSUE_REGISTER_REACTION`, typically
20-60 lines.

## Architecture

```
src/
  tissue/core/       Tissue (index-based topology, geometry, division/removal
                     surgery, sorting, connectivity checking), Matrix (flat
                     CSR-backed state tables), legacy-exact ran3 RNG, registry
  tissue/reactions/  Reaction base + self-registering reaction library
  tissue/compartment/ CompartmentChange base + division/removal rules
  tissue/solvers/    RK5Adaptive, RK4, Euler, HeunIto + shared print machinery
  tissue/io/         comment-filtered readers, init/model grammar, VTK output,
                     CLI config
  tissue/parallel/   thread pool, deterministic parallel-for / scatter helpers
  apps/simulator.cpp CLI entry point (legacy-compatible)
```

Key design points versus legacy:

- Topology is index-based (no pointer webs); the background is the sentinel
  `kBackground` = `size_t(-1)`, matching the legacy `int(-1)` cast.
- The 970-line if/else factory chains are replaced by self-registering
  factories; adding a reaction touches exactly one file.
- The legacy `Cell::sortWallAndVertex` cyclic-ordering algorithm, the
  `divideCell` topology surgery, and swap-with-last removal semantics are
  ported faithfully (observable state matches, index assignment included).
- Errors are exceptions with actionable messages instead of `exit()` calls
  sprinkled through 100k lines.
```
