# Tissue v2

A ground-up rewrite of the Tissue vertex-based tissue simulator in modern
C++20, with hardware acceleration and full backwards compatibility with the
legacy file formats and command line.

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
spring mechanics, RK5Adaptive, identical input files, Apple M1):

| job | legacy | v2 | speedup |
|---|---|---|---|
| t=0..5, eps=1e-5 | 137.2 s (2502 steps) | 0.28 s (14 steps) | ~490x |
| t=0..50, eps=1e-8 | killed after 865 s, unfinished | 1.2-1.4 s | >600x |
| meristem tutorial (divisions, ~400 cells) | 4.7 s | 1.5 s | ~3x |

The dominant factor at scale is the legacy RK5 wall-derivative bug (below):
its corrupted embedded error estimate forces ~180x more solver steps at the
same tolerance on any model with wall dynamics. On top of that, v2's
per-step cost is ~4x lower at 40k cells (flat memory layout + fused kernels),
with threading adding ~20% more on this memory-bandwidth-bound machine.
Models without wall dynamics take the same steps as legacy and see only the
per-step gains.

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

## What is ported so far

Reactions: `Creation::{Zero,One,Two,Three,SpatialSphere}`,
`Degradation::{One,Two,N}`, `MassAction::{TwoToOne,General}`,
`DiffusionSimple`, `MoveVertexRadially`,
`WallGrowth::{Constant,Stress,Strain}`,
`WallGrowth::CenterTriangulation::Constant`, `WallMechanics::Spring`,
`CenterTriangulation::{Initiate,EdgeSpring}`, `CenterCOM`,
`CenterCOMcenterTriangulation`, `InitiateWallLength`, `Initiation::Random`,
`SisterVertex::{InitiateFromDistance,CombineDerivatives,Spring}` — plus their
legacy aliases.

Compartment changes: `Division::VolumeRandomDirection`,
`Division::ShortestPath2D`, `RemovalOutsideRadius`.

Solvers: `RK5Adaptive`, `RK4`, `Euler`, `HeunIto`.

Print flags: 0, 1, 2 (VTK), 3, 4, 5 (gnuplot), 77, 107.

This covers every model shipped in `examples/tutorials`. A model using an
unported reaction/rule fails fast with a clear message naming it. The legacy
catalog (~150 further reaction classes, direction machinery, PLY output,
special-purpose print flags) can be ported incrementally: each reaction is a
self-contained class registered with `TISSUE_REGISTER_REACTION`, typically
20-60 lines.

## Architecture

```
v2/src/
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
