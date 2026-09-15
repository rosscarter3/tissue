# Tissue v1 → v2

A ground-up rewrite of the vertex-based tissue simulator in C++20 — what
changed, what it is measured to do, and the integrator defects it corrects in
the original.

Branch `v2-modern-rewrite` (local). Benchmarks and validation re-measured
14 September 2026 on Apple M1, both binaries built from the same working tree.

| | |
|---|---|
| Simulator source | 102,905 lines → **7,010 lines** |
| 40k-cell mechanical model, identical inputs | **430× faster** |
| Worst relative difference vs v1 on the reference tutorial | **5.4e-15** |
| Defects found in v1 | **6**, two affecting published-model numerics |

## 1. What this is

A complete reimplementation of the simulator core: topology, state, solvers,
reactions and I/O. It reads the same model, init and solver files as v1, takes
the same command line, and writes byte-identical VTK. No v1 code was carried
over, but the observable semantics were reproduced deliberately — the cyclic
wall/vertex sorting, the cell-division surgery, swap-with-last removal
ordering, and the `ran3` generator down to its default seed.

`v2/` sits alongside the untouched original; both build and run from the same
checkout, which is what made the validation below possible.

## 2. Architecture

| Concern | v1 | v2 |
|---|---|---|
| Topology | Pointer webs between `Cell`/`Wall`/`Vertex` | Index-based; background is `size_t(-1)` |
| State storage | `vector<vector<double>>` | Flat CSR-backed `Matrix`, one allocation |
| Reaction factory | 839-line if/else chain, 399 name strings | Self-registering; 38 classes, ~57 names |
| Adding a reaction | Header + source + factory chain | One file, one registration macro |
| Error handling | `exit()` at ~1,400 sites; asserts compiled out | Exceptions with actionable messages |
| Build | Makefile injecting `$CONDA_PREFIX` | CMake + `ctest` regression suite |
| Parallelism | none | Thread pool, deterministic reductions |

| Metric | v1 | v2 |
|---|---|---|
| Simulator source | 102,905 lines / 128 files | 7,010 lines / 48 files |
| Clean build, 8 cores | 16.8 s | 12.6 s |
| Binary size | 2.75 MB | 0.70 MB |
| Automated tests | none | 5 / 5 in 0.63 s |

The line-count gap is not all deletion: v2 implements the reaction set the
shipped tutorials and the hypocotyl work need, not the full v1 catalogue (§7).

## 3. Compatibility, and how it was checked

| Check | Result |
|---|---|
| `mesh/growth`, Euler, 10,000 steps | **bit-identical** |
| VTK output, cell + wall `.vtu` | **byte-identical** |
| `mass_action`, RK5Adaptive, 15 time units | 5.4e-15 worst relative |
| Meristem growth + division (705 divisions / 557 removals) | consistent endpoints, topology re-validated every step |
| Init round trip, incl. `-centerTri_init` | reproduces v1 format and precision |

The mass-action case matches to machine precision rather than exactly because
it exercises RK5Adaptive, where v2 deliberately diverges — see §4.1. Cell and
vertex variables are unaffected by that bug, hence agreement at 1e-15.

## 4. Defects found in v1

Not carried across. For affected models v2 will not reproduce v1's numbers.

1. **RK5 Cash–Karp wall stage is clobbered** (`rungeKutta.cc:486`) —
   *affects results.* The stage-3 call writes wall derivatives into the
   stage-2 buffer (`ak2W` for `ak3W`): stage 2 is overwritten before use,
   stage 3 holds the previous step's values. Any model with wall dynamics is
   integrated with a corrupted tableau, and the embedded error estimate it
   feeds is wrong — hence the step explosion in §5.
2. **RK4 drops the third vertex stage** (`rungeKutta.cc:804`) —
   *affects results.* `=` instead of `+=`, so k3 never enters the weighted
   sum for vertices.
3. **`Division::ShortestPath2D` can hang** — when a division point lies on a
   wall's supporting line the orientation loop cycles forever. This hangs v1
   on its own shipped `shortestPath.model` tutorial (confirmed by
   stack-sampling a run spinning for 17 CPU-minutes). v2 bounds the loop.
4. **`WallGrowth::Strain` reads out of range** — the 4-parameter form reads
   `parameter(4)` and `variableIndex(0,1)`, which it does not have.
5. **Redundant derivative evaluation every step** in Euler and RK4 — half of
   the 3.6× per-step gain is simply removing it.
6. **Validation compiled out** — `-DNDEBUG` disables every `assert`, and the
   readers rely on them exclusively, so malformed init files silently
   misparse in the shipped binary.

## 5. Benchmarks

Apple M1, same machine and input files. Test tissue: 200×200 grid = 40,000
cells, 80,400 walls, 40,401 vertices, wall growth + spring mechanics.

| Job | v1 | v2 | Gain |
|---|---|---|---|
| Euler, 1,000 fixed steps (identical work) | 18.4 s | 5.2 s | 3.6× |
| RK5Adaptive, t=0..5, eps 1e-5 | 104.7 s · 2,501 steps | 0.24 s · 13 steps | **430×** |
| RK5Adaptive, t=0..50, eps 1e-8 | killed at 865 s | 1.3 s | >600× |
| Meristem tutorial, ~400 cells with divisions | 4.7 s | 1.5 s | 3.1× |

Two separable effects. **Throughput** — flat memory layout, fused RK stage
kernels, no duplicate evaluation — is ~1.8× per derivative call and 3.6× per
Euler step. **Step count** produces the large factors: the v1 error estimate
is corrupted by the Cash–Karp defect, so it takes 192× more steps at the same
tolerance. Models without wall dynamics see only the throughput gain.

**Threading** is currently parity (5.2 vs 5.4 s at 40k cells; 122 vs 123 s at
350 cells): these kernels are memory-bandwidth-bound on this machine. The
pool's grain thresholds keep small models on the serial path. Independent runs
parallelise perfectly as separate processes.

### Tuning a real workload

3D apical-hook pipeline (2.5 h equilibration + 10 h light + 10 h dark, 350-cell
pressurised shell):

| Change | Pipeline | Per simulated hour |
|---|---|---|
| Starting point | 75 min | 390 s |
| Refresh anisotropic material every 0.02 h, not every step | 56 min | 292 s |
| Cache the Hill redistribution factor per cell | 47 min | 244 s |
| Solver tolerance 1e-4 → 1e-3 | **94 s** | **7 s** |

Profiling put 37% of runtime in the per-step stress-tensor and Jacobi pass, and
most of the rest in steps resolving elastic transients the growth-driven
trajectory never depends on. Verified rather than assumed: 1e-4, 1e-3 and 1e-2
give identical hook angles (162.7° → 116.1°) and fold changes over the full
run. 1e-3 ships, for margin.

Remaining structural cost: explicit integration needs steps proportional to
stiffness over growth rate — the very ratio that must be large for
quasi-staticity. A Newton or FIRE minimisation between growth steps (as the
published hook model used) is the principled next step, worth an estimated
further 5–20×.

## 6. New capability

- **3D TRBS shell mechanics** — `VertexFromTRBScenterTriangulation` (triangular
  biquadratic springs on centre-triangulated cells, giving walls a real 2D
  constitutive law), ported from the v1 force kernel with guarded `acos`/Heron
  evaluations so rejected trial steps cannot NaN the state. Paired with
  `Pressure3D::CenterTriangulation`.
- **Dynamic anisotropic wall material** — per-cell Cauchy stress tensor
  computed, diagonalised and stored; `WallMechanics::FiberSpring` redistributes
  the cellulose/CMT stiffness component between principal stress directions
  following Eq. 3 of Walia, Carter et al. (2024). Evaluated between steps, not
  inside derivative calls.
- **Supporting reactions** — `WallMechanics::BendingChain`,
  `CenterTriangulation::EdgeRelaxation`, `WallGrowth::StrainSaturationHill`,
  `WallGrowth::AcidGrowth`, `Pressure2D::AreaPotential`,
  `VertexNoUpdateFromIndex`.

These support a cell-resolved apical hook model (`v2/hypocotyl/`, see its own
README for the full validation). The 514-cell shell opens 159° → 13° over 10 h
under light against a maintained dark control, with a hook-angle RMSE of 11.9°
against the paper's measurements and inner/outer tissue extension of 2.16× /
1.11× (measured 2.15× / 1.06×). It independently reproduces the paper's
stress-anisotropy signature — circumferential principal stress on both flanks
in darkness, switching to longitudinal on the outer side within the first hour
under light while the inner side never switches — and the auxin- and
pH-dependent perturbations (YUC6-OX, low light). It does not reproduce the
microtubule and cellulose perturbations: the fibre stiffness carries too
little load relative to the matrix modulus for Eq. 3 to affect the mechanics,
so anisotropy is a faithful readout there rather than a driver.

## 7. Not ported

A model using any of these fails immediately with a message naming the missing
rule, rather than silently doing something else.

- The long tail of the v1 reaction catalogue (~150 classes, largely
  auxin-transport and PIN-polarity variants). Each is a self-contained 20–60
  line class in v2; porting is mechanical.
- The direction machinery (`UpdateMTDirection`, direction update/division).
  Models declaring a direction block are rejected at read time.
- The microtubule-feedback TRBS variants.
- PLY output, the optimiser, and paper-specific print flags. Flags 0–5, VTK
  1–2, 77 and 107 are implemented.

## 8. Building and running

```sh
cd v2
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # 5 tutorial regressions

build/simulator model.txt init.txt solver.rk5 \
    [-init_output out.init] [-init_output_format tissue|centerTriTissue] \
    [-vtk_output dir] [-centerTri_init] [-verbose 0|1]
```

Command line and file formats are unchanged from v1. On macOS with conda
active, keep conda's `CXXFLAGS`/`LDFLAGS` out of the environment — they mix
incompatible libc++ versions; this affects the v1 Makefile too, which injects
`$CONDA_PREFIX` itself.
