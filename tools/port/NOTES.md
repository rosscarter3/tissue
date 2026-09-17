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
