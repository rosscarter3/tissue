#!/usr/bin/env python3
"""Time the legacy simulator against this one on real published models.

Discovers runnable (model, init, solver) triples in a directory of cloned
model repositories and runs each through both binaries, reporting wall time.

    python3 tools/bench/benchmark.py BENCHDIR [--timeout 120] [--out r.json]
                                     [--keep-output] [--filter SUBSTR]

Both binaries get the same solver file, with the print flag forced to 0 and
the print count to 2 unless --keep-output is given: most of these models are
configured to write VTK, and that I/O would otherwise dominate and drown the
compute difference we are trying to measure. Use --keep-output to time the
shipped configuration instead.

A run that exceeds the timeout is reported as such rather than compared: if
legacy times out and this one does not, the speedup is a lower bound, which
is still worth recording. Models whose reactions are not ported yet fail at
read time and are reported as unsupported, which doubles as a coverage check.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LEGACY = os.path.join(ROOT, "bin", "simulator")
NEW = os.path.join(ROOT, "build", "simulator")
SOLVER_EXT = (".rk5", ".rk4", ".rk", ".solver", ".euler")


def rewrite_solver(path, dest, fraction, quiet=True):
    """Copy the solver file, optionally silenced and shortened.

    Line 2 of a solver file is "startTime endTime" and line 3 is
    "printFlag numPrint". Scaling the simulated duration by `fraction` is what
    makes this a benchmark rather than an endurance test: several of these
    models run for minutes, and both binaries need to do the *same* work for a
    ratio to mean anything, so truncating one on a timeout is useless.
    """
    out, seen = [], 0
    for line in open(path, errors="ignore"):
        body = line.split("#", 1)[0].strip()
        if body:
            seen += 1
            if seen == 2 and fraction != 1.0:
                parts = body.split()
                try:
                    t0, t1 = float(parts[0]), float(parts[1])
                    out.append(f"{t0} {t0 + (t1 - t0) * fraction}\n")
                    continue
                except (ValueError, IndexError):
                    pass
            if seen == 3 and quiet:
                out.append("0 2\n")
                continue
        out.append(line)
    open(dest, "w").write("".join(out))
    return dest


def cell_count(init):
    """First integer on the first non-comment line is the cell count."""
    for line in open(init, errors="ignore"):
        body = line.split("#", 1)[0].split()
        if body:
            try:
                return int(body[0])
            except ValueError:
                return None
    return None


def pair_init(model, inits):
    """Prefer an init sharing a name prefix with the model, else the first."""
    stem = os.path.splitext(os.path.basename(model))[0].lower()
    best, score = inits[0], -1
    for i in inits:
        istem = os.path.splitext(os.path.basename(i))[0].lower()
        n = len(os.path.commonprefix([stem, istem]))
        if n > score:
            best, score = i, n
    return best


def discover(benchdir):
    """(model, init, solver) triples, one init and one solver per model."""
    triples = []
    for dirpath, _, filenames in os.walk(benchdir):
        if os.sep + ".git" in dirpath:
            continue
        models = sorted(f for f in filenames if f.endswith(".model"))
        if not models:
            continue
        # Inits and solvers may sit beside the models or one level up.
        for where in (dirpath, os.path.dirname(dirpath)):
            inits = sorted(f for f in os.listdir(where) if f.endswith(".init"))
            solvers = sorted(f for f in os.listdir(where)
                             if f.endswith(SOLVER_EXT))
            if inits and solvers:
                break
        else:
            continue
        inits = [os.path.join(where, f) for f in inits]
        solvers = [os.path.join(where, f) for f in solvers]
        for m in models:
            triples.append((os.path.join(dirpath, m),
                            pair_init(m, inits), solvers[0]))
    return triples


def run(binary, model, init, solver, cwd, timeout):
    """Wall time in seconds, or ('timeout'|'fail', detail)."""
    args = [binary, model, init, solver]
    start = time.monotonic()
    try:
        p = subprocess.run(args, cwd=cwd, capture_output=True, text=True,
                           timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, "timeout"
    elapsed = time.monotonic() - start
    if p.returncode != 0:
        err = (p.stderr or "").strip().replace("\n", " ")[-160:]
        return None, err or f"exit {p.returncode}"
    return elapsed, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("benchdir")
    ap.add_argument("--timeout", type=float, default=120.0)
    ap.add_argument("--target", type=float, default=1.5,
                    help="aim for runs of at least this many seconds")
    ap.add_argument("--out", default=None)
    ap.add_argument("--keep-output", action="store_true")
    ap.add_argument("--filter", default="")
    ap.add_argument("--classify", action="store_true",
                    help="when this build cannot run a model, try legacy too, "
                         "so the failure is attributed to one side or both")
    ap.add_argument("--classify-timeout", type=float, default=20.0,
                    help="budget for those legacy attempts (default 20s)")
    args = ap.parse_args()

    for b, n in ((LEGACY, "bin/simulator"), (NEW, "build/simulator")):
        if not os.path.exists(b):
            sys.exit(f"missing {n}; build it first")

    triples = [t for t in discover(args.benchdir) if args.filter in t[0]]
    print(f"{len(triples)} candidate models under {args.benchdir}\n")
    results = []
    tmp = tempfile.mkdtemp(prefix="tissuebench")
    for model, init, solver in triples:
        rel = os.path.relpath(model, args.benchdir)
        cwd = os.path.dirname(model)
        row = {"model": rel, "init": os.path.basename(init),
               "solver": os.path.basename(solver), "cells": cell_count(init)}
        # Calibrate: shorten the run until this build finishes it quickly,
        # then time both binaries on that same shortened solver.
        used, fraction, tnew, enew = solver, 1.0, None, None
        for fraction in (0.01, 0.05, 0.25, 1.0):
            used = rewrite_solver(solver, os.path.join(tmp, "solver.q"),
                                  fraction, not args.keep_output)
            tnew, enew = run(NEW, model, init, used, cwd, args.timeout)
            if tnew is None:
                break          # failed or timed out; no point lengthening it
            if tnew >= args.target or fraction == 1.0:
                break
        row["fraction"] = fraction
        if tnew is None and enew != "timeout":
            if args.classify:
                # Run legacy anyway, briefly, so the failure can be attributed.
                # Several of these models are stale for legacy too - templates
                # with unsubstituted variable names, or reactions legacy itself
                # withdrew - and calling those coverage gaps overstates what is
                # left to port.
                told, eold = run(LEGACY, model, init, used, cwd,
                                 args.classify_timeout)
                if told is not None:
                    told, eold = None, "ran (this build could not)"
            else:
                told, eold = None, "skipped (this build cannot run the model)"
        else:
            told, eold = run(LEGACY, model, init, used, cwd, args.timeout)
        row.update(new=tnew, new_err=enew, legacy=told, legacy_err=eold)
        if tnew and told:
            row["speedup"] = told / tnew
        results.append(row)
        status = (f"{told:7.2f}s -> {tnew:7.2f}s  x{told/tnew:5.2f}"
                  if tnew and told else
                  f"legacy={eold or f'{told:.2f}s'}  new={enew or f'{tnew:.2f}s'}")
        print(f"  {rel[:56]:56s} {str(row['cells'] or '?'):>5} cells  {status}")
    shutil.rmtree(tmp, ignore_errors=True)

    ok = [r for r in results if r.get("speedup")]
    print(f"\n{len(ok)} of {len(results)} ran on both.")
    if ok:
        sp = sorted(r["speedup"] for r in ok)
        total_old = sum(r["legacy"] for r in ok)
        total_new = sum(r["new"] for r in ok)
        print(f"  speedup  min {sp[0]:.2f}  median {sp[len(sp)//2]:.2f}  "
              f"max {sp[-1]:.2f}")
        print(f"  total    {total_old:.1f}s -> {total_new:.1f}s  "
              f"(x{total_old/total_new:.2f} overall)")
    if args.out:
        json.dump(results, open(args.out, "w"), indent=1)
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
