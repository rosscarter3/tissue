#!/usr/bin/env python3
"""Run a model through the legacy simulator and this one, and compare.

Every reaction ported from `legacy/` should be checked against the behaviour it
is replacing, not just against whether it compiles. This drives both binaries
over the same model/init/solver and diffs the print-flag-0 output numerically.

    python3 tools/port/compare.py MODEL INIT SOLVER [--tol 1e-9] [--block=vertex]

--block=vertex compares only the vertex-position blocks. Use it to check a
reaction's derivs() independently of an update() that writes cell variables:
port and validate the forces first, then the bookkeeping.

Exits non-zero if the two disagree by more than the tolerance, printing the
worst offending value and where it is.

Reactions that are stochastic, or that legacy integrates through a path with a
known defect (see "Deliberate fixes over legacy" in the README), will not match
to machine precision and should be compared with a loose tolerance and a note.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LEGACY = os.path.join(ROOT, "bin", "simulator")
NEW = os.path.join(ROOT, "build", "simulator")


def text_solver(solver):
    """A copy of the solver file that prints numbers to stdout.

    Tutorials often select VTK or gnuplot output, which writes files and leaves
    stdout empty - nothing to compare. The third non-comment line of every
    solver file is `printFlag numPrint`, so forcing the flag to 0 gives the
    plain text dump for any of them without touching the trajectory.
    """
    out, seen = [], 0
    for line in open(solver):
        body = line.split("#", 1)[0].strip()
        if body:
            seen += 1
            if seen == 3:
                parts = body.split()
                n = parts[1] if len(parts) > 1 else "20"
                out.append(f"0 {n}\n")
                continue
        out.append(line)
    path = os.path.join(os.path.dirname(NEW), "_compare_solver.rk5")
    open(path, "w").write("".join(out))
    return path


def run(binary, model, init, solver, extra=()):
    # No DYLD_LIBRARY_PATH: the legacy binary in bin/ is built against the
    # system libc++ and pointing it at conda's makes it fail to resolve
    # __ZNSt3__117bad_function_callD1Ev. Rebuild it with `CONDA_PREFIX= make`
    # if it ever does need one.
    p = subprocess.run([binary, model, init, solver, *extra],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"{os.path.basename(binary)} failed ({p.returncode}):\n"
                 f"{p.stderr[-800:]}")
    return p.stdout


def numbers(text):
    out = []
    for tok in text.split():
        try:
            out.append(float(tok))
        except ValueError:
            pass
    return out


def vertex_numbers(text):
    """Just the vertex-position blocks.

    Each printed state is [vertex block, cell block, wall block], every block
    a "<rows> <cols>" header followed by its rows. The vertex block is the
    first, so its header gives nVertex and the dimension; every later block
    with that exact header is another print's vertices.
    """
    lines = text.splitlines()
    hdr = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and all(p.isdigit() for p in parts):
            hdr, start = parts, i
            break
    if hdr is None:
        return []
    nv, dim = int(hdr[0]), int(hdr[1])
    out = []
    for i, line in enumerate(lines):
        if line.split() == hdr:
            for row in lines[i + 1:i + 1 + nv]:
                out += [float(t) for t in row.split()[:dim]]
    return out


def compare(a, b, tol, extract=numbers):
    na, nb = extract(a), extract(b)
    if len(na) != len(nb):
        return None, f"output shape differs: {len(na)} values vs {len(nb)}"
    worst, where = 0.0, -1
    for i, (x, y) in enumerate(zip(na, nb)):
        d = abs(x - y) / max(1.0, abs(x), abs(y))   # relative, floored at 1
        if d > worst:
            worst, where = d, i
    return worst, where


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tol = 1e-9
    for a in sys.argv[1:]:
        if a.startswith("--tol"):
            tol = float(a.split("=", 1)[1]) if "=" in a else float(
                sys.argv[sys.argv.index(a) + 1])
    if len(args) < 3:
        sys.exit(__doc__)
    model, init, solver = args[:3]
    extract = vertex_numbers if any(
        a.startswith("--block=vertex") for a in sys.argv[1:]) else numbers

    for b, name in ((LEGACY, "legacy/bin/simulator"), (NEW, "build/simulator")):
        if not os.path.exists(b):
            sys.exit(f"missing {name}; build it first")

    solver = text_solver(solver)
    old = run(LEGACY, model, init, solver)
    new = run(NEW, model, init, solver)
    worst, where = compare(old, new, tol, extract)
    if worst is None:
        print(f"MISMATCH  {where}")
        sys.exit(1)
    verdict = "match" if worst <= tol else "MISMATCH"
    print(f"{verdict}  worst relative difference {worst:.3e} "
          f"(tolerance {tol:.0e}) over {len(extract(new))} values"
          + (f", first at index {where}" if worst > tol else ""))
    sys.exit(0 if worst <= tol else 1)


if __name__ == "__main__":
    main()
