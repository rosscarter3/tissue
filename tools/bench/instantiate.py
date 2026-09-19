#!/usr/bin/env python3
"""Turn tissue parameter-scan templates into runnable model files.

Several published model directories ship a template rather than a model: the
.model file has a bare name like `young` or `FORCE` where a number belongs,
and a sibling .auto file describes the sweep a scan driver was meant to run:

    VAR young FROM 0 TO 160 STEPS 61 LINEAR
    FILE model.model

Neither simulator can read the template - legacy fails on it too - so these
directories contribute nothing to a benchmark despite holding real models.
This substitutes one value per variable and writes a runnable copy, which is
enough to time the model even though it is not the published sweep.

    python3 tools/bench/instantiate.py BENCHDIR [--at 0.5] [--suffix .inst]

--at picks where in each variable's range to sample, 0 = FROM and 1 = TO;
the default is the midpoint. A range starting at zero often means "the
feature is off at one end", so the midpoint is usually the interesting case.
"""
import argparse
import os
import re
import sys


def parse_auto(path):
    """{varname: (lo, hi)} and the model filename the template lives in."""
    variables, model = {}, None
    for line in open(path, errors="ignore"):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "VAR" and len(parts) >= 5 and parts[2] == "FROM":
            try:
                variables[parts[1]] = (float(parts[3]), float(parts[5]))
            except (ValueError, IndexError):
                pass
        elif parts[0] == "FILE" and len(parts) >= 2:
            model = parts[1]
    return variables, model


def substitute(text, values):
    """Replace each bare variable name with its value. Returns (text, n)."""
    total = 0
    for name, value in values.items():
        # Whole word only: `young` must not match inside `young_fiber`, and a
        # name appearing in a trailing comment is left alone by doing the
        # replacement on the code part of each line only.
        pattern = re.compile(rf"(?<![\w.]){re.escape(name)}(?![\w.])")
        out = []
        for line in text.split("\n"):
            code, sep, comment = line.partition("#")
            code, n = pattern.subn(repr(value), code)
            total += n
            out.append(code + sep + comment)
        text = "\n".join(out)
    return text, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("benchdir")
    ap.add_argument("--at", type=float, default=0.5,
                    help="where in each range to sample, 0..1 (default 0.5)")
    ap.add_argument("--suffix", default=".inst",
                    help="inserted before .model in the output name")
    args = ap.parse_args()

    written = 0
    for dirpath, _, files in os.walk(args.benchdir):
        for name in files:
            if not name.endswith(".auto"):
                continue
            auto = os.path.join(dirpath, name)
            variables, model = parse_auto(auto)
            if not variables or not model:
                print(f"  skip {auto}: no VAR/FILE lines")
                continue
            src = os.path.join(dirpath, model)
            if not os.path.exists(src):
                print(f"  skip {auto}: {model} not found")
                continue
            values = {v: lo + (hi - lo) * args.at
                      for v, (lo, hi) in variables.items()}
            text, n = substitute(open(src, errors="ignore").read(), values)
            if n == 0:
                print(f"  skip {src}: no occurrences of "
                      f"{', '.join(variables)} to substitute")
                continue
            stem = model[:-len(".model")] if model.endswith(".model") else model
            dest = os.path.join(dirpath, stem + args.suffix + ".model")
            open(dest, "w").write(text)
            written += 1
            shown = ", ".join(f"{k}={v:g}" for k, v in values.items())
            print(f"  {os.path.relpath(dest, args.benchdir)}  ({shown}, "
                  f"{n} substitution{'s' if n != 1 else ''})")
    print(f"\nwrote {written} runnable model(s)")
    return 0 if written else 1


if __name__ == "__main__":
    sys.exit(main())
