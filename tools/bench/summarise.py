#!/usr/bin/env python3
"""Summarise benchmark.py output: speedups by repository and by mesh size.

    python3 tools/bench/summarise.py results.json [--md]

The "distinct cases" count is the one to read first. A published model
repository typically holds one model swept over many parameter sets or many
input meshes, so a run count says nothing about how varied the benchmark is:
grouping by (repository, model filename, element count) collapses those
sweeps and shows how many genuinely different shapes of work were timed.
Treat the per-size table as descriptive of whichever families happen to sit
in each bucket, not as a scaling curve.
"""
import json
import statistics
import sys


def bucket(n):
    if n is None:
        return "?"
    for lo, label in ((10000, "10k+"), (1000, "1k-10k"), (100, "100-1k")):
        if n >= lo:
            return label
    return "<100"


def main():
    path = sys.argv[1]
    md = "--md" in sys.argv
    rows = json.load(open(path))
    ok = [r for r in rows if r.get("speedup")]
    skipped = [r for r in rows if not r.get("speedup")]

    def table(title, groups):
        print(f"\n{title}")
        head = f"{'group':34s} {'n':>4} {'median':>7} {'range':>14} {'legacy s':>9} {'new s':>8}"
        print(head if not md else "| " + " | ".join(head.split()) + " |")
        for name, rs in groups:
            sp = sorted(r["speedup"] for r in rs)
            lo, hi = sp[0], sp[-1]
            print(f"{name:34s} {len(rs):4d} {statistics.median(sp):6.2f}x "
                  f"{lo:5.2f}-{hi:5.2f}x  {sum(r['legacy'] for r in rs):8.1f} "
                  f"{sum(r['new'] for r in rs):8.1f}")

    byrepo = {}
    for r in ok:
        byrepo.setdefault(r["model"].split("/")[0], []).append(r)
    table("By repository:", sorted(byrepo.items()))

    bysize = {}
    for r in ok:
        bysize.setdefault(bucket(r.get("cells")), []).append(r)
    order = ["<100", "100-1k", "1k-10k", "10k+", "?"]
    table("By mesh size (elements in the init file):",
          [(k, bysize[k]) for k in order if k in bysize])

    cases = {}
    for r in ok:
        key = (r["model"].split("/")[0], r["model"].split("/")[-1],
               r.get("cells"))
        cases.setdefault(key, []).append(r)
    print("\nDistinct cases (repository, model file, element count):")
    print(f"{'case':46s} {'runs':>5} {'median':>7}")
    for (repo, name, cells), rs in sorted(cases.items()):
        label = f"{repo}/{name} @ {cells}"
        med = statistics.median([r["speedup"] for r in rs])
        print(f"{label:46.46s} {len(rs):5d} {med:6.2f}x")

    sp = sorted(r["speedup"] for r in ok)
    print(f"\nOverall: {len(ok)} runs over {len(cases)} distinct cases.")
    print(f"  median speedup {statistics.median(sp):.2f}x, "
          f"range {sp[0]:.2f}x-{sp[-1]:.2f}x")
    tl = sum(r["legacy"] for r in ok)
    tn = sum(r["new"] for r in ok)
    print(f"  total wall time {tl:.0f}s -> {tn:.0f}s ({tl/tn:.2f}x)")

    if skipped:
        reasons = {}
        for r in skipped:
            why = r.get("new_err") or r.get("legacy_err") or "?"
            why = why.split(":")[-1].strip()[:70]
            reasons[why] = reasons.get(why, 0) + 1
        print(f"\n{len(skipped)} not compared:")
        for why, n in sorted(reasons.items(), key=lambda x: -x[1]):
            print(f"  {n:4d}  {why}")


if __name__ == "__main__":
    main()
