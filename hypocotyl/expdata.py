#!/usr/bin/env python3
"""Experimental data from Walia, Carter et al. (2024) Dev Cell 59:3245,
parsed straight from the authors' published analysis package
(zenodo.org/records/13379829, hook_model/data.py + length_measurements.csv)
so the numbers plotted are theirs, not transcribed."""
import ast, csv, os

_HERE = os.path.dirname(os.path.abspath(__file__))

def _find(name):
    """Locate one of the authors' data files.

    Default is the copy vendored under refdata/walia2024 (see its SOURCE.md);
    set WALIA_PKG to the unpacked Zenodo package to read that instead."""
    pkg = os.environ.get("WALIA_PKG")
    for cand in ([os.path.join(pkg, "hook_model", name)] if pkg else []) + [
            os.path.join(_HERE, "refdata", "walia2024", name)]:
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(f"{name}: not in refdata/walia2024 or $WALIA_PKG")


def _load():
    """Return {ClassName: {attr: literal}} from data.py without importing it
    (it pulls in pandas, which we do not need for the literals)."""
    src = open(_find("data.py")).read()
    out = {}
    for node in ast.parse(src).body:
        if not isinstance(node, ast.ClassDef):
            continue
        attrs = {}
        for st in node.body:
            if isinstance(st, ast.Assign) and isinstance(st.targets[0], ast.Name):
                try:
                    attrs[st.targets[0].id] = ast.literal_eval(st.value)
                except ValueError:
                    pass            # comprehensions / references, re-derived below
        if "t" not in attrs:
            # some classes build t with list(range(...)) or a comprehension
            for st in node.body:
                if (isinstance(st, ast.Assign)
                        and getattr(st.targets[0], "id", "") == "t"):
                    attrs["t"] = eval(compile(ast.Expression(st.value),
                                              "<t>", "eval"),
                                      {"list": list, "range": range}, {})
        if "t" in attrs and "t_h" not in attrs and max(attrs["t"], default=0) > 30:
            attrs["t_h"] = [v / 60.0 for v in attrs["t"]]   # minutes -> hours
        out[node.name] = attrs
    return out

DATA = _load()

def lengths_csv():
    """Tissue-length time course (Fig 1): fold change and growth rate."""
    rows = list(csv.DictReader(open(_find("length_measurements.csv"))))
    f = lambda k: [float(r[k]) if r[k] else None for r in rows]
    return {"t": f("time"), "inner": f("relative_inner"),
            "outer_short": f("relative_outer_short"),
            "outer_long": f("relative_outer_long"),
            "inner_rate": f("inner_growth_rate"),
            "outer_rate": f("outer_short_growth_rate")}

def series(d, key):
    """(t, y) for one series, truncated to their common length.

    The authors' series are not all the same length within a dataset -
    LightData carries dark to 28 h but the light conditions only to 10 h -
    so pair them explicitly rather than assuming a shared time axis."""
    y = d[key]
    n = min(len(d["t"]), len(y))
    return d["t"][:n], y[:n]

OPENING = DATA["OpeningData"]                 # angle mean + SE, 0-10 h
ISOX = DATA["IsoxabenData"]                   # cellulose synthesis inhibitor
ORYZALIN = DATA["OryzalinData"]               # microtubule depolymeriser
FOLD = DATA["LengthFoldChangeData"]           # inner/middle/outer +- std
LIGHT = DATA["LightData"]                     # dark / low light / normal light
YUC6 = DATA["YUC6Data"]                       # auxin overproduction (vs Col-0)
AUXINOLE = DATA["AuxinoleData"]               # auxin receptor antagonist

if __name__ == "__main__":
    for name, d in DATA.items():
        keys = [k for k in d if k != "t"]
        print(f"{name:22s} n_t={len(d.get('t', [])):3d}  series={keys}")
    L = lengths_csv()
    print(f"\nlength_measurements.csv: {len(L['t'])} timepoints, "
          f"inner {L['inner'][0]:.2f} -> {max(x for x in L['inner'] if x):.2f}")
    print("FOLD t_h:", [round(x, 2) for x in FOLD.get("t_h", [])])
