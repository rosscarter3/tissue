"""Build a tissue model file, and find out what reactions exist.

A model file is a list of reaction blocks, each one

    Name  numParameter  numLevel  countPerLevel...
    parameter values, one per line
    index values, one line per level

and the counts have to match what the reaction expects. They are not
documented anywhere outside the C++, so the usual way to find out is to run
it and read the error, which is a poor way to start.

    m = Model()
    m.add("Pressure2D::AreaElastic", [3.0, 30.0, 0.0134], [[2]])
    m.add("WallMechanics::Spring", [600.0], [[0]])
    m.write("my.model")

`reactions()` lists what the simulator actually has, read from its source,
and `describe()` gives one reaction's parameter names where the C++ names
them -- which it does, in the `configure` call, for most of them.
"""
from __future__ import annotations

import re

from .paths import SRC


class Model:
    """Reaction blocks in, model file out."""

    def __init__(self, comment: str = ""):
        self.comment = comment
        self.blocks: list[tuple[str, list[float], list[list[int]]]] = []

    def add(self, name: str, parameters=(), levels=()) -> "Model":
        """Add one reaction.

        `levels` is a list of index lists, one per level -- most reactions
        take one level naming the wall or cell variables they read and write.
        Passing a flat list of ints is taken as a single level, since that is
        what is meant nine times in ten.
        """
        levels = list(levels)
        if levels and all(isinstance(x, int) for x in levels):
            levels = [list(levels)]
        self.blocks.append((name, [float(p) for p in parameters],
                            [[int(i) for i in lv] for lv in levels]))
        return self

    def to_text(self) -> str:
        out = []
        for line in (self.comment or "").splitlines():
            out.append(f"# {line}")
        out.append(f"{len(self.blocks)} 0 0")
        for name, params, levels in self.blocks:
            counts = " ".join(str(len(lv)) for lv in levels)
            head = f"{name} {len(params)} {len(levels)}"
            out.append("")
            out.append(f"{head} {counts}".rstrip())
            out += [f"{p}" for p in params]
            out += [" ".join(str(i) for i in lv) for lv in levels]
        return "\n".join(out) + "\n"

    def write(self, path) -> str:
        with open(path, "w") as f:
            f.write(self.to_text())
        return str(path)

    def __repr__(self) -> str:
        return f"Model({len(self.blocks)} reactions)"


def _sources():
    if not SRC.is_dir():
        return []
    return sorted(SRC.rglob("*.cpp"))


def reactions(match: str = "") -> list[str]:
    """Every reaction the simulator registers, optionally filtered.

        reactions("Pressure")  ->  ['Pressure2D::AreaElastic', ...]
    """
    found = set()
    pat = re.compile(r'TISSUE_REGISTER_REACTION\(\s*\w+\s*,\s*"([^"]+)"')
    for p in _sources():
        found |= set(pat.findall(p.read_text(errors="ignore")))
    names = sorted(found)
    if match:
        low = match.lower()
        names = [n for n in names if low in n.lower()]
    return names


def _split_args(text: str) -> list[str]:
    """Top-level comma split, ignoring commas inside (), {} or "" ."""
    out, depth, cur, instr = [], 0, [], False
    for ch in text:
        if instr:
            cur.append(ch)
            if ch == '"':
                instr = False
            continue
        if ch == '"':
            instr = True
        elif ch in "({[":
            depth += 1
        elif ch in ")}]":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if cur:
        out.append("".join(cur).strip())
    return out


def describe(name: str) -> dict:
    """What one reaction expects: parameter names and index-level sizes.

    Read out of its `configure` call, which is where the C++ states both.
    The call is not always literal -- some reactions size their index levels
    from the file rather than fixing them -- so that field comes back as the
    expression when it cannot be reduced to numbers. Returns {} if the
    reaction is not registered.
    """
    needle = 'configure("' + name + '"'
    for p in _sources():
        text = p.read_text(errors="ignore")
        at = text.find(needle)
        if at < 0:
            continue
        start = text.index("(", at)
        depth, end = 0, None
        for k in range(start, len(text)):
            if text[k] == "(":
                depth += 1
            elif text[k] == ")":
                depth -= 1
                if depth == 0:
                    end = k
                    break
        if end is None:
            continue
        args = _split_args(text[start + 1:end])
        if len(args) < 4:
            continue
        nparam = args[3].strip()
        levels_raw = args[4].strip() if len(args) > 4 else ""
        names_raw = args[5].strip() if len(args) > 5 else ""
        nums = re.findall(r"\d+", levels_raw)
        levels = ([int(x) for x in nums]
                  if levels_raw.startswith("{") and levels_raw.endswith("}")
                  else levels_raw)
        pnames = re.findall(r'"([^"]+)"', names_raw)
        return {"name": name,
                "num_parameter": int(nparam) if nparam.isdigit() else nparam,
                "index_levels": levels,
                "parameters": pnames,
                "usage": _usage(name),
                "source": str(p.relative_to(SRC))}
    return {"name": name, "usage": _usage(name)} if _usage(name) else {}


def _usage(name: str) -> list[str]:
    """The reaction's own error messages, which say what it expects.

    Not every reaction states its parameters in `configure`; many check them
    in the constructor and throw with a sentence naming them. That sentence
    is the only documentation there is, and it is a good one, so it is worth
    reading out rather than leaving people to trigger it.
    """
    cls = None
    reg = re.compile(r'TISSUE_REGISTER_REACTION\(\s*(\w+)\s*,\s*"'
                     + re.escape(name) + r'"')
    for p in _sources():
        text = p.read_text(errors="ignore")
        m = reg.search(text)
        if not m:
            continue
        cls = m.group(1)
        at = text.find(f"class {cls} ")
        if at < 0:
            at = text.find(f"class {cls}:")
        body = text[at:m.start()] if at >= 0 else text
        out = []
        for msg in re.findall(r'runtime_error\(\s*((?:"[^"]*"\s*)+)\)',
                              body):
            joined = "".join(re.findall(r'"([^"]*)"', msg))
            if name.split("::")[-1].lower() in joined.lower() or \
                    "parameter" in joined.lower() or "index" in joined.lower():
                out.append(" ".join(joined.split()))
        return out
    return []
