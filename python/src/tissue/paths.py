"""Finding the simulator and the source tree it was built from.

Both are looked for relative to this package, so an in-tree checkout works
with no configuration, and both can be overridden for a build somewhere else.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

#: python/src/tissue/paths.py -> tissue -> src -> python -> repository root
ROOT = Path(__file__).resolve().parents[3]

#: The C++ source, used to list what reactions exist.
SRC = Path(os.environ.get("TISSUE_SRC", ROOT / "src" / "tissue"))


def simulator() -> Path:
    """The simulator binary.

    TISSUE_SIM wins; otherwise the usual build directories in the checkout,
    then anything named `simulator` on PATH.
    """
    env = os.environ.get("TISSUE_SIM")
    if env:
        return Path(env)
    for rel in ("build/simulator", "build/bin/simulator",
                "bin/simulator", "build/Release/simulator"):
        p = ROOT / rel
        if p.is_file() and os.access(p, os.X_OK):
            return p
    found = shutil.which("simulator")
    if found:
        return Path(found)
    raise FileNotFoundError(
        "cannot find the tissue simulator. Build it, or set TISSUE_SIM to "
        f"the binary. Looked under {ROOT}.")
