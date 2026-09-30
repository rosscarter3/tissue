"""Build, run and read back in one call.

    frames = simulate(tissue, model, end=24.0, prints=12)

The files still exist -- they are what the simulator reads -- but they go in
a temporary directory unless somewhere is named, so trying something out does
not mean choosing filenames for it.
"""
from __future__ import annotations

import os
import subprocess
import tempfile

from .paths import simulator
from .read import read_frames


def write_all(tissue, model, end: float, prints: int = 10,
              directory=None, name: str = "run", solver: str = "rk5",
              max_step: float = 0.01, eps: float = 1e-3) -> dict:
    """Write the three files the simulator needs, and say where they are."""
    d = directory or tempfile.mkdtemp(prefix="pavesim-")
    os.makedirs(d, exist_ok=True)
    paths_out = {"dir": d,
                 "init": os.path.join(d, f"{name}.init"),
                 "model": os.path.join(d, f"{name}.model"),
                 "solver": os.path.join(d, f"{name}.solver")}
    tissue.write(paths_out["init"])
    model.write(paths_out["model"])
    with open(paths_out["solver"], "w") as f:
        if solver == "qs":
            f.write(f"QuasiStatic\n0 {end}\n0 {prints}\n0.25 1e-3 20000\n")
        else:
            f.write(f"RK5Adaptive\n0 {end}\n0 {prints}\n{max_step} {eps}\n")
    return paths_out


def simulate(tissue, model, end: float, prints: int = 10, directory=None,
             name: str = "run", timeout: float | None = None, **kw):
    """Run and return the frames.

    Raises with the simulator's own last line on failure, which is where it
    says what it did not like about the model file -- usually the count of
    parameters or of indices, and usually in a sentence naming them.
    """
    f = write_all(tissue, model, end, prints, directory, name, **kw)
    p = subprocess.run([str(simulator()), f["model"], f["init"],
                        f["solver"]], cwd=f["dir"], capture_output=True,
                       text=True, timeout=timeout)
    if p.returncode != 0:
        tail = (p.stderr or "").strip().splitlines()
        raise RuntimeError(tail[-1] if tail else
                           f"simulator exited {p.returncode}")
    return read_frames(p.stdout, is_text=True)


def outlines(frame, interior_only: bool = True):
    """[(cell index, Nx2 outline)] for plotting.

    `interior_only` drops the boundary ring, whose shape is set by the edge
    of the patch rather than by the model.
    """
    import numpy as np
    out = []
    for i, (ids, cvars) in enumerate(zip(frame.cell_vertices,
                                         frame.cell_vars)):
        if interior_only and cvars and cvars[0] != 0:
            continue
        out.append((i, np.array([frame.vertices[j][:2] for j in ids], float)))
    return out
