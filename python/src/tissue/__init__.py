"""Build, run and read tissue simulations from Python.

    from tissue import TissueInit, Model, simulate, outlines

    t = TissueInit()
    t.add_cell([(0, 0), (10, 0), (10, 10), (0, 10)])
    t.add_cell([(10, 0), (20, 0), (20, 10), (10, 10)])
    t.set_cell_variables(["cellType", "A0"], defaults={"A0": "area"})

    m = (Model()
         .add("Pressure2D::AreaElastic", [3.0, 30.0, 0.0], [1])
         .add("WallMechanics::Spring", [600.0, 1.0], [0]))

    frames = simulate(t, m, end=4.0, prints=4)

The simulator reads three files: an init holding the tissue, a model holding
the reactions, and a solver saying how long to run. Writing them by hand is
the barrier, and not because they are long -- the init states its topology
explicitly, so a wall listed against the wrong cell loads, runs, and is
simply a different tissue from the one intended.

`TissueInit` derives the whole topology from the cell outlines, `Model`
fills in the counts each reaction block declares, and `reactions()` and
`describe()` say what the simulator actually has and what each one wants.
"""
from .build import TissueInit, hex_grid
from .paths import ROOT, SRC, simulator
from .read import Frame, read_frames
from .reactions import Model, describe, reactions
from .run import outlines, simulate, write_all

__all__ = ["TissueInit", "hex_grid", "Model", "reactions", "describe",
           "simulate", "write_all", "outlines", "read_frames", "Frame",
           "simulator", "ROOT", "SRC"]
