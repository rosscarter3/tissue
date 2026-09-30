# Driving tissue from Python

```python
from tissue import TissueInit, Model, simulate, outlines

t = TissueInit()
t.add_cell([(0, 0), (10, 0), (10, 10), (0, 10)])
t.add_cell([(10, 0), (20, 0), (20, 10), (10, 10)])
t.set_cell_variables(["cellType", "A0"], defaults={"A0": "area"})

m = (Model()
     .add("Pressure2D::AreaElastic", [3.0, 30.0, 0.0], [1])
     .add("WallMechanics::Spring", [600.0, 1.0], [0]))

frames = simulate(t, m, end=4.0, prints=4)
```

Start with [`notebooks/tissue_quickstart.ipynb`](notebooks/tissue_quickstart.ipynb),
which builds a tissue, runs it, plots it and then makes it lobe.

## What this is for

The simulator reads three files: an init holding the tissue, a model holding
the reactions, and a solver saying how long to run. Writing them by hand is
the barrier to starting, and not because they are long.

An init states its topology explicitly — every wall names the two cells it
separates and the two vertices it joins, and four blocks have to agree on how
many of everything there is. A wall listed against the wrong cell loads,
runs, and is simply a different tissue from the one you meant. All of it
follows from the cell outlines, so `TissueInit` derives it: vertices shared
between cells are recognised by position, so an interface between two cells
is one wall with a cell on each side rather than two walls back to back.

A model block declares its own parameter and index counts, and what those
should be is not written down outside the C++. `reactions()` lists the 206
the simulator has and `describe(name)` says what each wants, reading the
`configure` call or, for reactions that check their arguments in the
constructor instead, the sentence they throw.

## Installing

The package is pure Python and needs nothing to build an init or a model.
`simulate` needs a built simulator, found under the checkout or named by
`TISSUE_SIM`. Plotting in the notebook needs matplotlib.

```sh
pip install -e python              # or: uv pip install -e python
python -m pytest python/tests
```

## The three traps in the output format

`read_frames` handles these, and they are worth knowing if you ever parse the
output yourself, because each returns plausible numbers when read wrong:

- A cell row starts with its vertex **count**, then that many vertex indices,
  and only then the variables.
- The cell block's declared column count is **two more** than the number of
  variables: the last two fields are volume and wall count.
- The wall block's declared count **excludes** the two vertex columns and
  **includes** the wall index and three geometric values after the variables.

## Related tools

- `tools/pyinit` converts a PyVista-readable mesh into an init, for geometry
  that comes from a file rather than from polygons you have in hand.
- `tools/tissue-init-generator` and `tools/createInit` are the older
  generators for specific starting geometries.
