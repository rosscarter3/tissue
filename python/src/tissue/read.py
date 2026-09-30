"""Reading what the simulator prints.

The print-flag-0 format is simple and has three traps that are easy to get
wrong and give plausible-looking numbers when you do. They are handled here
so they are handled once:

- **A cell row does not start with its variables.** It starts with the vertex
  count, then that many vertex indices, and only then the variables. So cell
  variable `v` is at field `1 + numVertex + v`.
- **The cell block's declared column count is two more than the number of
  variables**: the last two fields of every row are the cell's volume and its
  wall count.
- **The wall block's declared count excludes the two vertex columns** and
  includes the wall index and three geometric values after the variables.

Reading any of them wrong shifts every variable by a fixed offset, and a run
read that way still produces numbers.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Frame:
    """One printed state."""
    vertices: list[list[float]]            # nVertex x dim
    cell_vertices: list[list[int]]         # per cell, its vertex indices
    cell_vars: list[list[float]]           # per cell, its variables
    cell_volume: list[float]
    walls: list[tuple[int, int]]           # (v1, v2)
    wall_vars: list[list[float]]           # index 0 is the resting length
    wall_length: list[float]               # current length, as printed
    wall_strain: list[float]

    @property
    def dim(self) -> int:
        return len(self.vertices[0]) if self.vertices else 0

    def var(self, index: int) -> list[float]:
        """Cell variable `index` across every cell."""
        return [c[index] for c in self.cell_vars]

    def interior(self, type_index: int = 0) -> list[int]:
        """Cells whose type variable is 0.

        The boundary ring's shape is set by the edge of the patch rather than
        by the model, so most measurements want to leave it out.
        """
        return [i for i, v in enumerate(self.cell_vars)
                if v and v[type_index] == 0]


def read_frames(path_or_text: str, is_text: bool = False) -> list[Frame]:
    """Parse a print-flag-0 file, or the simulator's stdout, into frames."""
    text = (path_or_text if is_text
            else open(path_or_text, errors="ignore").read())
    tok = text.split()
    pos = 1 if tok else 0

    def take(n: int) -> list[str]:
        nonlocal pos
        out = tok[pos:pos + n]
        pos += n
        return out

    frames: list[Frame] = []
    while pos + 1 < len(tok):
        try:
            nv, dim = int(tok[pos]), int(tok[pos + 1])
        except ValueError:
            break
        pos += 2
        raw = take(nv * dim)
        if len(raw) < nv * dim:
            break
        verts = [[float(raw[i * dim + d]) for d in range(dim)]
                 for i in range(nv)]

        ncell, ncols = int(take(1)[0]), int(take(1)[0])
        nvars = ncols - 2                      # the last two are volume, nWall
        cverts, cvars, cvol = [], [], []
        for _ in range(ncell):
            n = int(take(1)[0])
            cverts.append([int(x) for x in take(n)])
            cvars.append([float(x) for x in take(nvars)])
            vol, _nwall = take(2)
            cvol.append(float(vol))

        # The declared width covers the wall variables, the wall index and
        # three geometric values, but not the two vertex columns. Splitting
        # them here keeps `wall_vars` meaning what its name says: a
        # hard-coded width silently desynchronised the whole parse when the
        # model gained a variable.
        nwall, wcols = int(take(1)[0]), int(take(1)[0])
        nwallvar = wcols - 4
        walls, wvars, wlen, wstrain = [], [], [], []
        for _ in range(nwall):
            v1, v2 = int(take(1)[0]), int(take(1)[0])
            walls.append((v1, v2))
            wvars.append([float(x) for x in take(nwallvar)])
            take(1)                                   # the wall's own index
            d, _dl, strain = (float(x) for x in take(3))
            wlen.append(d)
            wstrain.append(strain)

        frames.append(Frame(verts, cverts, cvars, cvol, walls, wvars,
                            wlen, wstrain))
    return frames
