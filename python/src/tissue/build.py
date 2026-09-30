"""Build a tissue init file from polygons, without writing topology by hand.

A tissue init file states its topology explicitly: every wall names the two
cells it separates and the two vertices it joins, cells are implied by the
walls that reference them, and four separate blocks have to agree on how many
of everything there is. Writing that by hand is where the mistakes happen,
and they are quiet ones -- a wall listed against the wrong cell still loads
and still runs.

All of it follows from the cell outlines, so it is derived here instead:

    t = TissueInit()
    t.add_cell([(0, 0), (10, 0), (10, 10), (0, 10)])
    t.add_cell([(10, 0), (20, 0), (20, 10), (10, 10)])
    t.write("two_cells.init")

Vertices shared between cells are recognised by position, so the wall down
the middle is one wall with a cell on each side rather than two walls back to
back. That is not a detail: a shared wall is what makes one cell's lobe its
neighbour's indentation.
"""
from __future__ import annotations

import math
from collections import OrderedDict


class TissueInit:
    """Cells in, init file out.

    Parameters
    ----------
    dimension:
        2 for a flat sheet, 3 for a surface in space.
    tol:
        How close two points must be to count as the same vertex, in the same
        units as the coordinates. Too small and a shared wall silently
        becomes two, which removes the mechanical coupling between the cells
        without removing anything visible.
    """

    def __init__(self, dimension: int = 2, tol: float = 1e-6):
        if dimension not in (2, 3):
            raise ValueError("dimension must be 2 or 3")
        self.dimension = dimension
        self.tol = tol
        self._verts: OrderedDict[tuple, int] = OrderedDict()
        self._pos: list[tuple[float, ...]] = []
        self._cells: list[list[int]] = []
        self.wall_variables: list[str] = []
        self.cell_variables: list[str] = []
        self._wall_defaults: dict[str, float] = {}
        self._cell_defaults: dict[str, float] = {}

    # ---------------------------------------------------------------- cells

    def _vertex(self, p) -> int:
        p = tuple(float(x) for x in p)
        if len(p) != self.dimension:
            raise ValueError(f"point {p} is not {self.dimension}-dimensional")
        key = tuple(round(x / self.tol) for x in p)
        if key not in self._verts:
            self._verts[key] = len(self._pos)
            self._pos.append(p)
        return self._verts[key]

    def add_cell(self, outline) -> int:
        """Add one cell from its outline, in order around the boundary.

        The outline is treated as closed, so the last point joins the first;
        do not repeat it.
        """
        pts = [tuple(float(x) for x in p) for p in outline]
        if len(pts) > 1 and _close(pts[0], pts[-1], self.tol):
            pts = pts[:-1]
        if len(pts) < 3:
            raise ValueError("a cell needs at least three distinct points")
        ids = [self._vertex(p) for p in pts]
        if len(set(ids)) != len(ids):
            raise ValueError("outline visits the same vertex twice")
        self._cells.append(ids)
        return len(self._cells) - 1

    def subdivide(self, n: int) -> "TissueInit":
        """Split every cell edge into `n` segments, returning a new tissue.

        This is the discretisation choice that decides whether the model can
        lobe at all. With one segment per cell-cell interface the interface
        is a straight line and stays one however it is pulled; with twenty,
        its interior vertices are free to move and it can undulate, while the
        original corners stay as genuine three-way junctions.
        """
        if n < 1:
            raise ValueError("n must be at least 1")
        out = TissueInit(self.dimension, self.tol)
        out.wall_variables = list(self.wall_variables)
        out.cell_variables = list(self.cell_variables)
        out._wall_defaults = dict(self._wall_defaults)
        out._cell_defaults = dict(self._cell_defaults)
        for cell in self._cells:
            pts = []
            for k, a in enumerate(cell):
                b = cell[(k + 1) % len(cell)]
                pa, pb = self._pos[a], self._pos[b]
                for j in range(n):
                    t = j / n
                    pts.append(tuple(pa[d] + t * (pb[d] - pa[d])
                                     for d in range(self.dimension)))
            out.add_cell(pts)
        return out

    # ------------------------------------------------------------ variables

    def set_wall_variables(self, names, defaults=None) -> None:
        """Name the wall variables, in order from index 1.

        Index 0 is always the resting length and is written for you, so the
        first name here is variable 1. Getting that off by one is how a
        diagnostic came to overwrite the variable next to the one it meant.
        """
        self.wall_variables = list(names)
        self._wall_defaults = dict(defaults or {})

    def set_cell_variables(self, names, defaults=None) -> None:
        self.cell_variables = list(names)
        self._cell_defaults = dict(defaults or {})

    # ----------------------------------------------------------- derivation

    def walls(self) -> list[tuple[int, int, int, int]]:
        """[(v1, v2, cell1, cell2)] with -1 for background, in index order."""
        seen: OrderedDict[frozenset, list] = OrderedDict()
        for ci, cell in enumerate(self._cells):
            for k, a in enumerate(cell):
                b = cell[(k + 1) % len(cell)]
                if a == b:
                    raise ValueError(
                        f"cell {ci} has an edge from vertex {a} to itself; "
                        "a zero-length wall has no direction and the "
                        "mechanics divide by its length")
                key = frozenset((a, b))
                if key not in seen:
                    seen[key] = [a, b, ci, -1]
                elif seen[key][3] == -1:
                    seen[key][3] = ci
                else:
                    raise ValueError(
                        f"edge {sorted(key)} is shared by more than two "
                        "cells; a wall separates exactly two")
        return [tuple(v) for v in seen.values()]

    def check(self) -> None:
        """Raise if the tissue is not one the simulator will accept.

        Deriving the walls from the outlines already rules out most of what
        can go wrong, so what is left is the geometry: a wall shorter than
        the tolerance is two vertices the builder should have merged, and
        the mechanics divide by wall length.
        """
        if not self._cells:
            raise ValueError("no cells")
        walls = self.walls()
        for v1, v2, _, _ in walls:
            d = _dist(self._pos[v1], self._pos[v2])
            if d <= self.tol:
                raise ValueError(
                    f"wall between vertices {v1} and {v2} is {d:g} long, "
                    f"at or below the {self.tol:g} merge tolerance")
        per_cell = {}
        for v1, v2, c1, c2 in walls:
            for c in (c1, c2):
                if c >= 0:
                    per_cell[c] = per_cell.get(c, 0) + 1
        for ci, cell in enumerate(self._cells):
            if per_cell.get(ci, 0) != len(cell):
                raise ValueError(
                    f"cell {ci} has {len(cell)} edges but {per_cell.get(ci, 0)}"
                    " walls; its outline is not closed")

    def area(self, ci: int) -> float:
        """Polygon area of one cell, for seeding a resting area."""
        pts = [self._pos[i] for i in self._cells[ci]]
        s = 0.0
        for k, p in enumerate(pts):
            q = pts[(k + 1) % len(pts)]
            s += p[0] * q[1] - q[0] * p[1]
        return abs(s) / 2.0

    # --------------------------------------------------------------- output

    def to_text(self) -> str:
        self.check()
        walls = self.walls()
        out = [f"{len(self._cells)} {len(walls)} {len(self._pos)}"]
        for w, (v1, v2, c1, c2) in enumerate(walls):
            out.append(f"{w} {c1} {c2} {v1} {v2}")
        out.append("")
        out.append(f"{len(self._pos)} {self.dimension}")
        for p in self._pos:
            out.append(" ".join(f"{x:.6f}" for x in p))
        out.append("")
        out.append(f"{len(walls)} 1 {len(self.wall_variables)}")
        for v1, v2, _, _ in walls:
            length = _dist(self._pos[v1], self._pos[v2])
            vals = [length] + [self._wall_defaults.get(n, 0.0)
                               for n in self.wall_variables]
            out.append(" ".join(f"{x:.6f}" for x in vals))
        out.append("")
        out.append(f"{len(self._cells)} {len(self.cell_variables)}")
        for ci in range(len(self._cells)):
            vals = []
            for n in self.cell_variables:
                d = self._cell_defaults.get(n, 0.0)
                vals.append(self.area(ci) if d == "area" else float(d))
            out.append(" ".join(f"{x:.6f}" for x in vals))
        return "\n".join(out) + "\n"

    def write(self, path) -> str:
        text = self.to_text()
        with open(path, "w") as f:
            f.write(text)
        return str(path)

    def __repr__(self) -> str:
        try:
            nw = len(self.walls())
        except ValueError:
            nw = -1
        return (f"TissueInit({len(self._cells)} cells, {nw} walls, "
                f"{len(self._pos)} vertices, {self.dimension}D)")


def _dist(a, b) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _close(a, b, tol) -> bool:
    return _dist(a, b) <= tol


def hex_grid(rows: int, cols: int, radius: float = 10.0) -> TissueInit:
    """A honeycomb: every interior junction joins three cells.

    That topology is why a vertex model of an epidermis lobes at all. Every
    cell-cell interface is shared by exactly two cells, so one cell's lobe is
    the other's indentation by construction, with nothing passing between
    them.
    """
    t = TissueInit(2)
    dx = radius * math.sqrt(3)
    dy = radius * 1.5
    for r in range(rows):
        for c in range(cols):
            cx = c * dx + (dx / 2 if r % 2 else 0.0)
            cy = r * dy
            t.add_cell([(cx + radius * math.sin(a), cy + radius * math.cos(a))
                        for a in [math.pi / 3 * k for k in range(6)]])
    return t
