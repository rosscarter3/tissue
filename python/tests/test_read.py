"""Reading what the simulator prints.

The format has three traps that return plausible numbers when misread, which
is the worst kind, so each one is pinned here.
"""
import pytest

from tissue import read_frames


def frame_text(cell_vars, wall_vars, nvert=4):
    verts = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]][:nvert]
    out = ["1", f"{len(verts)} 2"]
    out += [" ".join(f"{c:.6f}" for c in v) for v in verts]
    out.append(f"1 {len(cell_vars) + 2}")
    out.append(f"{len(verts)} " + " ".join(str(i) for i in range(len(verts)))
               + " " + " ".join(f"{v:.6f}" for v in cell_vars)
               + f" 2.5 {len(verts)}")
    out.append(f"{len(verts)} {len(wall_vars) + 4}")
    for i in range(len(verts)):
        out.append(f"{i} {(i + 1) % len(verts)} "
                   + " ".join(f"{v:.6f}" for v in wall_vars)
                   + f" {i} 1.0 1.0 0.0")
    return "\n".join(out) + "\n"


@pytest.fixture
def frame():
    return read_frames(frame_text([7.0, 8.0, 9.0], [0.5, 1.5]),
                       is_text=True)[0]


def test_cell_variables_are_read_past_the_vertex_list(frame):
    """Trap one: a cell row starts with its vertex count and indices."""
    assert frame.cell_vars[0] == [7.0, 8.0, 9.0]
    assert frame.cell_vertices[0] == [0, 1, 2, 3]


def test_the_declared_cell_width_counts_two_trailing_fields(frame):
    """Trap two: volume and wall count are not variables."""
    assert len(frame.cell_vars[0]) == 3
    assert frame.cell_volume[0] == pytest.approx(2.5)


def test_wall_variables_are_read_past_the_vertex_columns(frame):
    """Trap three: the declared wall width excludes v1 and v2, and includes
    the wall index and three geometric values after the variables."""
    assert frame.wall_vars[0] == [0.5, 1.5]
    assert frame.walls[0] == (0, 1)


def test_the_geometry_is_kept_out_of_the_variables(frame):
    """Otherwise wall_vars silently gains four columns that are not variables,
    and anything indexing by variable number reads the wrong one."""
    assert len(frame.wall_vars[0]) == 2
    assert frame.wall_length[0] == pytest.approx(1.0)
    assert frame.wall_strain[0] == pytest.approx(0.0)


def test_wall_variable_zero_is_the_resting_length(frame):
    assert frame.wall_vars[0][0] == 0.5


def test_a_truncated_frame_is_dropped_rather_than_half_read():
    """A killed run leaves a partial frame, which must not be returned."""
    text = frame_text([1.0], [1.0]) + "4 2\n0.0 0.0\n1.0"
    assert len(read_frames(text, is_text=True)) == 1


def test_the_interior_helper_drops_the_boundary_ring():
    two = read_frames(frame_text([0.0], [1.0]), is_text=True)[0]
    two.cell_vars.append([1.0])
    assert two.interior() == [0]


def test_variables_can_be_read_across_cells(frame):
    assert frame.var(1) == [8.0]


def test_the_dimension_comes_from_the_frame(frame):
    assert frame.dim == 2
