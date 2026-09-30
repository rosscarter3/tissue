"""Building the files the simulator reads.

The point of these helpers is to remove the two ways a file can be accepted
and still be wrong: a wall listed against the wrong cell, and a reaction
given the wrong number of parameters or indices. So the tests are mostly
about the topology being derived correctly rather than about formatting.
"""
import pytest

from tissue import Model, TissueInit, describe, hex_grid, reactions

SQ = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]


def two_squares():
    t = TissueInit()
    t.add_cell(SQ)
    t.add_cell([(10.0, 0.0), (20.0, 0.0), (20.0, 10.0), (10.0, 10.0)])
    return t


def test_a_single_cell_has_every_wall_on_the_background():
    t = TissueInit()
    t.add_cell(SQ)
    walls = t.walls()
    assert len(walls) == 4
    assert all(c2 == -1 for _, _, _, c2 in walls)


def test_two_cells_sharing_an_edge_share_one_wall():
    """The whole reason a vertex model of an epidermis lobes.

    One wall with a cell on each side makes one cell's lobe the other's
    indentation by construction. Two walls back to back would leave the
    cells mechanically independent and nothing would interdigitate.
    """
    walls = two_squares().walls()
    assert len(walls) == 7                       # 4 + 4 - 1 shared
    shared = [w for w in walls if w[3] != -1]
    assert len(shared) == 1
    assert set(shared[0][2:]) == {0, 1}


def test_vertices_are_shared_by_position():
    t = two_squares()
    assert len(t._pos) == 6                      # not 8


def test_a_repeated_closing_point_is_ignored():
    t = TissueInit()
    t.add_cell(SQ + [SQ[0]])
    assert len(t.walls()) == 4


def test_an_outline_that_revisits_a_vertex_is_refused():
    with pytest.raises(ValueError):
        TissueInit().add_cell([(0, 0), (1, 0), (0, 0), (1, 1)])


def test_three_cells_on_one_edge_is_refused():
    """A wall separates exactly two cells; more is not representable."""
    t = TissueInit()
    for dy in (0.0, 10.0, 20.0):
        t.add_cell([(0, dy), (10, dy), (10, dy + 10), (0, dy + 10)])
    t.walls()                                    # stacked, fine
    bad = TissueInit()
    edge = [(0.0, 0.0), (10.0, 0.0)]
    for k in range(3):
        bad.add_cell(edge + [(10.0, 10.0 + k), (0.0, 10.0 + k)])
    with pytest.raises(ValueError, match="more than two"):
        bad.walls()


def test_subdivide_keeps_the_cells_and_multiplies_the_walls():
    t = two_squares()
    s = t.subdivide(5)
    assert len(s._cells) == len(t._cells)
    assert len(s.walls()) == 5 * len(t.walls())


def test_subdivide_keeps_the_shared_wall_shared():
    """The interface must stay one interface after refinement.

    If refinement split the shared edge into two independent chains the
    cells would come apart, and the model would look finer while having lost
    the thing that makes it a tissue.
    """
    s = two_squares().subdivide(4)
    shared = [w for w in s.walls() if w[3] != -1]
    assert len(shared) == 4


def test_area_is_the_polygon_area():
    t = TissueInit()
    t.add_cell(SQ)
    assert t.area(0) == pytest.approx(100.0)


def test_cell_variables_can_be_seeded_from_the_area():
    t = TissueInit()
    t.add_cell(SQ)
    t.set_cell_variables(["cellType", "A0"], defaults={"A0": "area"})
    assert t.to_text().strip().splitlines()[-1].split()[1] == "100.000000"


def test_the_written_header_counts_match_the_blocks():
    t = two_squares()
    t.set_wall_variables(["m", "force"])
    t.set_cell_variables(["cellType"])
    lines = [l for l in t.to_text().splitlines() if l.strip()]
    ncell, nwall, nvert = (int(x) for x in lines[0].split())
    assert (ncell, nwall, nvert) == (2, 7, 6)
    # wall block header: numWall 1 numWallVar, and rows are 1 + numWallVar wide
    hdr = next(i for i, l in enumerate(lines) if l.split() == ["7", "1", "2"])
    assert len(lines[hdr + 1].split()) == 3


def test_resting_length_is_the_first_wall_column():
    """Index 0 is the length; named variables start at 1.

    Being off by one here is how a diagnostic came to overwrite the variable
    next to the one it meant to write.
    """
    t = TissueInit()
    t.add_cell(SQ)
    t.set_wall_variables(["m"])
    lines = [l for l in t.to_text().splitlines() if l.strip()]
    hdr = next(i for i, l in enumerate(lines) if l.split() == ["4", "1", "1"])
    assert float(lines[hdr + 1].split()[0]) == pytest.approx(10.0)


def test_a_degenerate_edge_is_refused():
    """A zero-length wall has no direction and the mechanics divide by it."""
    t = TissueInit()
    t.add_cell(SQ)
    t._cells[0] = [t._cells[0][0]] + t._cells[0]
    with pytest.raises(ValueError, match="itself"):
        t.walls()


def test_points_closer_than_the_tolerance_become_one_vertex():
    """Which then shows up as a degenerate outline rather than a tiny wall.

    Merging is what makes a shared edge shared, so it has to happen; the
    consequence is that an outline with two near-coincident corners is
    refused at the point it is added, where the message can still name the
    cell, rather than becoming a zero-length wall later.
    """
    t = TissueInit(tol=1e-6)
    with pytest.raises(ValueError, match="same vertex twice"):
        t.add_cell([(0.0, 0.0), (10.0, 0.0), (10.0, 1e-9), (0.0, 10.0)])


def test_the_honeycomb_shares_walls_between_neighbours():
    g = hex_grid(2, 2)
    assert len(g._cells) == 4
    assert sum(1 for w in g.walls() if w[3] != -1) > 0


def test_a_model_block_states_its_own_counts():
    m = Model().add("Pressure2D::AreaElastic", [3.0, 30.0, 0.0134], [[2]])
    lines = m.to_text().splitlines()
    assert lines[0].split() == ["1", "0", "0"]
    head = next(l for l in lines if l.startswith("Pressure2D"))
    assert head.split()[1:] == ["3", "1", "1"]


def test_a_flat_index_list_is_taken_as_one_level():
    a = Model().add("X", [1.0], [0, 1]).to_text()
    b = Model().add("X", [1.0], [[0, 1]]).to_text()
    assert a == b


def test_reactions_are_discovered_from_the_simulator_source():
    names = reactions()
    assert len(names) > 100
    assert "Pressure2D::AreaElastic" in names
    assert reactions("AreaElastic") == ["Pressure2D::AreaElastic"]


def test_a_reaction_reports_its_parameter_names():
    d = describe("Pressure2D::AreaElastic")
    assert d["num_parameter"] == 3
    assert d["parameters"] == ["P_turgor", "K_area", "k_areaGrowth"]
    assert d["index_levels"] == [1]


def test_a_reaction_that_checks_in_its_constructor_still_reports_usage():
    """Not every reaction names its parameters in configure.

    Many throw a sentence instead, and that sentence is the only
    documentation there is.
    """
    d = describe("WallMechanics::Spring")
    assert any("K_force" in u for u in d["usage"])
