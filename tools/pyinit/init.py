from collections import defaultdict
from scipy.spatial.distance import euclidean #type: ignore
from itertools import repeat
from typing import Any, Callable, List, Dict, Tuple
from pyvista import PolyData #type: ignore
from numpy import ndarray


class InitWall:
    wid: int
    c1: int
    c2: int
    p1: int
    p2: int
    vars: List[float]
    
    def __init__(self, wid, c1, c2, p1, p2, vars):
        self.wid = wid
        self.c1 = c1
        self.c2 = c2
        self.p1 = p1
        self.p2 = p2
        self.vars = vars

    def to_init_conn(self) -> str:
        vals = [self.wid, self.c1, self.c2, self.p1, self.p2]
        return " ".join([str(v) for v in vals])


    def to_init_vars(self) -> str:
        return " ".join([str(v) for v in self.vars])


class InitVertex:
    x: float
    y: float
    z: float
    
    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def to_init(self) -> str:
        return " ".join([str(v) for v in [self.x, self.y, self.z]])


class InitCell:
    vars: List[float]
    
    def __init__(self, vars: List[float]):
        self.vars = vars

    def to_init(self) -> str:
        return " ".join([str(v) for v in self.vars])

    
class InitTissue:
    ws: List[InitWall]
    cs: List[InitCell]
    vs: List[InitVertex]
    
    def __init__(self,
                 ws: List[InitWall],
                 cs: List[InitCell],
                 vs: List[InitVertex]):
        self.ws = ws
        self.cs = cs
        self.vs = vs

        
    def to_init(self) -> str:
        ncells = len(self.cs)
        nwalls = len(self.ws)
        npoints = len(self.vs)

        wall_conn = "\n".join([w.to_init_conn() for w in self.ws])
        ndims = str(3)
        point_repr = "\n".join([v.to_init() for v in self.vs])
        nwvars = str(len(self.ws[0].vars))
        wall_vars = "\n".join([w.to_init_vars() for w in self.ws])
        ncvars = len(self.cs[0].vars)
        cell_vars = "\n".join([c.to_init() for c in self.cs])
        
        return f"""
{ncells} {nwalls} {npoints}
{wall_conn}
        
{npoints} {ndims}
{point_repr}

{nwalls} {nwvars} 0
{wall_vars}

{ncells} {ncvars}
{cell_vars}
"""


def reshape_fs(fs: ndarray) -> List[List[int]]:
    rfs: List[List[int]] = []
    while True:
        if len(fs) == 0:
            return rfs

        n_vs = fs[0]
        vs = fs[1:(n_vs+1)]

        rfs.append(vs)
        fs = fs[(n_vs+1):]


def add_region(c, m, i, attr="region", pos=13):
    r = m.cell_data[attr][i]

    c.vars[pos] = r

    return c

    
def ident(c: InitCell, m: PolyData, i: int) -> InitCell:
    return c


def to_init(m: PolyData,
            f: Callable[[InitCell, PolyData, int], InitCell]=ident,
            nvars=39) -> InitTissue:
    vs = m.points
    fs = reshape_fs(m.faces)
    
    w_to_cids = defaultdict(list)
    w_lens = dict()
    
    for i, vinds in enumerate(fs):
        vinds_ = list(vinds) + list([vinds[0]])
        ws_ = list(zip(vinds_, vinds_[1:]))
        for w in ws_:
            w_ = frozenset(w)
            w_to_cids[w_].append(i)
            w_lens[w_] = euclidean(vs[w[0]], vs[w[1]])

    ws: List[InitWall] = list()
    for i, k in enumerate(w_to_cids.keys()):
        cs = w_to_cids[k]
        if len(cs) == 1:
            cs.append(-1)
            
        wall: Tuple[int, ...] = tuple(k)
        ws.append(InitWall(i, cs[0], cs[1], wall[0], wall[1], [w_lens[k]]))

    vs = [InitVertex(v[0], v[1], v[2]) for v in vs]

    icells = [f(InitCell(list(repeat(0.0, nvars))), m, i)
              for i, c in enumerate(fs)]

    return InitTissue(ws, icells, vs)


def print_init(m: PolyData, fn: str) -> None:
    t = to_init(m)

    with open(fn, "w") as fout:
        fout.write(t.to_init())
