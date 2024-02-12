import pyvista as pv  # type: ignore
from pyvista import UnstructuredGrid  # type: ignore
import os  # type: ignore
from plyfile import PlyData, PlyElement #type: ignore
from typing import List, Dict, Any, Callable, Tuple
import pandas as pd #type: ignore
import flower.eply as eply


class Vertex3:
    x: float
    y: float
    z: float

    def __init__(self, x: float, y: float, z: float):
        self.x = x
        self.y = y
        self.z = z

    def to_dat(self):
        return " ".join([str(self.x),
                         str(self.y),
                         str(self.z)])


class CCell:
    ctr: Vertex3
    varss: Dict[str, Any]

    def __init__(self, ctr: Vertex3, varss: Dict[str, Any]):
        self.ctr = ctr
        self.varss = varss

    def _rm_brackets(self, s):
        return s.replace("[", "").replace("]", "")

    def _varss_to_dat(self) -> str:
        return " ".join([self._rm_brackets(str(v))
                         for _, v in self.varss.items()])

    def to_dat(self, i: int) -> str:
        dvol = str(10.0)
        return " ".join([self.ctr.to_dat(),
                         dvol,
                         str(i),
                         self._varss_to_dat()])


def cells_to_dat(ccells: List[CCell]) -> str:
    ntpoints = 1
    ncells = len(ccells)
    nvars = len(ccells[0].varss.keys()) + 7
    tPointIndex = 0

    header = " ".join([str(ntpoints),
                       str(ncells),
                       str(nvars),
                       str(tPointIndex)])

    return "\n".join([header] + [c.to_dat(i) for i, c in enumerate(ccells)])


class TissueOutput:
    tcells: List[UnstructuredGrid]
    twalls: List[UnstructuredGrid]

    def __init__(self, tcells, twalls):
        self.tcells = tcells
        self.twalls = twalls

    @classmethod
    def from_dir(cls, out_dir: str, everyN=1):
        ind_fnames_cells = list()
        ind_fnames_walls = list()
        print(out_dir)
        for f in os.listdir(out_dir):
            if 'cells' in f:
                i = cls._extract_ind(f, 'cells')
                ind_fnames_cells.append((i, f))
            elif 'walls' in f:
                i = cls._extract_ind(f, 'walls')
                ind_fnames_walls.append((i, f))

        tcells = [pv.UnstructuredGrid(os.path.join(out_dir, fn))
                  for i, fn in sorted(ind_fnames_cells,
                                      key=lambda x: x[0])[::everyN]]

        twalls = [pv.UnstructuredGrid(os.path.join(out_dir, fn))
                  for i, fn in sorted(ind_fnames_walls,
                                      key=lambda x: x[0])[::everyN]]

        return cls(tcells, twalls)

    @classmethod
    def from_dir_cells(cls, out_dir: str, everyN=1):
        ind_fnames_cells = list()
        print(out_dir)
        for f in os.listdir(out_dir):
            if 'cells' in f:
                i = cls._extract_ind(f, 'cells_')
                ind_fnames_cells.append((i, f))

        tcells = [pv.read(os.path.join(out_dir, fn))
                  for i, fn in sorted(ind_fnames_cells,
                                      key=lambda x: x[0])[::everyN]]

        return cls(tcells, list())

    @classmethod
    def from_dir_ply(cls, out_dir: str, everyN: int=1):
        def _extract_ind_ply(fn: str) -> int:
            f, _ = os.path.splitext(fn)
            return int(f.split("_")[1])
        
        ind_fnames: List[Tuple[str, str]] = list()
        print(out_dir)
        for f in os.listdir(out_dir):
            if 'output' in f:
                i = _extract_ind_ply(fn=f)
                ind_fnames.append((i, f))

        tcells = []
        for i, fn in sorted(ind_fnames, key=lambda x: x[0])[::everyN]:
            print(i)
            tcells.append(eply.get_mesh_vars(os.path.join(out_dir, fn)))

        return cls(tcells, None)
    
    @staticmethod
    def _extract_ind(fn: str, word: str) -> int:
        ind_fn = fn.split(word)[1]
        ind_s, _ = os.path.splitext(ind_fn)

        return int(ind_s)

    def filter_with_key(self,
                        f: Callable[[int], bool]):
        tcells_ = list()
        twalls_ = list()
        for i, (tcells, twalls) in enumerate(zip(self.tcells, self.twalls)):
            if f(i):
                tcells_.append(tcells)
                twalls_.append(twalls)

        return TissueOutput(tcells_, twalls_)

    def _cell_centres_ugrid(self, ugrid: UnstructuredGrid) -> List[CCell]:
        ccells = list()
        ugrid = ugrid.cell_centers()

        for i in range(ugrid.n_cells):
            cc = Vertex3(ugrid.points[i][0],
                         ugrid.points[i][1],
                         ugrid.points[i][2])
            vnames = set(ugrid.array_names)
            cvars = {var_name: ugrid.cell_arrays[var_name][i]
                     for var_name in vnames}
            ccells.append(CCell(cc, cvars))

        return ccells

    def cell_centres(self) -> List[List[CCell]]:
        ccs = list()
        for i, ugrid in enumerate(self.tcells):
            print(i)
            ccs.append(self._cell_centres_ugrid(ugrid))

        return ccs

    def save_cells(self, tag="cells-"):
        tag_ = tag + "{i}.vtu"
        for i, tc in enumerate(self.tcells):
            tc.save(tag_.format(i=i))
        
    def fmap_cells(self, f: Callable[[UnstructuredGrid], UnstructuredGrid]):
        self.tcells = [f(tc) for tc in self.tcells]

    def for_cells(self, f: Callable[[UnstructuredGrid], None]):
        for tc in self.tcells:
            f(tc)

    def fmap_walls(self, f: Callable[[UnstructuredGrid], UnstructuredGrid]):
        self.twalls = [f(tc) for tc in self.twalls]

    def to_pandas_cells(self):
        ds = list()
        for tc in self.tcells:
            var_names = tc.array_names
            ds.append(pd.DataFrame({vnm: tc.cell_data[vnm]
                                    for vnm in var_names
                                    if len(tc.cell_data[vnm].shape) == 1}))

        return ds

    def save_cells_dir(self, dir: str):
        from os.path import join
        
        for i, tc in enumerate(self.tcells):
            tc.save(join(dir, f"cells{i}.vtk"))

        return


class SolverOutputStep():
    step: int
    time: float
    ncells: int
    nwalls: int
    nvertex: int
    n_succ_steps: int
    n_unsucc_steps: int
    t_diff: float
    n_n_succ_steps: int
    n_n_unsucc_steps: int

    def __init__(self, step, time, ncells, nwalls, nvertex,
                 n_succ_steps, n_unsucc_steps, t_diff,
                 n_n_succ_steps, n_n_unsucc_steps):
        self.step = step
        self.time = time
        self.ncells =  ncells
        self.nwalls = nwalls
        self.nvertex = nvertex
        self.n_succ_steps = n_succ_steps
        self.n_unsucc_steps = n_unsucc_steps
        self.t_diff = t_diff
        self.n_n_succ_steps = n_n_succ_steps
        self.n_n_unsucc_steps = n_n_unsucc_steps

    @classmethod
    def from_stdout(cls, s):
        cols = s.strip().split()
        return cls(int(cols[0]),
                   float(cols[1]),
                   int(cols[2]),
                   int(cols[3]),
                   int(cols[4]),
                   int(cols[5]),
                   int(cols[6]),
                   float(cols[7]),
                   int(cols[9]),
                   int(cols[10]))

    def __repr__(self):
        return " ".join([str(self.step),
                         str(self.time),
                         str(self.ncells),
                         str(self.nwalls),
                         str(self.nvertex),
                         str(self.n_succ_steps),
                         str(self.n_unsucc_steps),
                         str(self.t_diff),
                         str(self.n_n_succ_steps),
                         str(self.n_n_unsucc_steps)])


class SolverOutput():
    steps: List[SolverOutputStep]
    
    def __init__(self, steps):
        self.steps = steps

    @classmethod
    def from_stdout(cls, solver_stdout):
        def _is_solver_step(ln):
            cols = ln.strip().split()
            return cols[0].isdigit()
        
        steps = [SolverOutputStep.from_stdout(ln)
                 for ln in solver_stdout if _is_solver_step(ln)]

        return cls(steps)
