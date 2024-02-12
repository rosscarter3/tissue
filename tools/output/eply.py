import pyvista as pv  # type: ignore
from pyvista import UnstructuredGrid  # type: ignore
import os  # type: ignore
from plyfile import PlyData, PlyElement #type: ignore
import numpy as np
from collections import defaultdict
import sys
import vtk

def get_obb(m):
    corner = np.array([0.]*3)
    max = np.array([0.]*3)
    mid = np.array([0.]*3)
    min = np.array([0.]*3)
    size = np.array([0.]*3)
    obb_tree = vtk.vtkOBBTree()
    obb_tree.ComputeOBB(
        m,
        corner,
        max,
        mid,
        min,
        size)

    return (corner, max, mid, min)

def get_mesh_vars(fn):
    with open(fn, "rb") as fin:
        plydata = PlyData.read(fin)

    
    vs = np.array([plydata['vertex'].data['x'],
                   plydata['vertex'].data['y'],
                   plydata['vertex'].data['z']]).T

    fdata = plydata['face'].data['vertex_index']

    fs = []
    rem_ids = []
    for i, vinds in enumerate(fdata):
        if len(vinds) > 2:
            fs.append([len(vinds)] + list(vinds))
        else:
            rem_ids.append(i)
            
    fs = np.hstack(fs)
    rem_ids = set(rem_ids)

    m = pv.PolyData(vs, fs)

    vdata = plydata['face'].data['var']
    vars = np.stack(np.array([v for i, v in enumerate(vdata)
                              if i not in rem_ids]))
    ncols = vars.shape[1]

    base_nm = "cell variable {i}"

    for i in range(ncols):
        m.cell_data[base_nm.format(i=str(i))] = vars[:, i]

    return m


def get_faces(m):
    face_vs = {}
    fs = list(m.faces)
    i=0
    while fs:
        n = fs[0]
        vs = fs[1:(n+1)]
        face_vs[i] = vs
        i += 1
        fs = fs[(n+1):]

    return face_vs


def get_edges(vs):
    es = list(zip(vs, vs[1:]))
    es.append((vs[-1], vs[0]))

    return es


def get_neighbours(m):
    neighs = []
    face_vs = get_faces(m)
    d = defaultdict(list)
    for fid in face_vs:
        es = get_edges(face_vs[fid])
        for e in es:
            d[frozenset(e)].append(fid)

    for e, cs in d.items():
        if len(cs) == 2:
            neighs.append(tuple(cs))

    return neighs


def v_length(ax):
    return np.sqrt(ax[0]**2 + ax[1]**2)


def aspect_ratio(ax1, ax2):
    max_length = v_length(ax1)
    min_length = v_length(ax2)
    return (max_length / min_length,
            max_length,
            min_length)
            


def write_neighbours_info(m):
    import networkx as nx
    G = nx.Graph()
    G.add_edges_from(get_neighbours(m))

    neighs_ns = {}
    for n in G.nodes:
        neighs_ns[n] = len(list(G.neighbors(n)))

    ccs = m.cell_centers().points
    areas = m.compute_cell_sizes().cell_data["Area"]

    rows = []
    for i in range(m.n_cells):
        print(i)
        _, ax1, ax2, _ = get_obb(m.extract_cells([i]))
        ar, max_length, min_length = aspect_ratio(ax1, ax2)
        ctr = ccs[i, :]
        if i not in neighs_ns:
            continue
        nns = neighs_ns.get(i, 0)
        area = areas[i]
        sides = len(m.cell_points(i))


        rows.append(",".join([str(elem) for elem in [ctr[0],
                                                     ctr[1],
                                                     ctr[2],
                                                     nns,
                                                     area,
                                                     sides,
                                                     ar,
                                                     max_length,
                                                     min_length]]))

    header = "x,y,z,neighbours,area,sides,aspect_ratio,max_length,min_length"
    rows = [header] + rows
    
    with open("neighbour_info.csv", "w") as fout:
        fout.write("\n".join(rows))

    return


def go(fn):
    m = get_mesh_vars(fn)
    write_neighbours_info(m)
        

if __name__ == "__main__":
    go(sys.argv[1])
