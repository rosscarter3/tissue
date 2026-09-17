#include "tissue/io/vtk.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace tissue::io {

namespace {

constexpr double kWallRelativeThickness = 0.1;

std::string numberedName(const std::string &base, size_t count) {
  auto dot = base.find_last_of('.');
  std::ostringstream fileid;
  fileid.fill('0');
  fileid.width(6);
  fileid << count;
  return base.substr(0, dot) + fileid.str() + base.substr(dot);
}

double safeCellVar(const Matrix &cellState, size_t c, size_t j) {
  // Legacy indexes variables 0..3 unconditionally (UB when fewer exist);
  // v2 pads with zeros instead.
  return j < cellState.rowSize(c) ? cellState[c][j] : 0.0;
}

struct Point3 {
  double x = 0, y = 0, z = 0;
  void displaceTowards(const Point3 &p, double d) {
    x += d * (p.x - x);
    y += d * (p.y - y);
    z += d * (p.z - z);
  }
};

Point3 vertexPoint(const Tissue &T, size_t v) {
  const Matrix &pos = T.vertexState();
  Point3 p;
  p.x = pos[v][0];
  p.y = pos[v][1];
  if (pos.cols() > 2)
    p.z = pos[v][2];
  return p;
}

Point3 cellCenter(const Tissue &T, size_t c) {
  Vec3 center = T.cellPosition(c, T.vertexState());
  return Point3{center[0], center[1], center[2]};
}

void writeHeader(std::ostream &os) {
  os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n"
     << "<UnstructuredGrid>\n";
}
void writeFooter(std::ostream &os) {
  os << "</UnstructuredGrid>\n"
     << "</VTKFile>\n";
}
void writePieceHeader(std::ostream &os, size_t nPts, size_t nCell) {
  os << "<Piece  NumberOfPoints=\"";
  os.width(8);
  os << nPts;
  os << "\" NumberOfCells=\"";
  os.width(8);
  os << nCell << "\">";
  os << "\n";
}

void writeCellData(std::ostream &os, const Tissue &T) {
  const Matrix &cs = T.cellState();
  const size_t nvars = T.numCellVariable();
  os << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"cell "
        "vector\" format=\"ascii\">\n";
  for (size_t c = 0; c < T.numCell(); ++c)
    os << safeCellVar(cs, c, 0) << " " << safeCellVar(cs, c, 1) << " "
       << safeCellVar(cs, c, 2) << "\n";
  os << "</DataArray>\n";
  os << "<DataArray type=\"Float64\" Name=\"cell vector length\" "
        "format=\"ascii\">\n";
  for (size_t c = 0; c < T.numCell(); ++c)
    os << safeCellVar(cs, c, 3) << " ";
  os << "\n"
     << "</DataArray>\n";
  for (size_t j = 4; j < nvars; ++j) {
    os << "<DataArray type=\"Float64\" Name=\"cell variable " << j
       << "\" format=\"ascii\">\n";
    for (size_t c = 0; c < T.numCell(); ++c)
      os << cs[c][j] << " ";
    os << "\n"
       << "</DataArray>\n";
  }
}

void writeCells2(std::ostream &os, const Tissue &T) {
  size_t npts = 0;
  bool triangles = true;
  for (size_t c = 0; c < T.numCell(); ++c) {
    size_t nvrt = T.cell(c).numVertex();
    if (nvrt != 3)
      triangles = false;
    npts += nvrt;
  }
  writePieceHeader(os, npts, T.numCell());
  const size_t dim = T.dimension();
  // Points: each cell's vertices displaced towards the cell center.
  os << "<Points>\n"
     << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";
  for (size_t c = 0; c < T.numCell(); ++c) {
    Point3 center = cellCenter(T, c);
    for (size_t v : T.cell(c).vertices) {
      Point3 p = vertexPoint(T, v);
      p.displaceTowards(center, kWallRelativeThickness);
      if (dim == 2)
        os << p.x << " " << p.y << " 0.0\n";
      else
        os << p.x << " " << p.y << " " << p.z << "\n";
    }
  }
  os << "</DataArray>\n"
     << "</Points>\n";
  // Connectivity/offsets/types.
  os << "<Cells>\n"
     << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  size_t count = 0;
  for (size_t c = 0; c < T.numCell(); ++c) {
    for (size_t k = 0; k < T.cell(c).numVertex(); ++k, ++count)
      os << count << " ";
    os << "\n";
  }
  os << "</DataArray>\n"
     << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  size_t totalOffset = 0;
  for (size_t c = 0; c < T.numCell(); ++c) {
    totalOffset += T.cell(c).numVertex();
    os << totalOffset << " ";
  }
  os << "\n"
     << "</DataArray>\n"
     << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  const int cellType = triangles ? 5 : 7;
  for (size_t c = 0; c < T.numCell(); ++c)
    os << cellType << " ";
  os << "\n";
  os << "</DataArray>\n"
     << "</Cells>\n";
  os << "<CellData Scalars=\"cell variable 4\" Vectors=\"cell vector\">\n";
  writeCellData(os, T);
  os << "</CellData>\n";
  os << "</Piece>\n";
}

// Builds, per cell, its vertex cycle in the order induced by consecutive
// walls (legacy write_wall_point_geometry2 traversal).
std::vector<size_t> cellWallVertexOrder(const Tissue &T, size_t c) {
  const CellTopo &cell = T.cell(c);
  std::vector<size_t> verts;
  verts.reserve(cell.numWall());
  const Wall &w1 = T.wall(cell.walls[0]);
  const Wall &w2 = T.wall(cell.walls[1]);
  size_t w1v1 = w1.vertex1, w1v2 = w1.vertex2;
  size_t chain;
  if (w1v1 == w2.vertex1) {
    verts.push_back(w1v2);
    verts.push_back(w1v1);
    chain = w2.vertex2;
  } else if (w1v1 == w2.vertex2) {
    verts.push_back(w1v2);
    verts.push_back(w1v1);
    chain = w2.vertex1;
  } else if (w1v2 == w2.vertex1) {
    verts.push_back(w1v1);
    verts.push_back(w1v2);
    chain = w2.vertex2;
  } else if (w1v2 == w2.vertex2) {
    verts.push_back(w1v1);
    verts.push_back(w1v2);
    chain = w2.vertex1;
  } else {
    std::cerr << "vtk: consecutive cell walls do not share a vertex.\n";
    std::exit(EXIT_FAILURE);
  }
  for (size_t k = 2; k < cell.numWall(); ++k) {
    const Wall &w = T.wall(cell.walls[k]);
    if (chain == w.vertex1) {
      verts.push_back(w.vertex1);
      chain = w.vertex2;
    } else if (chain == w.vertex2) {
      verts.push_back(w.vertex2);
      chain = w.vertex1;
    } else {
      std::cerr << "vtk: consecutive cell walls do not share a vertex.\n";
      std::exit(EXIT_FAILURE);
    }
  }
  return verts;
}

void writeWalls(std::ostream &os, const Tissue &T, bool pairedVariables) {
  size_t ncell = 0;
  for (size_t c = 0; c < T.numCell(); ++c)
    ncell += T.cell(c).numWall();
  writePieceHeader(os, T.numVertex() + ncell, ncell);
  const size_t dim = T.dimension();
  const Matrix &pos = T.vertexState();
  const Matrix &ws = T.wallState();

  // Per-cell wall-chain vertex orders.
  std::vector<std::vector<size_t>> order(T.numCell());
  for (size_t c = 0; c < T.numCell(); ++c)
    order[c] = cellWallVertexOrder(T, c);

  // Points: all vertices, then per cell the displaced chain vertices.
  os << "<Points>\n"
     << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
        "format=\"ascii\">\n";
  for (size_t v = 0; v < T.numVertex(); ++v) {
    for (size_t d = 0; d < dim; ++d)
      os << pos[v][d] << " ";
    if (dim < 3)
      os << "0 ";
    os << "\n";
  }
  for (size_t c = 0; c < T.numCell(); ++c) {
    Point3 center = cellCenter(T, c);
    for (size_t v : order[c]) {
      Point3 p = vertexPoint(T, v);
      p.displaceTowards(center, kWallRelativeThickness);
      if (dim == 3)
        os << p.x << " " << p.y << " " << p.z << "\n";
      else
        os << p.x << " " << p.y << " 0.0\n";
    }
  }
  os << "</DataArray>\n"
     << "</Points>\n";

  // Quads: original vertex pair + the two displaced counterparts.
  os << "<Cells>\n"
     << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  size_t count = 0;
  const size_t offset = T.numVertex();
  for (size_t c = 0; c < T.numCell(); ++c) {
    const auto &verts = order[c];
    size_t nwall = verts.size() - 1;
    size_t temp = offset + count;
    for (size_t k = 0; k < nwall; ++k) {
      os << verts[k] << " ";
      os << verts[k + 1] << " " << temp + k + 1 << " " << temp + k << " ";
    }
    os << verts[nwall] << " " << verts[0] << " " << temp << " "
       << temp + nwall << " ";
    count += verts.size();
  }
  os << "\n";
  os << "</DataArray>\n"
     << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  size_t totalOffset = 0;
  for (size_t k = 0; k < count; ++k) {
    totalOffset += 4;
    os << totalOffset << " ";
  }
  os << "\n"
     << "</DataArray>\n"
     << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (size_t k = 0; k < count; ++k)
    os << 7 << " ";
  os << "\n";
  os << "</DataArray>\n"
     << "</Cells>\n";

  os << "<CellData Scalars=\"wall variable 0\">\n";
  // Wall lengths per cell wall.
  os << "<DataArray type=\"Float64\" Name=\"wall length\" format=\"ascii\">\n";
  for (size_t c = 0; c < T.numCell(); ++c)
    for (size_t w : T.cell(c).walls)
      os << ws[w][0] << " ";
  os << "\n"
     << "</DataArray>\n";
  const size_t nvars = ws.cols() - 1;
  if (!pairedVariables) {
    for (size_t j = 0; j < nvars; ++j) {
      os << "<DataArray type=\"Float64\" Name=\"wall variable " << j
         << "\" format=\"ascii\">\n";
      for (size_t c = 0; c < T.numCell(); ++c)
        for (size_t w : T.cell(c).walls)
          os << ws[w][1 + j] << " ";
      os << "\n"
         << "</DataArray>\n";
    }
  } else {
    // Paired structure: variable j/2 read from index j (cell1 side) or j+1
    // (cell2 side).
    for (size_t j = 0; j < nvars; j += 2) {
      os << "<DataArray type=\"Float64\" Name=\"wall variable " << j / 2
         << "\" format=\"ascii\">\n";
      for (size_t c = 0; c < T.numCell(); ++c) {
        for (size_t w : T.cell(c).walls) {
          const Wall &wall = T.wall(w);
          size_t jj;
          if (wall.cell1 == c)
            jj = j;
          else if (wall.cell2 == c)
            jj = j + 1;
          else {
            std::cerr << "vtk: wall does not report connection to the cell "
                         "which was accessed through\n";
            std::exit(EXIT_FAILURE);
          }
          os << ws[w][1 + jj] << " ";
        }
      }
      os << "\n"
         << "</DataArray>\n";
    }
  }
  os << "</CellData>\n";
  os << "</Piece>\n";
}

} // namespace

void writeFullPvd(const std::string &pvdFile, const std::string &cellFile,
                  const std::string &wallFile, size_t n) {
  const std::string bases[2] = {cellFile, wallFile};
  std::string filestart[2], extension[2];
  for (size_t i = 0; i < 2; ++i) {
    std::string fname = bases[i].substr(bases[i].find_last_of('/') + 1);
    auto dot = fname.find_last_of('.');
    filestart[i] = fname.substr(0, dot);
    extension[i] = fname.substr(dot);
  }
  std::ofstream pvd(pvdFile);
  pvd << "<?xml version=\"1.0\"?> " << std::endl;
  pvd << "<VTKFile type=\"Collection\" version=\"0.1\" > " << std::endl;
  pvd << "<Collection> " << std::endl;
  for (size_t i = 0; i < n; ++i) {
    std::ostringstream fileid;
    fileid.fill('0');
    fileid.width(6);
    fileid << i;
    for (size_t j = 0; j < 2; ++j) {
      pvd << "<DataSet timestep=\"" << i << "\" part=\"" << j << "\" file=\""
          << filestart[j] + fileid.str() + extension[j] << "\"/>" << std::endl;
    }
  }
  pvd << "</Collection> " << std::endl;
  pvd << "</VTKFile> " << std::endl;
}

void writeVtu(const Tissue &T, const std::string &cellFile,
              const std::string &wallFile, size_t count, bool singleWall) {
  {
    std::ofstream os(numberedName(cellFile, count));
    writeHeader(os);
    writeCells2(os, T);
    writeFooter(os);
  }
  {
    std::ofstream os(numberedName(wallFile, count));
    writeHeader(os);
    writeWalls(os, T, !singleWall);
    writeFooter(os);
  }
}

} // namespace tissue::io
