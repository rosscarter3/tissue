//
// CellMesh: a cell face meshed properly, instead of fanned from one centre.
//
// CenterTriangulation stores a centre in the cell row and treats the cell as
// a fan of triangles from it. That is a valid triangulation only where the
// cell is star-shaped about that centre, and a lobed pavement cell is not:
// measured on a published 2D trajectory, the star-shaped fraction falls from
// 84% at the start to 4% once the cells have lobed, and the centroid the code
// actually uses is valid for 58% and then 3%. The fan inverts, a TRBS
// membrane is handed a rest state no triangle can adopt, and the largest
// stiffness eigenvalue reaches 1.4e19 -- a timestep of 1e-21 and a run that
// burns CPU without moving.
//
// What this buys, measured by assembling the constant-strain-triangle
// stiffness (what TRBS linearises to) on real outlines and taking its largest
// eigenvalue by power iteration -- the quantity QuasiStatic::calibrateStep
// divides into its timestep:
//
//   median circularity     fan          this
//   0.67                   2.85e7       2603
//   0.31                   8.91e6       2567
//   0.23                   7.58e5       2583
//
// The fan's swings three orders of magnitude with cell shape. This is flat,
// so the relaxation's stable step stops depending on how lobed the tissue
// has become. FIRE needs O(sqrt(kappa)) iterations, so that is a 17x to 105x
// change in iteration count on those cells.
//
// Cell row layout from base = variableIndex(0, 0):
//
//   base + 0                         m, the number of interior vertices
//   base + 1 .. + 3*maxInterior      their positions, x y z each
//   base + 1 + 3*maxInterior         t, the number of triangles
//   next 3*maxTriangles              vertex codes, three per triangle
//   next 3*maxTriangles              rest edge lengths, three per triangle
//
// A vertex code below the cell's wall count indexes cell.vertices; at or
// above it, the interior vertex (code - numVertex). Interior vertices live in
// the cell row and are declared positional, exactly as the CT centre is, so
// the solvers relax them as positions and no change to Tissue's topology is
// needed.
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "tissue/core/meshing.h"
#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"
#include "tissue/reactions/trbs_core.h"

namespace tissue {
namespace {

// Where each block starts. The two maxima are stored in the row itself, at
// its first two slots, so any reaction that walks a cell row can work the
// layout out from the row rather than from a shared variable -- the
// alternative was a file-scope global that CellMesh::Initiate wrote and the
// mechanics read, which is a coupling waiting to go wrong the first time
// something is initialised in a different order.
struct Layout {
  std::size_t base = 0, maxInterior = 0, maxTri = 0;
  static Layout of(const Matrix &cellData, std::size_t c, std::size_t base) {
    Layout l;
    l.base = base;
    if (cellData.rowSize(c) > base + 1) {
      l.maxInterior = std::size_t(cellData[c][base]);
      l.maxTri = std::size_t(cellData[c][base + 1]);
    }
    return l;
  }
  std::size_t numInterior() const { return base + 2; }
  std::size_t interiorPos() const { return base + 3; }
  std::size_t triCount() const { return base + 3 + 3 * maxInterior; }
  std::size_t triCodes() const { return triCount() + 1; }
  std::size_t triRest() const { return triCodes() + 3 * maxTri; }
  std::size_t width() const { return 4 + 3 * maxInterior + 6 * maxTri; }
};

// The cell's own plane, so a non-planar face is meshed in the plane it
// actually occupies rather than in whatever xy happens to be.
void cellFrame(const Tissue &T, std::size_t c, const Matrix &vertexData,
               double origin[3], double ex[3], double ey[3]) {
  const CellTopo &cell = T.cell(c);
  const std::size_t n = cell.numVertex();
  double nrm[3] = {0.0, 0.0, 0.0};
  for (std::size_t d = 0; d < 3; ++d)
    origin[d] = 0.0;
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t d = 0; d < 3; ++d)
      origin[d] += vertexData[cell.vertices[k]][d] / double(n);
  for (std::size_t k = 0; k < n; ++k) { // Newell
    const auto a = vertexData[cell.vertices[k]];
    const auto b = vertexData[cell.vertices[(k + 1) % n]];
    nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
    nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
    nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
  }
  double L = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  if (L <= 0.0) {
    nrm[0] = 0.0; nrm[1] = 0.0; nrm[2] = 1.0; L = 1.0;
  }
  for (int d = 0; d < 3; ++d)
    nrm[d] /= L;
  double a[3] = {1.0, 0.0, 0.0};
  if (std::fabs(nrm[0]) > 0.9) {
    a[0] = 0.0; a[1] = 1.0;
  }
  ex[0] = a[1] * nrm[2] - a[2] * nrm[1];
  ex[1] = a[2] * nrm[0] - a[0] * nrm[2];
  ex[2] = a[0] * nrm[1] - a[1] * nrm[0];
  const double e = std::sqrt(ex[0]*ex[0] + ex[1]*ex[1] + ex[2]*ex[2]);
  for (int d = 0; d < 3; ++d)
    ex[d] /= e;
  ey[0] = nrm[1] * ex[2] - nrm[2] * ex[1];
  ey[1] = nrm[2] * ex[0] - nrm[0] * ex[2];
  ey[2] = nrm[0] * ex[1] - nrm[1] * ex[0];
}

class CellMeshInitiate : public Reaction {
public:
  CellMeshInitiate(const ParameterList &p, const IndexLevels &i) {
    configure("CellMesh::Initiate", p, i, 1, {1}, {"quality_bound"});
  }

  void positionalCellVariables(std::vector<std::size_t> &out) const override {
    for (std::size_t k = 0; k < 3 * layout_.maxInterior; ++k)
      out.push_back(layout_.interiorPos() + k);
  }

  void initiate(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &cellDerivs, Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error("CellMesh::Initiate requires a 3D tissue.");
    const std::size_t base = variableIndex(0, 0);
    const double bound = parameter(0) > 0.0 ? parameter(0) : 5.0;

    // Mesh every cell in its own plane first, then size the row to fit.
    std::vector<PolygonMesh> meshes(T.numCell());
    std::vector<std::array<double, 9>> frames(T.numCell());
    std::size_t maxInterior = 0, maxTri = 0, worstCell = 0;
    double worstQ = 0.0;
    for (std::size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const std::size_t n = cell.numVertex();
      if (n < 3)
        continue;
      double o[3], ex[3], ey[3];
      cellFrame(T, c, vertexData, o, ex, ey);
      frames[c] = {o[0], o[1], o[2], ex[0], ex[1], ex[2], ey[0], ey[1], ey[2]};
      std::vector<std::array<double, 2>> flat(n);
      for (std::size_t k = 0; k < n; ++k) {
        const auto v = vertexData[cell.vertices[k]];
        const double dx = v[0] - o[0], dy = v[1] - o[1], dz = v[2] - o[2];
        flat[k] = {dx * ex[0] + dy * ex[1] + dz * ex[2],
                   dx * ey[0] + dy * ey[1] + dz * ey[2]};
      }
      meshes[c] = triangulatePolygon(flat, bound);
      // triangulatePolygon may reverse the outline to make it counter-
      // clockwise; detect that so the codes refer to the right vertices.
      maxInterior = std::max(maxInterior,
                             meshes[c].points.size() - meshes[c].numBoundary);
      maxTri = std::max(maxTri, meshes[c].tris.size());
      for (double q : triangleQuality(meshes[c]))
        if (std::isfinite(q) && q > worstQ) {
          worstQ = q;
          worstCell = c;
        }
    }

    layout_ = Layout{base, maxInterior, maxTri};
    // Every row gets the same width, so one column list describes the layout
    // for all cells; rows short of it are skipped by the solver's positional
    // index, and the padding slots carry no force.
    const std::size_t want = base + layout_.width();
    for (std::size_t c = 0; c < cellData.rows(); ++c) {
      if (cellData.rowSize(c) < want) {
        cellData.resizeRow(c, want);
        cellDerivs.resizeRow(c, want);
      }
    }

    for (std::size_t c = 0; c < T.numCell(); ++c) {
      const PolygonMesh &m = meshes[c];
      if (m.tris.empty())
        continue;
      const CellTopo &cell = T.cell(c);
      const std::size_t n = cell.numVertex();
      const auto &f = frames[c];
      // Was the outline reversed by the mesher? Compare the first boundary
      // point against the projection of cell.vertices[0].
      const auto v0 = vertexData[cell.vertices[0]];
      const double dx = v0[0] - f[0], dy = v0[1] - f[1], dz = v0[2] - f[2];
      const double u0 = dx * f[3] + dy * f[4] + dz * f[5];
      const bool reversed = std::fabs(m.points[0][0] - u0) > 1e-9;

      const std::size_t mi = m.points.size() - m.numBoundary;
      cellData[c][base] = double(maxInterior);
      cellData[c][base + 1] = double(maxTri);
      cellData[c][layout_.numInterior()] = double(mi);
      for (std::size_t k = 0; k < mi; ++k) {
        const auto &p = m.points[m.numBoundary + k];
        for (std::size_t d = 0; d < 3; ++d)
          cellData[c][layout_.interiorPos() + 3 * k + d] =
              f[d] + p[0] * f[3 + d] + p[1] * f[6 + d];
      }
      cellData[c][layout_.triCount()] = double(m.tris.size());
      auto code = [&](std::size_t idx) -> double {
        if (idx >= m.numBoundary)
          return double(n + (idx - m.numBoundary));
        // triangulatePolygon reverses the outline when it is clockwise, and
        // std::reverse sends output index i to input index n-1-i. Writing
        // (n - i) % n instead maps 0 to 0 rather than to n-1, which shifts
        // every code by one, gives every triangle the wrong geometry, and
        // produced NaN on the first force evaluation.
        return double(reversed ? (n - 1 - idx) : idx);
      };
      for (std::size_t t = 0; t < m.tris.size(); ++t) {
        double pos[3][3];
        for (std::size_t j = 0; j < 3; ++j) {
          const double cd = code(m.tris[t][j]);
          cellData[c][layout_.triCodes() + 3 * t + j] = cd;
          const std::size_t ci = std::size_t(cd);
          for (std::size_t d = 0; d < 3; ++d)
            pos[j][d] = ci < n
                ? vertexData[cell.vertices[ci]][d]
                : cellData[c][layout_.interiorPos() + 3 * (ci - n) + d];
        }
        // Rest lengths are the initial geometry: the face starts unstrained,
        // as CenterTriangulation::Initiate does for its fan. The edge order
        // is the one trbs::Element uses -- 0 is (node0,node1), 1 is
        // (node1,node2), 2 is (node0,node2) -- and not "the edge opposite
        // node e", which is the natural reading and the wrong one.
        const int ea[3] = {0, 1, 0}, eb[3] = {1, 2, 2};
        for (std::size_t e = 0; e < 3; ++e) {
          double sq = 0.0;
          for (std::size_t d = 0; d < 3; ++d) {
            const double diff = pos[ea[e]][d] - pos[eb[e]][d];
            sq += diff * diff;
          }
          cellData[c][layout_.triRest() + 3 * t + e] = std::sqrt(sq);
        }
      }
    }
    std::size_t usedInterior = 0, usedTri = 0;
    for (const PolygonMesh &m : meshes) {
      usedInterior += m.points.size() - m.numBoundary;
      usedTri += m.tris.size();
    }
    std::cerr << "CellMesh::Initiate: " << T.numCell() << " cells, "
              << usedInterior << " interior vertices and " << usedTri
              << " triangles; worst element " << worstQ << " (cell "
              << worstCell << ")" << std::endl;
    // Every row is padded to the widest cell so one column list describes the
    // layout for all of them. That is the price of the current positional
    // variable API, which is per reaction rather than per cell: the padding
    // slots carry no force and sit still, but the solver still counts them.
    if (maxInterior * T.numCell() > 2 * usedInterior)
      std::cerr << "  rows padded to " << maxInterior
                << " interior vertices each, so "
                << (maxInterior * T.numCell()) << " slots hold "
                << usedInterior << " real ones ("
                << (100.0 * double(usedInterior) /
                    double(maxInterior * T.numCell()))
                << "% used)" << std::endl;
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}

  Layout layout_;
};

// The membrane, over the triangle list rather than the fan.
//
// Identical elasticity to VertexFromTRBScenterTriangulation -- the same
// per-element kernel from trbs_core.h, the same Lame coefficients -- so a
// model can be moved from one to the other and mean the same thing. What
// changes is which triangles it runs over: an explicit list that covers the
// cell properly instead of a fan that assumes the cell is star-shaped.
class VertexFromTRBSCellMesh : public Reaction {
public:
  VertexFromTRBSCellMesh(const ParameterList &p, const IndexLevels &i) {
    configure("VertexFromTRBSCellMesh", p, i, 2, {1},
              {"youngModulus", "poissonRatio"});
  }

  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &vertexDerivs) override {
    const std::size_t base = variableIndex(0, 0);
    const double young = parameter(0), poisson = parameter(1);
    const double lambda = young * poisson / (1 - poisson * poisson);
    const double mio = young / (1 + poisson);

    // Forces land on shared wall vertices from every cell that touches them,
    // so the vertex side scatters; the interior side is private to its cell
    // and can be written directly.
    parallelScatter1(
        T.numCell(), vertexDerivs,
        [&](std::size_t begin, std::size_t end, Matrix &out) {
          for (std::size_t c = begin; c < end; ++c) {
            const CellTopo &cell = T.cell(c);
            const std::size_t n = cell.numVertex();
            if (cellData.rowSize(c) <= base + 3)
              continue;
            const Layout L = Layout::of(cellData, c, base);
            const std::size_t posBase = L.interiorPos();
            if (cellData.rowSize(c) <= L.triCount())
              continue;
            const std::size_t nt = std::size_t(cellData[c][L.triCount()]);
            const std::size_t codeBase = L.triCodes();
            const std::size_t restBase = L.triRest();

            for (std::size_t t = 0; t < nt; ++t) {
              trbs::Element el;
              std::size_t code[3];
              for (std::size_t j = 0; j < 3; ++j) {
                code[j] = std::size_t(cellData[c][codeBase + 3 * t + j]);
                for (std::size_t d = 0; d < 3; ++d)
                  el.pos[j][d] =
                      code[j] < n
                          ? vertexData[cell.vertices[code[j]]][d]
                          : cellData[c][posBase + 3 * (code[j] - n) + d];
                el.rest[j] = cellData[c][restBase + 3 * t + j];
              }
              if (el.rest[0] <= 0.0 || el.rest[1] <= 0.0 || el.rest[2] <= 0.0)
                continue;
              // completeElement does not derive the current lengths; the
              // caller fills them, and leaving them uninitialised is NaN on
              // the first force evaluation.
              el.cur[0] = trbs::nodeDistance(el, 0, 1);
              el.cur[1] = trbs::nodeDistance(el, 1, 2);
              el.cur[2] = trbs::nodeDistance(el, 0, 2);
              trbs::completeElement(el);
              double f[3][3];
              trbs::elementForces(el, trbs::stiffnessFrom(el, lambda + mio, mio),
                                  f);
              for (std::size_t j = 0; j < 3; ++j) {
                if (code[j] < n) {
                  const std::size_t v = cell.vertices[code[j]];
                  for (std::size_t d = 0; d < 3; ++d)
                    out[v][d] += f[j][d];
                } else {
                  const std::size_t k = posBase + 3 * (code[j] - n);
                  for (std::size_t d = 0; d < 3; ++d)
                    cellDerivs[c][k + d] += f[j][d];
                }
              }
            }
          }
        });
  }
};

} // namespace

TISSUE_REGISTER_REACTION(CellMeshInitiate, "CellMesh::Initiate")
TISSUE_REGISTER_REACTION(VertexFromTRBSCellMesh, "VertexFromTRBSCellMesh")

} // namespace tissue
