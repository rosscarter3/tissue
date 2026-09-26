//
// Diagnostic::SurfaceVTK: write the tissue's face as a 3D surface, with the
// fields that make its failures visible.
//
// The print format carries the cell variables a tissue was initiated with and
// nothing a reaction adds afterwards, so a center-triangulated centre or a
// CellMesh triangle list never reaches the output at all. Every tool
// downstream therefore sees cell outlines and infers the rest, which is how a
// shell that was wrinkling at the mesh scale, and a fan that had inverted in
// every lobed cell, both went unnoticed while their summary statistics looked
// reasonable. This writes what is actually there.
//
// Legacy VTK PolyData, one file per interval, numbered so ParaView loads the
// set as a time series. Fields:
//
//   quality    radius ratio of the triangle: 1 equilateral, large is a
//              sliver, and the inverted ones are flagged separately because
//              they have no meaningful ratio
//   inverted   1 where the triangle has folded over
//   stretch    sqrt(area/restArea), so a face carrying load is visible
//   cellId     to colour cells apart
//   height     signed distance of a vertex from its own cell's plane, which
//              is corrugation before anything averages it away
//   interior   0 on the outline, 1 for a vertex the face mesh added
//
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class SurfaceVTK : public Reaction {
public:
  SurfaceVTK(const ParameterList &p, const IndexLevels &i) {
    configure("Diagnostic::SurfaceVTK", p, i, 2, {1},
              {"interval", "mode"}); // mode 0 = center-triangulation fan,
                                     //      1 = CellMesh triangle list
  }

  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              double h) override {
    time_ += h;
    if (wrote_ && time_ < next_)
      return;
    next_ = time_ + std::max(parameter(0), 1e-9);
    wrote_ = true;

    const std::size_t base = variableIndex(0, 0);
    const bool cellMesh = parameter(1) >= 0.5;
    const std::size_t dim = vertexData.cols();

    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<std::size_t, 3>> tris;
    std::vector<double> height, pcell, interior;
    std::vector<double> quality, inverted, stretch, ccell;

    for (std::size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const std::size_t n = cell.numVertex();
      if (n < 3 || cellData.rowSize(c) <= base + 2)
        continue;

      // The cell's own plane, so "height" means out of its surface rather
      // than out of whatever global axis happens to be up.
      double cen[3] = {0, 0, 0}, nrm[3] = {0, 0, 0};
      for (std::size_t k = 0; k < n; ++k)
        for (std::size_t d = 0; d < 3 && d < dim; ++d)
          cen[d] += vertexData[cell.vertices[k]][d] / double(n);
      for (std::size_t k = 0; k < n; ++k) {
        const auto a = vertexData[cell.vertices[k]];
        const auto b = vertexData[cell.vertices[(k + 1) % n]];
        const double az = dim > 2 ? a[2] : 0.0, bz = dim > 2 ? b[2] : 0.0;
        nrm[0] += (a[1] - b[1]) * (az + bz);
        nrm[1] += (az - bz) * (a[0] + b[0]);
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
      }
      double nl = std::sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2]);
      if (!(nl > 0.0)) { nrm[0] = nrm[1] = 0; nrm[2] = 1; nl = 1; }
      for (int d = 0; d < 3; ++d)
        nrm[d] /= nl;

      const std::size_t first = pts.size();
      std::vector<std::array<double, 3>> local;
      std::vector<double> isInner;
      for (std::size_t k = 0; k < n; ++k) {
        const auto v = vertexData[cell.vertices[k]];
        local.push_back({v[0], v[1], dim > 2 ? v[2] : 0.0});
        isInner.push_back(0.0);
      }

      std::vector<std::array<std::size_t, 3>> localTris;
      std::vector<std::array<double, 3>> localRest;
      if (cellMesh) {
        const std::size_t maxI = std::size_t(cellData[c][base]);
        const std::size_t maxT = std::size_t(cellData[c][base + 1]);
        const std::size_t m = std::size_t(cellData[c][base + 2]);
        const std::size_t ipos = base + 3;
        const std::size_t tcnt = ipos + 3 * maxI;
        if (cellData.rowSize(c) <= tcnt)
          continue;
        const std::size_t nt = std::size_t(cellData[c][tcnt]);
        const std::size_t codes = tcnt + 1, rests = codes + 3 * maxT;
        for (std::size_t k = 0; k < m; ++k) {
          local.push_back({cellData[c][ipos + 3*k], cellData[c][ipos + 3*k+1],
                           cellData[c][ipos + 3*k+2]});
          isInner.push_back(1.0);
        }
        for (std::size_t t = 0; t < nt; ++t) {
          localTris.push_back({std::size_t(cellData[c][codes + 3*t]),
                               std::size_t(cellData[c][codes + 3*t + 1]),
                               std::size_t(cellData[c][codes + 3*t + 2])});
          localRest.push_back({cellData[c][rests + 3*t],
                               cellData[c][rests + 3*t + 1],
                               cellData[c][rests + 3*t + 2]});
        }
      } else {
        local.push_back({cellData[c][base], cellData[c][base + 1],
                         cellData[c][base + 2]});
        isInner.push_back(1.0);
        for (std::size_t k = 0; k < n; ++k) {
          localTris.push_back({n, k, (k + 1) % n});
          localRest.push_back({0.0, 0.0, 0.0});
        }
      }

      for (std::size_t k = 0; k < local.size(); ++k) {
        pts.push_back(local[k]);
        double dz = 0.0;
        for (int d = 0; d < 3; ++d)
          dz += (local[k][d] - cen[d]) * nrm[d];
        height.push_back(dz);
        pcell.push_back(double(c));
        interior.push_back(isInner[k]);
      }
      for (std::size_t t = 0; t < localTris.size(); ++t) {
        const auto &tr = localTris[t];
        if (tr[0] >= local.size() || tr[1] >= local.size() ||
            tr[2] >= local.size())
          continue;
        tris.push_back({first + tr[0], first + tr[1], first + tr[2]});
        const auto &a = local[tr[0]], &b = local[tr[1]], &d2 = local[tr[2]];
        double u[3], v[3], cr[3];
        for (int d = 0; d < 3; ++d) { u[d] = b[d]-a[d]; v[d] = d2[d]-a[d]; }
        cr[0] = u[1]*v[2]-u[2]*v[1];
        cr[1] = u[2]*v[0]-u[0]*v[2];
        cr[2] = u[0]*v[1]-u[1]*v[0];
        const double a2 = std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        const double area = 0.5 * a2;
        auto len = [](const std::array<double,3> &p,
                      const std::array<double,3> &q) {
          return std::sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) +
                           (p[2]-q[2])*(p[2]-q[2]));
        };
        const double la = len(b, d2), lb = len(d2, a), lc = len(a, b);
        const double sp = 0.5 * (la + lb + lc);
        double qy = 1e6;
        bool inv = area <= 1e-14;
        if (!inv) {
          const double R = la * lb * lc / (4.0 * area);
          const double rin = area / sp;
          qy = rin > 0.0 ? std::min(R / (2.0 * rin), 1e6) : 1e6;
        }
        quality.push_back(qy);
        inverted.push_back(inv ? 1.0 : 0.0);
        const auto &r = localRest[t];
        double st = 0.0;
        if (r[0] > 0 && r[1] > 0 && r[2] > 0) {
          const double s = 0.5 * (r[0] + r[1] + r[2]);
          const double val = s * (s-r[0]) * (s-r[1]) * (s-r[2]);
          if (val > 0)
            st = std::sqrt(area / std::sqrt(val));
        }
        stretch.push_back(st);
        ccell.push_back(double(c));
      }
    }

    std::ostringstream name;
    name << "surface_" << std::setw(4) << std::setfill('0') << index_++
         << ".vtk";
    std::ofstream f(name.str());
    if (!f) {
      std::fprintf(stderr, "Diagnostic::SurfaceVTK: cannot write %s\n",
                   name.str().c_str());
      return;
    }
    f << "# vtk DataFile Version 3.0\npavement tissue t=" << time_
      << "\nASCII\nDATASET POLYDATA\nPOINTS " << pts.size() << " float\n";
    for (const auto &p : pts)
      f << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    f << "POLYGONS " << tris.size() << ' ' << 4 * tris.size() << '\n';
    for (const auto &t : tris)
      f << "3 " << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
    auto scalars = [&](const char *lab, const std::vector<double> &v) {
      f << "SCALARS " << lab << " float 1\nLOOKUP_TABLE default\n";
      for (double x : v)
        f << (std::isfinite(x) ? x : 0.0) << '\n';
    };
    f << "POINT_DATA " << pts.size() << '\n';
    scalars("height", height);
    scalars("cellId", pcell);
    scalars("interior", interior);
    f << "CELL_DATA " << tris.size() << '\n';
    scalars("quality", quality);
    scalars("inverted", inverted);
    scalars("stretch", stretch);
    scalars("cellId", ccell);
    std::fprintf(stderr, "Diagnostic::SurfaceVTK: %s (%zu points, %zu "
                 "triangles)\n", name.str().c_str(), pts.size(), tris.size());
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}

private:
  double time_ = 0.0, next_ = 0.0;
  bool wrote_ = false;
  std::size_t index_ = 0;
};

} // namespace

TISSUE_REGISTER_REACTION(SurfaceVTK, "Diagnostic::SurfaceVTK")

} // namespace tissue
