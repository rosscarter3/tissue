//
// 2D turgor pressure: force on each cell's vertices along the gradient of the
// cell area (shoelace), i.e. an internal pressure inflating the polygon.
// Ported from legacy pressure2D.cc (Pressure2D::AreaPotential) with one
// robustness fix: the force is oriented by the cell's signed area, so P > 0
// always inflates regardless of the tissue's sorting orientation (legacy
// depended on the arbitrary orientation picked during sorting).
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class Pressure2DAreaPotential : public Reaction {
public:
  Pressure2DAreaPotential(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 || (p[1] != 0.0 && p[1] != 1.0))
      throw std::runtime_error(
          "Pressure2D::AreaPotential: uses two parameters, P_force and "
          "normalizeVolumeFlag (0 or 1).");
    configure("Pressure2D::AreaPotential", p, i, 2, {}, {"P_force", "f_V_norm"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error("Pressure2D::AreaPotential requires 2D.");
    const double pForce = parameter(0);
    const bool normalize = parameter(1) == 1.0;
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numVertex();
            double signedArea = T.cellVolume(c, vertexData, true);
            double factor = 0.5 * pForce;
            if (normalize)
              factor /= std::fabs(signedArea);
            if (signedArea < 0.0)
              factor = -factor; // orientation fix: P>0 always inflates
            for (size_t k = 0; k < n; ++k) {
              size_t v = cell.vertices[k];
              size_t vPlus = cell.vertices[(k + 1) % n];
              size_t vMinus = cell.vertices[(k + n - 1) % n];
              out[v][0] += factor * (vertexData[vPlus][1] - vertexData[vMinus][1]);
              out[v][1] += factor * (vertexData[vMinus][0] - vertexData[vPlus][0]);
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure2DAreaPotential, "Pressure2D::AreaPotential")

} // namespace

// The remaining three pressure rules share a different discretisation: instead
// of the shoelace gradient, each wall contributes the gradient of the area of
// the triangle it forms with the cell centre (base b = wall length, height h =
// distance from centre to the wall line). Summed over walls this is the same
// area, but the force is distributed differently - each vertex is pushed away
// from the centre along the local outward normal rather than along the
// perpendicular to its two neighbours.
//
// Unlike Pressure2D::AreaPotential this form is sign-safe already: swapping a
// wall's two vertices maps the v1 expression onto the v2 one exactly, and the
// only use of the cell area is |A| as a normaliser, so no orientation fix is
// needed here.
namespace {

// Gradient of the centre-triangle area for one wall, accumulated onto both of
// its vertices. Transcribed from legacy pressure2D.cc; `fac` carries the
// pressure and any 1/|A| normalisation.
inline void triWallForce(const Matrix &vertexData, const Vec3 &xCenter,
                         size_t v1I, size_t v2I, double fac, Matrix &out) {
  double n[2], dx[2];
  double b = 0.0;
  for (size_t d = 0; d < 2; ++d) {
    n[d] = vertexData[v2I][d] - vertexData[v1I][d];
    b += n[d] * n[d];
    dx[d] = xCenter[d] - 0.5 * (vertexData[v1I][d] + vertexData[v2I][d]);
  }
  b = std::sqrt(b);
  for (size_t d = 0; d < 2; ++d)
    n[d] /= b;
  const double bInv = 1.0 / b;
  const double ndx = n[0] * dx[0] + n[1] * dx[1];
  const double h = std::sqrt(dx[0] * dx[0] + dx[1] * dx[1] - ndx * ndx);
  const double hInv = 1.0 / h;

  out[v1I][0] += fac * (0.5 * b * hInv *
                            (-dx[0] - 2.0 * ndx *
                                          (-n[1] * n[1] * bInv * dx[0] -
                                           0.5 * n[0] + n[0] * n[1] * bInv * dx[1])) -
                        h * n[0]);
  out[v2I][0] += fac * (0.5 * b * hInv *
                            (-dx[0] - 2.0 * ndx *
                                          (n[1] * n[1] * bInv * dx[0] -
                                           0.5 * n[0] - n[0] * n[1] * bInv * dx[1])) +
                        h * n[0]);
  out[v1I][1] += fac * (0.5 * b * hInv *
                            (-dx[1] - 2.0 * ndx *
                                          (-n[0] * n[0] * bInv * dx[1] -
                                           0.5 * n[1] + n[0] * n[1] * bInv * dx[0])) -
                        h * n[1]);
  out[v2I][1] += fac * (0.5 * b * hInv *
                            (-dx[1] - 2.0 * ndx *
                                          (n[0] * n[0] * bInv * dx[1] -
                                           0.5 * n[1] - n[0] * n[1] * bInv * dx[0])) +
                        h * n[1]);
}

bool touchesBackground(const Tissue &T, size_t c) {
  for (size_t w : T.cell(c).walls)
    if (Tissue::isBackground(T.wall(w).otherCell(c)))
      return true;
  return false;
}


class Pressure2DAreaPotentialTri : public Reaction {
public:
  Pressure2DAreaPotentialTri(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error(
          "Pressure2D::AreaPotentialTri: uses two or three parameters, "
          "P_force, flag_Vnorm (0 or 1), [flag_internalCellsOnly (0 or 1)].");
    if (p[1] != 0.0 && p[1] != 1.0)
      throw std::runtime_error(
          "Pressure2D::AreaPotentialTri: flag_Vnorm must be 0 or 1.");
    if (p.size() > 2 && p[2] != 0.0 && p[2] != 1.0)
      throw std::runtime_error("Pressure2D::AreaPotentialTri: "
                               "flag_internalCellsOnly must be 0 or 1.");
    std::vector<std::string> names{"P_force", "f_V_norm"};
    if (p.size() > 2)
      names.push_back("f_internal");
    configure("Pressure2D::AreaPotentialTri", p, i, p.size(), {}, names);
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error("Pressure2D::AreaPotentialTri requires 2D.");
    const bool internalOnly = numParameter() > 2 && parameter(2) == 1.0;
    const bool normalize = parameter(1) == 1.0;
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            if (internalOnly && touchesBackground(T, c))
              continue;
            const Vec3 xCenter = T.cellPosition(c, vertexData);
            double fac = 0.5 * parameter(0);
            if (normalize)
              fac /= T.cellVolume(c, vertexData);
            for (size_t w : T.cell(c).walls)
              triWallForce(vertexData, xCenter, T.wall(w).vertex1,
                           T.wall(w).vertex2, fac, out);
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure2DAreaPotentialTri,
                         "Pressure2D::AreaPotentialTri")

// As AreaPotentialTri, but pressure acts only within X_th of the tissue's
// leading edge along one axis - a growing apex with a quiescent base.
class Pressure2DAreaPotentialTriSpatialThreshold : public Reaction {
public:
  Pressure2DAreaPotentialTriSpatialThreshold(const ParameterList &p,
                                             const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "Pressure2D::AreaPotentialTriSpatialThreshold: one index "
          "(the axis the threshold is measured along).");
    configure("Pressure2D::AreaPotentialTriSpatialThreshold", p, i, 2, {1},
              {"K_force", "X_th"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t maxDim = variableIndex(0, 0);
    if (maxDim >= vertexData.cols())
      throw std::runtime_error("Pressure2D::AreaPotentialTriSpatialThreshold: "
                               "axis index outside the tissue dimension.");
    double max = vertexData[0][maxDim];
    for (size_t v = 1; v < T.numVertex(); ++v)
      max = std::max(max, vertexData[v][maxDim]);

    const double fac = 0.5 * parameter(0);
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const Vec3 xCenter = T.cellPosition(c, vertexData);
            if (max - xCenter[maxDim] >= parameter(1))
              continue;
            for (size_t w : T.cell(c).walls)
              triWallForce(vertexData, xCenter, T.wall(w).vertex1,
                           T.wall(w).vertex2, fac, out);
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure2DAreaPotentialTriSpatialThreshold,
                         "Pressure2D::AreaPotentialTriSpatialThreshold")

// Pressure proportional to how far the cell is below a target area held in a
// cell variable (which another reaction typically grows), rather than a
// constant turgor. With f_no_contraction set, a cell over its target is left
// alone instead of being pulled back in.
class Pressure2DAreaPotentialTargetArea : public Reaction {
public:
  Pressure2DAreaPotentialTargetArea(const ParameterList &p,
                                    const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Pressure2D::AreaPotentialTargetArea: one "
                               "index (the target-area cell variable).");
    configure("Pressure2D::AreaPotentialTargetArea", p, i, 2, {1},
              {"P", "f_no_contraction"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "Pressure2D::AreaPotentialTargetArea requires 2D.");
    const size_t target = variableIndex(0, 0);
    const double pressure = parameter(0);
    const bool noContraction = parameter(1) != 0.0;
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numVertex();
            const double signedArea = T.cellVolume(c, vertexData, true);
            const double sa = signedArea < 0.0 ? -1.0 : 1.0;
            const double area = std::fabs(signedArea);
            if (noContraction && area > cellData[c][target])
              continue;
            const double f = pressure * (1.0 - area / cellData[c][target]);
            for (size_t k = 0; k < n; ++k) {
              const size_t v = cell.vertices[k];
              auto prev = vertexData[cell.vertices[(k + n - 1) % n]];
              auto next = vertexData[cell.vertices[(k + 1) % n]];
              out[v][0] += f * sa * 0.5 * (next[1] - prev[1]);
              out[v][1] += f * sa * 0.5 * (prev[0] - next[0]);
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure2DAreaPotentialTargetArea,
                         "Pressure2D::AreaPotentialTargetArea")

} // namespace
} // namespace tissue
