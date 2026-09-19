//
// Turgor pressure on a centre-triangulated cell, ported from legacy
// mechanical.cc (the CenterTriangulation::VertexFromCellPressure pair).
//
// Unlike the 2D forms in pressure.cpp, which push along the gradient of the
// polygon area, these treat each cell as a fan of triangles around a stored
// centre and push each edge outwards, away from that centre and within the
// triangle's plane. 3D only.
//
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// The outward push on one edge of a centre-triangulated cell.
//
// `dx` comes out as the unit vector from the cell centre to the edge midpoint
// with the component along the edge removed - so it is perpendicular to the
// edge and lies in the triangle's plane, which is the direction an internal
// pressure acts in. `wallLength` is the edge length that scales the force.
//
// Returns false when the geometry is degenerate (a zero-length edge, or a
// midpoint sitting on the centre so the direction is undefined). Legacy exits
// the process in both cases; the callers here throw instead, naming the cell.
bool edgePushDirection(const Matrix &vertexData, size_t v1, size_t v2,
                       const std::vector<double> &center,
                       std::vector<double> &dx, double &wallLength) {
  std::vector<double> wallVector(3);
  wallLength = 0.0;
  for (size_t d = 0; d < 3; ++d) {
    dx[d] = 0.5 * (vertexData[v1][d] + vertexData[v2][d]) - center[d];
    wallVector[d] = vertexData[v1][d] - vertexData[v2][d];
    wallLength += wallVector[d] * wallVector[d];
  }
  if (!(wallLength > 0.0))
    return false;
  wallLength = std::sqrt(wallLength);

  // Project out the along-edge component, leaving the in-plane normal.
  double along = 0.0;
  for (size_t d = 0; d < 3; ++d)
    along += dx[d] * wallVector[d];
  along /= wallLength * wallLength;

  double norm = 0.0;
  for (size_t d = 0; d < 3; ++d) {
    dx[d] -= along * wallVector[d];
    norm += dx[d] * dx[d];
  }
  if (!(norm > 0.0))
    return false;
  norm = 1.0 / std::sqrt(norm);
  for (size_t d = 0; d < 3; ++d)
    dx[d] *= norm;
  return true;
}

// Shared body of both reactions. `strength` gives the pressure for a cell,
// letting the two forms differ only in how that is computed.
template <class Strength>
void centerTriangulationPressure(Tissue &T, const Matrix &cellData,
                                 const Matrix &vertexData,
                                 Matrix &vertexDerivs, size_t centerIndex,
                                 bool normalizeVolume, const char *name,
                                 Strength strength) {
  if (vertexData.cols() != 3)
    throw std::runtime_error(std::string(name) +
                             ": assumes vertices in three dimensions.");
  std::vector<double> center(3), dx(3);
  for (size_t n = 0; n < T.numCell(); ++n) {
    const CellTopo &cell = T.cell(n);
    for (size_t d = 0; d < 3; ++d)
      center[d] = cellData[n][centerIndex + d];

    double cellFactor = 0.5 * strength(n);
    if (normalizeVolume)
      cellFactor /= std::fabs(T.cellVolume(n, vertexData));

    const size_t nv = cell.numVertex();
    for (size_t k = 0; k < nv; ++k) {
      const size_t v1 = cell.vertices[k];
      const size_t v2 = cell.vertices[(k + 1) % nv];
      double wallLength = 0.0;
      if (!edgePushDirection(vertexData, v1, v2, center, dx, wallLength))
        throw std::runtime_error(
            std::string(name) + ": degenerate edge in cell " +
            std::to_string(n) + " (zero length, or its midpoint coincides "
                                "with the stored cell centre).");
      // Each of the edge's two vertices gets the whole force, not half of
      // it - legacy's convention, kept.
      const double factor = cellFactor * wallLength;
      for (size_t d = 0; d < 3; ++d) {
        vertexDerivs[v1][d] += factor * dx[d];
        vertexDerivs[v2][d] += factor * dx[d];
      }
    }
  }
}

///
/// Constant turgor on a centre-triangulated cell, optionally scaled by a cell
/// concentration and/or divided by the cell volume.
///
///   CenterTriangulation::VertexFromCellPressure 2 1 2
///     K_force, normalizeVolumeFlag
///     center_index concentration_index
///
class CenterTriangulationVertexFromCellPressure : public Reaction {
public:
  CenterTriangulationVertexFromCellPressure(const ParameterList &p,
                                            const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 2)
      throw std::runtime_error(
          "CenterTriangulation::VertexFromCellPressure: level 0 = the start "
          "of the centre-triangulation cell variables, then the "
          "concentration index.");
    configure("CenterTriangulation::VertexFromCellPressure", p, i, 2, {2},
              {"K_force", "f_V_norm"});
    if (parameter(1) != 0.0 && parameter(1) != 1.0)
      throw std::runtime_error(
          "CenterTriangulation::VertexFromCellPressure: normalizeVolumeFlag "
          "must be 0 or 1.");
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    // Legacy reads a concentration index of 0 as "no concentration", so cell
    // variable 0 cannot be used as the concentration. Kept: a model written
    // against legacy that names index 0 means "constant".
    const size_t concIndex = variableIndex(0, 1);
    centerTriangulationPressure(
        T, cellData, vertexData, vertexDerivs, variableIndex(0, 0),
        parameter(1) == 1.0, id().c_str(), [&](size_t n) {
          return concIndex != 0 ? parameter(0) * cellData[n][concIndex]
                                : parameter(0);
        });
  }
};
TISSUE_REGISTER_REACTION(CenterTriangulationVertexFromCellPressure,
                         "CenterTriangulation::VertexFromCellPressure",
                         "VertexFromCellPressurecenterTriangulation")

///
/// The same force, ramped linearly from zero to full strength over deltaT -
/// inflating a tissue gradually rather than kicking it at t = 0.
///
///   CenterTriangulation::VertexFromCellPressureLinear 3 1 1
///     K_force, normalizeVolumeFlag, deltaT
///     center_index
///
class CenterTriangulationVertexFromCellPressureLinear : public Reaction {
public:
  CenterTriangulationVertexFromCellPressureLinear(const ParameterList &p,
                                                  const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "CenterTriangulation::VertexFromCellPressureLinear: level 0 = the "
          "start of the centre-triangulation cell variables.");
    configure("CenterTriangulation::VertexFromCellPressureLinear", p, i, 3,
              {1}, {"K_force", "f_V_norm", "deltaT"});
    if (parameter(1) != 0.0 && parameter(1) != 1.0)
      throw std::runtime_error(
          "CenterTriangulation::VertexFromCellPressureLinear: "
          "normalizeVolumeFlag must be 0 or 1.");
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    centerTriangulationPressure(T, cellData, vertexData, vertexDerivs,
                                variableIndex(0, 0), parameter(1) == 1.0,
                                id().c_str(),
                                [&](size_t) { return timeFactor_ * parameter(0); });
  }
  void update(Tissue &, Matrix &, Matrix &, Matrix &, double h) override {
    if (timeFactor_ < 1.0)
      timeFactor_ += h / parameter(2);
    if (timeFactor_ >= 1.0)
      timeFactor_ = 1.0;
  }

private:
  // Legacy keeps the elapsed time in a function-level `static`, shared by
  // every instance of the class; it only ever feeds a disabled branch, so
  // this keeps the ramp per-instance and drops the static.
  double timeFactor_ = 0.0;
};
TISSUE_REGISTER_REACTION(CenterTriangulationVertexFromCellPressureLinear,
                         "CenterTriangulation::VertexFromCellPressureLinear",
                         "VertexFromCellPressurecenterTriangulationLinear")

///
/// Pressure on a triangular face, pushing all three vertices along the
/// face normal. Ported from legacy mechanical.cc (Pressure3D::Triangular).
///
///   Pressure3D::Triangular 3 0
///     k_force, areaFlag (0 = no area weighting, 1 = weighted), Zplane
///   Pressure3D::Triangular 6 0
///     k_force, areaFlag (2 or 3), Zplane, V0, Vfactor, Pfactor
///
/// The six-parameter form was meant to raise the pressure as the enclosed
/// volume fell, but the line accumulating that volume is commented out in
/// legacy, so the volume term is identically zero and the form reduces to a
/// constant pressure scaled by (Vfactor - Pfactor) / (Vfactor - 1). V0 is
/// then unused, and Zplane is unused in every form. Reproduced as written:
/// these are the coefficients the published runs used.
class Pressure3DTriangular : public Reaction {
public:
  Pressure3DTriangular(const ParameterList &p, const IndexLevels &i) {
    if (!i.empty())
      throw std::runtime_error("Pressure3D::Triangular: uses no indices.");
    if (p.size() != 3 && p.size() != 6)
      throw std::runtime_error(
          "Pressure3D::Triangular: three parameters (k_force, areaFlag 0/1, "
          "Zplane) or six (k_force, areaFlag 2/3, Zplane, V0, Vfactor, "
          "Pfactor).");
    std::vector<std::string> ids{"k_force", "areaFlag", "Zplane"};
    if (p.size() == 6) {
      ids.push_back("V0");
      ids.push_back("Vfactor");
      ids.push_back("Pfactor");
    }
    configure("Pressure3D::Triangular", p, i, p.size(), {}, std::move(ids));
    const double flag = parameter(1);
    if (flag != 0.0 && flag != 1.0 && flag != 2.0 && flag != 3.0)
      throw std::runtime_error(
          "Pressure3D::Triangular: areaFlag must be 0 or 2 (no area "
          "weighting) or 1 or 3 (area weighted).");
    if ((flag >= 2.0) != (p.size() == 6))
      throw std::runtime_error(
          "Pressure3D::Triangular: areaFlag 2 and 3 need the six-parameter "
          "form, 0 and 1 the three-parameter one.");
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "Pressure3D::Triangular: only implemented for three dimensions.");
    const double flag = parameter(1);
    const bool areaWeighted = flag == 1.0 || flag == 3.0;
    // The volume term is dead in legacy (see above), so the six-parameter
    // pressure is this constant.
    const double pressure =
        numParameter() == 6
            ? (parameter(4) - parameter(5)) * parameter(0) / (parameter(4) - 1.0)
            : parameter(0);

    for (size_t n = 0; n < T.numCell(); ++n) {
      const CellTopo &cell = T.cell(n);
      if (cell.numVertex() != 3)
        throw std::runtime_error(
            "Pressure3D::Triangular: only implemented for triangular cells; "
            "cell " + std::to_string(n) + " has " +
            std::to_string(cell.numVertex()) + " vertices.");
      const size_t v0 = cell.vertices[0];
      const size_t v1 = cell.vertices[1];
      const size_t v2 = cell.vertices[2];

      // Normal as legacy builds it: (v1 - v0) x (v2 - v1).
      double e0[3], e1[3], normal[3];
      for (size_t d = 0; d < 3; ++d) {
        e0[d] = vertexData[v1][d] - vertexData[v0][d];
        e1[d] = vertexData[v2][d] - vertexData[v1][d];
      }
      normal[0] = e0[1] * e1[2] - e0[2] * e1[1];
      normal[1] = e0[2] * e1[0] - e0[0] * e1[2];
      normal[2] = e0[0] * e1[1] - e0[1] * e1[0];
      double norm = 0.0;
      for (size_t d = 0; d < 3; ++d)
        norm += normal[d] * normal[d];
      // Legacy skips the normalisation when the squared norm is already
      // exactly 1.0, which is the same result either way.
      if (norm != 1.0) {
        if (!(norm > 0.0))
          throw std::runtime_error(
              "Pressure3D::Triangular: degenerate triangle in cell " +
              std::to_string(n) + " (no well-defined normal).");
        const double inv = 1.0 / std::sqrt(norm);
        for (size_t d = 0; d < 3; ++d)
          normal[d] *= inv;
      }

      double coeff = pressure;
      if (areaWeighted)
        coeff *= T.cellVolume(n, vertexData) / cell.numVertex();
      // Each of the three vertices takes the whole coefficient.
      for (size_t k = 0; k < 3; ++k)
        for (size_t d = 0; d < 3; ++d)
          vertexDerivs[cell.vertices[k]][d] += coeff * normal[d];
    }
  }
};
TISSUE_REGISTER_REACTION(Pressure3DTriangular, "Pressure3D::Triangular",
                         "VertexFromCellPlaneTriangular")

} // namespace
} // namespace tissue
