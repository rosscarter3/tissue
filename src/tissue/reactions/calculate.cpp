//
// Diagnostic reactions ported from legacy calculate.cc: these compute nothing
// the dynamics depend on directly, but write quantities other reactions read
// as gates - most importantly the per-cell vertex velocity that
// UpdateMTDirectionEquilibrium and FiberModel use to decide whether the
// tissue is close enough to mechanical equilibrium to update.
//
// Every one of them writes into cellData from derivs() rather than
// contributing to cellDerivs, so they run once per derivative evaluation and
// differ from legacy at the t=0 print only; see NOTES.md.
//
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Normalises a 3-vector at `index` in place, leaving it alone if it is zero.
// Returns the length it had.
double normalizeInPlace(Matrix &cellData, size_t cell, size_t index) {
  double norm = 0.0;
  for (size_t d = 0; d < 3; ++d)
    norm += cellData[cell][index + d] * cellData[cell][index + d];
  norm = std::sqrt(norm);
  if (norm > 0.0)
    for (size_t d = 0; d < 3; ++d)
      cellData[cell][index + d] /= norm;
  return norm;
}

///
/// |cos| of the angle between two 3D directions stored in cellData, which is
/// how a model reports the agreement between (say) a stress axis and the
/// microtubule direction following it: 1 aligned, 0 perpendicular.
///
///   Calculate::AngleVectors 0 2 2 1
///     vector1_index vector2_index
///     store_index
///
/// Both input vectors are normalised in place, so this rewrites them.
class CalculateAngleVectors : public Reaction {
public:
  CalculateAngleVectors(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
      throw std::runtime_error(
          "Calculate::AngleVectors: level 0 = the two vector start indices; "
          "level 1 = where to store |cos| of the angle between them.");
    configure("Calculate::AngleVectors", p, i, 0, {2, 1}, {});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &, Matrix &,
              Matrix &, Matrix &) override {
    const size_t a = variableIndex(0, 0);
    const size_t b = variableIndex(0, 1);
    const size_t store = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      normalizeInPlace(cellData, n, a);
      normalizeInPlace(cellData, n, b);
      double dot = 0.0;
      for (size_t d = 0; d < 3; ++d)
        dot += cellData[n][a + d] * cellData[n][b + d];
      cellData[n][store] = std::abs(dot);
    }
  }
};
TISSUE_REGISTER_REACTION(CalculateAngleVectors, "Calculate::AngleVectors",
                         "CalculateAngleVectors")

///
/// Elevation of a 3D direction above the xy plane, in radians.
///
///   Calculate::AngleVectorXYplane 0 2 1 1
///     vector_index
///     store_index
///
/// The class documentation calls this "abs(cos(...)) of angle between a
/// vector and XY plane", but the code stores `atan(z / |xy|)` - an angle in
/// radians, not a cosine. The code is what models were written against.
class CalculateAngleVectorXYplane : public Reaction {
public:
  CalculateAngleVectorXYplane(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "Calculate::AngleVectorXYplane: level 0 = the vector start index; "
          "level 1 = where to store the angle.");
    configure("Calculate::AngleVectorXYplane", p, i, 0, {1, 1}, {});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &, Matrix &,
              Matrix &, Matrix &) override {
    constexpr double kPi = 3.14159265; // legacy's literal, local to this class
    const size_t v = variableIndex(0, 0);
    const size_t store = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      const double xy = std::sqrt(cellData[n][v] * cellData[n][v] +
                                  cellData[n][v + 1] * cellData[n][v + 1]);
      // A vector along z has no xy projection to take the angle from, so it
      // is called vertical rather than dividing by (almost) zero.
      cellData[n][store] =
          xy < 1e-8 ? 0.5 * kPi : std::atan(cellData[n][v + 2] / xy);
    }
  }
};
TISSUE_REGISTER_REACTION(CalculateAngleVectorXYplane,
                         "Calculate::AngleVectorXYplane",
                         "CalculateAngleVectorXYplane")

///
/// Angle in degrees between a 3D direction and the x axis, folded to [0, 180)
/// by the sign of the y component.
///
///   Calculate::AngleVector 1 2 1 1
///     axis (only 0, the x axis, is implemented)
///     vector_index
///     store_index
///
/// Legacy accepts axis 1 and 2 in its range check and then exits immediately
/// with "The code should be modified for 3d", so only axis 0 has ever run.
/// This rejects 1 and 2 at construction instead, which fails with a message
/// naming the model rule rather than exiting with status 0 mid-run.
class CalculateAngleVector : public Reaction {
public:
  CalculateAngleVector(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "Calculate::AngleVector: level 0 = the vector start index; level 1 "
          "= where to store the angle in degrees.");
    configure("Calculate::AngleVector", p, i, 1, {1, 1}, {"axis"});
    if (parameter(0) != 0)
      throw std::runtime_error(
          "Calculate::AngleVector: only axis 0 (x) is implemented. Legacy "
          "accepts 1 and 2 and then exits without running.");
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &, Matrix &,
              Matrix &, Matrix &) override {
    constexpr double kPi = 3.141592; // legacy spells pi with 7 digits here
    const size_t v = variableIndex(0, 0);
    const size_t store = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      // Unguarded, as in legacy: a zero direction gives a NaN angle rather
      // than being silently treated as pointing along x.
      double norm = 0.0;
      for (size_t d = 0; d < 3; ++d)
        norm += cellData[n][v + d] * cellData[n][v + d];
      norm = std::sqrt(norm);
      for (size_t d = 0; d < 3; ++d)
        cellData[n][v + d] /= norm;

      double teta = 180.0 * std::acos(cellData[n][v]) / kPi;
      if (cellData[n][v + 1] < 0.0)
        teta = 180.0 - teta;
      cellData[n][store] = teta;
    }
  }
};
TISSUE_REGISTER_REACTION(CalculateAngleVector, "Calculate::AngleVector",
                         "AngleVector")

///
/// Mean speed of a cell's vertices, as a measure of how far the tissue is
/// from mechanical equilibrium.
///
///   Calculate::VertexVelocity 0 1 1
///     store_index
///
/// Reads vertexDerivs directly, so it must appear *after* every reaction that
/// moves vertices, or it sees a partial force. (It was once called
/// Calculate::MaxVelocity / maxVelocity; legacy now refuses those names, and
/// so does this build, since the quantity is a mean and never was a maximum.)
class CalculateVertexVelocity : public Reaction {
public:
  CalculateVertexVelocity(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "Calculate::VertexVelocity: level 0 = where to store the mean "
          "vertex speed.");
    configure("Calculate::VertexVelocity", p, i, 0, {1}, {});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t store = variableIndex(0, 0);
    const size_t dim = vertexData.cols();
    for (size_t n = 0; n < T.numCell(); ++n) {
      const CellTopo &cell = T.cell(n);
      const size_t nv = cell.numVertex();
      double total = 0.0;
      for (size_t k = 0; k < nv; ++k) {
        const size_t v = cell.vertices[k];
        double speed = 0.0;
        for (size_t d = 0; d < dim; ++d)
          speed += vertexDerivs[v][d] * vertexDerivs[v][d];
        total += std::sqrt(speed);
      }
      cellData[n][store] = total / nv;
    }
  }
};
TISSUE_REGISTER_REACTION(CalculateVertexVelocity, "Calculate::VertexVelocity")

///
/// Total distance the tissue's vertices have moved from their starting
/// positions, the change in that since the last evaluation, and the summed
/// vertex speed. 3D only.
///
///   Calculate::TissueVolumeChange 0 1 6
///     cell1 col1  cell2 col2  cell3 col3
///
/// Unusually, the six indices are three (cell, column) *pairs* naming single
/// cellData entries - these are tissue-wide numbers parked wherever the model
/// author had room, not per-cell variables.
///
/// Despite the name none of the three is a volume: the first is a sum of
/// vertex displacements. Kept as the name models use.
class CalculateTissueVolumeChange : public Reaction {
public:
  CalculateTissueVolumeChange(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 6)
      throw std::runtime_error(
          "Calculate::TissueVolumeChange: level 0 = three (cell, column) "
          "pairs naming where to store the change, its increment and the "
          "summed vertex speed.");
    configure("Calculate::TissueVolumeChange", p, i, 0, {6}, {});
  }
  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "Calculate::TissueVolumeChange: only implemented for 3D models.");
    rest_.assign(T.numVertex(), {0.0, 0.0, 0.0});
    for (size_t v = 0; v < T.numVertex(); ++v)
      for (size_t d = 0; d < 3; ++d)
        rest_[v][d] = vertexData[v][d];
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    double change = 0.0;
    double totalDerivative = 0.0;
    for (size_t v = 0; v < T.numVertex(); ++v) {
      double moved = 0.0;
      double speed = 0.0;
      for (size_t d = 0; d < 3; ++d) {
        const double dx = vertexData[v][d] - rest_[v][d];
        moved += dx * dx;
        speed += vertexDerivs[v][d] * vertexDerivs[v][d];
      }
      change += std::sqrt(moved);
      totalDerivative += std::sqrt(speed);
    }
    // The increment is against whatever is in the store right now, so it is
    // the change since the previous *derivative evaluation*, not per step.
    const double delta =
        change - cellData[variableIndex(0, 0)][variableIndex(0, 1)];
    cellData[variableIndex(0, 0)][variableIndex(0, 1)] = change;
    cellData[variableIndex(0, 2)][variableIndex(0, 3)] = delta;
    cellData[variableIndex(0, 4)][variableIndex(0, 5)] = totalDerivative;
  }

private:
  std::vector<std::array<double, 3>> rest_;
};
TISSUE_REGISTER_REACTION(CalculateTissueVolumeChange,
                         "Calculate::TissueVolumeChange",
                         "TemplateVolumeChange")

} // namespace
} // namespace tissue
