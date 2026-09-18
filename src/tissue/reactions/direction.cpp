//
// Direction-update reactions ported from legacy directionReaction.cc: a cell
// carries a "target" direction computed elsewhere (a principal stress or
// strain axis, say) and a microtubule direction that follows it, either
// through the derivative or as an explicit step between solver steps.
//
// A direction occupies `dimension` consecutive cell variables starting at the
// given index, so "index 6" in a 3D model means cell variables 6, 7 and 8.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Legacy's myMath::pi() is this truncated literal rather than M_PI, and the
// angle arithmetic below is sensitive to it.
constexpr double kLegacyPi = 3.14159265;

// Folds an angle into (-pi/2, pi/2]: these are axes, not vectors, so a
// direction and its negation are the same thing.
double foldToAxis(double a) {
  while (a > 0.5 * kLegacyPi || a <= -0.5 * kLegacyPi) {
    if (a > 0.5 * kLegacyPi)
      a -= kLegacyPi;
    if (a <= -0.5 * kLegacyPi)
      a += kLegacyPi;
  }
  return a;
}

// Normalises the direction at `index`, in place. Legacy does not guard
// against a zero vector here, and neither does this: a zero direction is a
// broken model, and silently substituting one would hide that.
void normalizeAt(Matrix &cellData, size_t cell, size_t index, size_t dim) {
  double norm = 0.0;
  for (size_t d = 0; d < dim; ++d)
    norm += cellData[cell][index + d] * cellData[cell][index + d];
  norm = 1.0 / std::sqrt(norm);
  for (size_t d = 0; d < dim; ++d)
    cellData[cell][index + d] *= norm;
}

// Flips `out` if it points away from `in`, so the step that follows moves the
// axis the short way round rather than through a half turn.
void alignAxis(Matrix &cellData, size_t cell, size_t in, size_t out,
               size_t dim) {
  double inner = 0.0;
  for (size_t d = 0; d < dim; ++d)
    inner += cellData[cell][out + d] * cellData[cell][in + d];
  if (inner < 0.0)
    for (size_t d = 0; d < dim; ++d)
      cellData[cell][out + d] *= -1.0;
}

///
/// Rotates a 2D direction towards a target through the derivative, at a speed
/// that saturates with the angle between them.
///
///   ContinousMTDirection 1 2 1 1
///     k_rate
///     target_index
///     MT_index
///
class ContinousMTDirection : public Reaction {
public:
  ContinousMTDirection(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "ContinousMTDirection: level 0 = the target direction to follow; "
          "level 1 = the direction to update.");
    configure("ContinousMTDirection", p, i, 1, {1, 1}, {"k_rate"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "ContinousMTDirection: only implemented for two dimensions.");
    const size_t target = variableIndex(0, 0);
    const size_t real = variableIndex(1, 0);
    const double kRate = parameter(0);

    for (size_t n = 0; n < T.numCell(); ++n) {
      const double x = cellData[n][real + 0];
      const double y = cellData[n][real + 1];
      const double sigma = foldToAxis(std::atan2(y, x));
      const double dsigma =
          foldToAxis(std::atan2(cellData[n][target + 1],
                                cellData[n][target + 0]));

      double angle = dsigma - sigma;
      while (angle > kLegacyPi || angle <= -kLegacyPi) {
        if (angle > kLegacyPi)
          angle -= 2.0 * kLegacyPi;
        if (angle <= -kLegacyPi)
          angle += 2.0 * kLegacyPi;
      }

      // Saturating speed: proportional to the angle when small, tending to
      // k_rate when large.
      double speed = kRate * std::abs(angle) / (0.25 * kLegacyPi + std::abs(angle));
      speed *= (angle >= 0.0) ? 1.0 : -1.0;

      // A rotation of (x, y), so the direction turns without changing length.
      cellDerivs[n][real + 0] += -y * speed;
      cellDerivs[n][real + 1] += x * speed;
    }
  }
};
TISSUE_REGISTER_REACTION(ContinousMTDirection, "ContinousMTDirection")

///
/// The 3D counterpart: moves the direction towards the target along the
/// straight line between them rather than by rotating.
///
///   ContinousMTDirection3d 1 2 1 1
///     k_rate
///     target_index
///     MT_index
///
/// Note this one normalises both directions *in cellData* from inside
/// derivs(), so it rewrites state during a derivative evaluation rather than
/// only contributing a rate; see NOTES.md on why that shows at the t=0 print.
class ContinousMTDirection3d : public Reaction {
public:
  ContinousMTDirection3d(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "ContinousMTDirection3d: level 0 = the target direction to follow; "
          "level 1 = the direction to update.");
    configure("ContinousMTDirection3d", p, i, 1, {1, 1}, {"k_rate"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t target = variableIndex(0, 0);
    const size_t real = variableIndex(1, 0);
    const double kRate = parameter(0);

    for (size_t n = 0; n < T.numCell(); ++n) {
      // A zero direction is replaced by the x axis rather than dividing by
      // zero - legacy's one guard of this kind, kept.
      for (size_t index : {real, target}) {
        double norm = 0.0;
        for (size_t d = 0; d < 3; ++d)
          norm += cellData[n][index + d] * cellData[n][index + d];
        norm = std::sqrt(norm);
        if (norm == 0.0) {
          cellData[n][index + 0] = 1.0;
          norm = 1.0;
        }
        for (size_t d = 0; d < 3; ++d)
          cellData[n][index + d] /= norm;
      }

      alignAxis(cellData, n, target, real, 3);
      for (size_t d = 0; d < 3; ++d)
        cellDerivs[n][real + d] +=
            kRate * (cellData[n][target + d] - cellData[n][real + d]);
    }
  }
};
TISSUE_REGISTER_REACTION(ContinousMTDirection3d, "ContinousMTDirection3d")

///
/// Steps the direction towards the target between solver steps and
/// renormalises. Initiated to the target.
///
///   UpdateMTDirection 1 2 1 1
///     k_rate
///     target_index
///     MT_index
///
/// Unlike UpdateMTDirectionEquilibrium below, this one does *not* align the
/// two axes first, so a target pointing the other way drags the direction
/// through zero length rather than taking the short way round.
class UpdateMTDirection : public Reaction {
public:
  UpdateMTDirection(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "UpdateMTDirection: level 0 = the target direction to follow; "
          "level 1 = the direction to update.");
    configure("UpdateMTDirection", p, i, 1, {1, 1}, {"k_rate"});
  }
  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const size_t dim = vertexData.cols();
    for (size_t n = 0; n < cellData.rows(); ++n)
      for (size_t d = 0; d < dim; ++d)
        cellData[n][variableIndex(1, 0) + d] =
            cellData[n][variableIndex(0, 0) + d];
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
              double h) override {
    if (parameter(0) == 0.0)
      return;
    const size_t dim = vertexData.cols();
    const size_t in = variableIndex(0, 0);
    const size_t out = variableIndex(1, 0);
    for (size_t n = 0; n < cellData.rows(); ++n) {
      for (size_t d = 0; d < dim; ++d)
        cellData[n][out + d] +=
            parameter(0) * h * (cellData[n][in + d] - cellData[n][out + d]);
      normalizeAt(cellData, n, out, dim);
    }
  }
};
TISSUE_REGISTER_REACTION(UpdateMTDirection, "UpdateMTDirection")

///
/// As UpdateMTDirection, but the step is taken only in cells that are close
/// enough to mechanical equilibrium, and optionally only where two stress
/// measures disagree by more than a threshold - the direction of a principal
/// axis is meaningless while the two are nearly equal.
///
///   UpdateMTDirectionEquilibrium 1/2/3 2/3/4 1 1 [1] [2]
///     k_rate [velocity_threshold] [stressDiff_threshold]
///     target_index
///     MT_index
///     [velocity_index]
///     [stress1_index stress2_index]
///
class UpdateMTDirectionEquilibrium : public Reaction {
public:
  UpdateMTDirectionEquilibrium(const ParameterList &p, const IndexLevels &i) {
    if (p.size() < 1 || p.size() > 3)
      throw std::runtime_error(
          "UpdateMTDirectionEquilibrium: k_rate, plus an optional velocity "
          "threshold and an optional stress-difference threshold.");
    if (i.size() < 2 || i.size() > 4 || i[0].size() != 1 ||
        i[1].size() != 1 || (i.size() > 2 && i[2].size() != 1) ||
        (i.size() == 4 && i[3].size() != 2))
      throw std::runtime_error(
          "UpdateMTDirectionEquilibrium: level 0 = the target direction; "
          "level 1 = the direction to update; optional level 2 = the stored "
          "velocity; optional level 3 = the two stresses to compare.");
    std::vector<size_t> shape{1, 1};
    if (i.size() > 2)
      shape.push_back(1);
    if (i.size() == 4)
      shape.push_back(2);
    std::vector<std::string> names{"k_rate"};
    if (p.size() > 1)
      names.push_back("velocity_threshold");
    if (p.size() > 2)
      names.push_back("stressDiff_threshold");
    configure("UpdateMTDirectionEquilibrium", p, i, p.size(), shape,
              std::move(names));
  }
  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const size_t dim = vertexData.cols();
    for (size_t n = 0; n < cellData.rows(); ++n)
      for (size_t d = 0; d < dim; ++d)
        cellData[n][variableIndex(1, 0) + d] =
            cellData[n][variableIndex(0, 0) + d];
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
              double h) override {
    if (parameter(0) == 0.0)
      return;
    const size_t dim = vertexData.cols();
    const size_t in = variableIndex(0, 0);
    const size_t out = variableIndex(1, 0);
    const size_t np = numParameter();

    for (size_t n = 0; n < cellData.rows(); ++n) {
      // The gates are keyed off the *parameter* count, not the index count,
      // exactly as legacy: with three parameters both conditions must hold.
      bool doUpdate = (np == 1);
      if (np == 2)
        doUpdate = cellData[n][variableIndex(2, 0)] < parameter(1);
      else if (np == 3) {
        const double s1 = cellData[n][variableIndex(3, 0)];
        const double s2 = cellData[n][variableIndex(3, 1)];
        doUpdate = cellData[n][variableIndex(2, 0)] < parameter(1) &&
                   std::abs((s1 - s2) / s1) > parameter(2);
      }
      if (!doUpdate)
        continue;

      alignAxis(cellData, n, in, out, dim);
      for (size_t d = 0; d < dim; ++d)
        cellData[n][out + d] +=
            parameter(0) * h * (cellData[n][in + d] - cellData[n][out + d]);
      normalizeAt(cellData, n, out, dim);
    }
  }
};
TISSUE_REGISTER_REACTION(UpdateMTDirectionEquilibrium,
                         "UpdateMTDirectionEquilibrium")

///
/// As UpdateMTDirection, but the rate is scaled by a Hill function of a cell
/// concentration (in practice a stress or strain anisotropy), so the
/// direction only follows where the target is well defined.
///
///   UpdateMTDirectionConcenHill 3 3 1 1 1
///     k_rate, k_Hill, n_Hill
///     target_index
///     MT_index
///     concentration_index
///
class UpdateMTDirectionConcenHill : public Reaction {
public:
  UpdateMTDirectionConcenHill(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 3 || i[0].size() != 1 || i[1].size() != 1 ||
        i[2].size() != 1)
      throw std::runtime_error(
          "UpdateMTDirectionConcenHill: level 0 = the target direction; "
          "level 1 = the direction to update; level 2 = the concentration "
          "scaling the rate.");
    configure("UpdateMTDirectionConcenHill", p, i, 3, {1, 1, 1},
              {"k_rate", "k_Hill", "n_Hill"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &) override {
    if (parameter(0) == 0.0)
      return;
    const size_t dim = vertexData.cols();
    const size_t in = variableIndex(0, 0);
    const size_t out = variableIndex(1, 0);
    const size_t conc = variableIndex(2, 0);
    const double kh = parameter(1);
    const double nh = parameter(2);

    for (size_t n = 0; n < cellData.rows(); ++n) {
      const double c = cellData[n][conc];
      const double cN = std::pow(c, nh);
      const double concFactor = cN / (std::pow((1.0 - c) * kh, nh) + cN);

      // Legacy's alignment test is hard-coded to three components even though
      // everything around it loops over `dimension`, so in a 2D model it also
      // reads whatever cell variable follows each direction. Kept: changing it
      // would change 2D results, and no 2D model is known to use this.
      double inner = 0.0;
      for (size_t d = 0; d < 3; ++d)
        inner += cellData[n][in + d] * cellData[n][out + d];
      if (inner < 0.0)
        for (size_t d = 0; d < dim; ++d)
          cellData[n][out + d] *= -1.0;

      // No h and no cellDerivs: this writes the direction straight into
      // cellData once per *derivative evaluation*, so an adaptive solver
      // applies it as many times per step as it evaluates.
      for (size_t d = 0; d < dim; ++d)
        cellData[n][out + d] +=
            parameter(0) * concFactor * (cellData[n][in + d] -
                                         cellData[n][out + d]);
      normalizeAt(cellData, n, out, dim);
    }
  }
};
TISSUE_REGISTER_REACTION(UpdateMTDirectionConcenHill,
                         "UpdateMTDirectionConcenHill")

///
/// Rotates a direction at a constant rate in the xy plane, for testing what a
/// model does under a turning field.
///
///   RotatingDirection 1 1 1
///     k_rate
///     MT_index
///
class RotatingDirection : public Reaction {
public:
  RotatingDirection(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "RotatingDirection: level 0 = the direction to rotate.");
    configure("RotatingDirection", p, i, 1, {1}, {"k_rate"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t xIndex = variableIndex(0, 0);
    const size_t yIndex = variableIndex(0, 0) + 1;
    const double k = parameter(0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      const double x = cellData[n][xIndex];
      const double y = cellData[n][yIndex];
      cellDerivs[n][xIndex] -= k * y;
      cellDerivs[n][yIndex] += k * x;
    }
  }
};
TISSUE_REGISTER_REACTION(RotatingDirection, "RotatingDirection")

} // namespace
} // namespace tissue
