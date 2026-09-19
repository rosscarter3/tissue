//
// Force::VectorLinear, ported from legacy force.cc.
//
// A constant force applied to a named list of vertices, ramped linearly from
// zero to full strength over deltaT - pulling on a template's edge without
// kicking it at t = 0.
//
#include <stdexcept>
#include <string>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

///
///   Force::VectorLinear 2/3/4 1 N
///     F_x [F_y [F_z]] deltaT
///     vertex indices
///
class ForceVectorLinear : public Reaction {
public:
  ForceVectorLinear(const ParameterList &p, const IndexLevels &i) {
    if (p.size() < 2 || p.size() > 4)
      throw std::runtime_error(
          "Force::VectorLinear: a force vector in one (F_x), two (F_x F_y) "
          "or three (F_x F_y F_z) dimensions, plus the deltaT over which it "
          "ramps up from zero.");
    if (i.size() != 1 || i[0].empty())
      throw std::runtime_error(
          "Force::VectorLinear: level 0 lists the vertices to pull on.");
    std::vector<std::string> ids;
    for (size_t k = 0; k + 1 < p.size(); ++k)
      ids.push_back(std::string("F_") + static_cast<char>('x' + k));
    ids.push_back("deltaT");
    configure("Force::VectorLinear", p, i, p.size(), {kAnyCount},
              std::move(ids));
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dimension = vertexData.cols();
    for (size_t k = 0; k < numVariableIndex(0); ++k) {
      const size_t v = variableIndex(0, k);
      if (v >= vertexData.rows())
        throw std::runtime_error(
            "Force::VectorLinear: vertex index " + std::to_string(v) +
            " is outside the tissue.");
      // Legacy's guard is `numParameter() > d`, which counts deltaT as a
      // force component: with two parameters (F_x, deltaT) in a 2D tissue it
      // applies deltaT as F_y, and with three in a 3D tissue it applies
      // deltaT as F_z. Only the form whose component count matches the
      // tissue's dimension escapes it. Reproduced, because the published
      // runs were made with it and the affected models were tuned against
      // whatever force they actually got.
      for (size_t d = 0; d < dimension; ++d)
        if (numParameter() > d)
          vertexDerivs[v][d] += timeFactor_ * parameter(d);
    }
  }

  void update(Tissue &, Matrix &, Matrix &, Matrix &, double h) override {
    if (timeFactor_ < 1.0)
      timeFactor_ += h / parameter(numParameter() - 1);
    if (timeFactor_ > 1.0)
      timeFactor_ = 1.0;
  }

private:
  double timeFactor_ = 0.0;
};
TISSUE_REGISTER_REACTION(ForceVectorLinear, "Force::VectorLinear",
                         "VertexFromForceLinear")

} // namespace
} // namespace tissue
