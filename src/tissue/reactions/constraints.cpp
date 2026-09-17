//
// Vertex constraint reactions (ported from legacy adhocReaction.cc). Placed
// after the force reactions in the model file, they cancel accumulated
// derivatives, clamping vertices — e.g. anchoring the basal end of an organ
// (a seedling fixed at the soil line).
//
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Fully clamps the listed vertices (all coordinates).
class VertexNoUpdateFromIndex : public Reaction {
public:
  VertexNoUpdateFromIndex(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].empty())
      throw std::runtime_error(
          "VertexNoUpdateFromIndex: vertex indices in first level.");
    configure("VertexNoUpdateFromIndex", p, i, 0, {kAnyCount}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dimension = vertexData.cols();
    for (size_t i = 0; i < numVariableIndex(0); ++i) {
      size_t v = variableIndex(0, i);
      for (size_t d = 0; d < dimension; ++d)
        vertexDerivs[v][d] = 0.0;
    }
  }
};
TISSUE_REGISTER_REACTION(VertexNoUpdateFromIndex, "VertexNoUpdateFromIndex")

} // namespace
} // namespace tissue
