//
// Vertex constraint reactions (ported from legacy adhocReaction.cc). Placed
// after the force reactions in the model file, they cancel accumulated
// derivatives, clamping vertices — e.g. anchoring the basal end of an organ
// (a seedling fixed at the soil line).
//
#include <algorithm>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// A vertex is on the tissue boundary when any of its walls faces the
// background. Legacy spells this Vertex::isBoundary.
bool vertexOnBoundary(const Tissue &T, size_t v) {
  for (size_t w : T.vertex(v).walls)
    if (Tissue::isBackground(T.wall(w).cell1) ||
        Tissue::isBackground(T.wall(w).cell2))
      return true;
  return false;
}

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

// Clamps one coordinate of the listed vertices, leaving the others free - a
// vertex allowed to slide along an axis but not across it. Legacy ships three
// near-identical classes, one per axis, so the axis is a template argument
// rather than three copies of the loop.
template <size_t kAxis> class VertexNoUpdateFromIndexHold : public Reaction {
public:
  VertexNoUpdateFromIndexHold(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].empty())
      throw std::runtime_error(std::string(kName) +
                               ": vertex indices in first level.");
    configure(kName, p, i, 0, {kAnyCount}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    // Legacy indexes the axis unguarded, so holding Z in a 2D model reads off
    // the end of the row; this says so instead.
    if (kAxis >= vertexData.cols())
      throw std::runtime_error(std::string(kName) +
                               ": the tissue has too few dimensions for that "
                               "axis.");
    for (size_t i = 0; i < numVariableIndex(0); ++i)
      vertexDerivs[variableIndex(0, i)][kAxis] = 0.0;
  }

private:
  static constexpr const char *kName =
      kAxis == 0   ? "VertexNoUpdateFromIndexHoldX"
      : kAxis == 1 ? "VertexNoUpdateFromIndexHoldY"
                   : "VertexNoUpdateFromIndexHoldZ";
};
using VertexNoUpdateFromIndexHoldX = VertexNoUpdateFromIndexHold<0>;
using VertexNoUpdateFromIndexHoldY = VertexNoUpdateFromIndexHold<1>;
using VertexNoUpdateFromIndexHoldZ = VertexNoUpdateFromIndexHold<2>;
TISSUE_REGISTER_REACTION(VertexNoUpdateFromIndexHoldX,
                         "VertexNoUpdateFromIndexHoldX")
TISSUE_REGISTER_REACTION(VertexNoUpdateFromIndexHoldY,
                         "VertexNoUpdateFromIndexHoldY")
TISSUE_REGISTER_REACTION(VertexNoUpdateFromIndexHoldZ,
                         "VertexNoUpdateFromIndexHoldZ")

// Clamps every vertex whose position along one axis is past a threshold -
// freezing everything below a soil line, say, without listing the vertices.
//
//   VertexNoUpdateFromPosition 2 1 1
//     threshold, direction (+1 = freeze above, -1 = freeze below)
//     position_index
class VertexNoUpdateFromPosition : public Reaction {
public:
  VertexNoUpdateFromPosition(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "VertexNoUpdateFromPosition: the position index in first level.");
    configure("VertexNoUpdateFromPosition", p, i, 2, {1},
              {"threshold", "direction"});
    if (parameter(1) != 1.0 && parameter(1) != -1.0)
      throw std::runtime_error(
          "VertexNoUpdateFromPosition: direction must be 1 (freeze above the "
          "threshold) or -1 (freeze below it).");
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t posIndex = variableIndex(0, 0);
    const size_t dimension = vertexData.cols();
    if (posIndex >= dimension)
      throw std::runtime_error(
          "VertexNoUpdateFromPosition: position index outside the tissue's "
          "dimensions.");
    const bool above = parameter(1) > 0.0;
    for (size_t v = 0; v < T.numVertex(); ++v) {
      const double x = vertexData[v][posIndex];
      if (above ? x > parameter(0) : x < parameter(0))
        for (size_t d = 0; d < dimension; ++d)
          vertexDerivs[v][d] = 0.0;
    }
  }
};
TISSUE_REGISTER_REACTION(VertexNoUpdateFromPosition,
                         "VertexNoUpdateFromPosition")

// Freezes every vertex *except* those on the tissue boundary joining exactly
// two walls - the corner-free boundary vertices. The inverse sense is
// legacy's: the list names what may move.
//
//   VertexNoUpdateFromList 0 0
//
// The list is rebuilt in update(), which runs after the first derivative
// evaluation, so on the first step it is empty and *everything* is frozen.
// Kept: it is what models using this were tuned against, and the effect is
// one step of a held tissue.
class VertexNoUpdateFromList : public Reaction {
public:
  VertexNoUpdateFromList(const ParameterList &p, const IndexLevels &i) {
    if (!i.empty())
      throw std::runtime_error(
          "VertexNoUpdateFromList: uses no variable indices.");
    configure("VertexNoUpdateFromList", p, i, 0, {}, {});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dimension = vertexData.cols();
    for (size_t v = 0; v < T.numVertex(); ++v)
      if (std::find(mobile_.begin(), mobile_.end(), v) == mobile_.end())
        for (size_t d = 0; d < dimension; ++d)
          vertexDerivs[v][d] = 0.0;
  }
  void update(Tissue &T, Matrix &, Matrix &, Matrix &, double) override {
    mobile_.clear();
    for (size_t v = 0; v < T.numVertex(); ++v) {
      const VertexTopo &vertex = T.vertex(v);
      if (vertex.walls.size() != 2)
        continue;
      bool bothOnBoundary = true;
      for (size_t w : vertex.walls)
        if (!Tissue::isBackground(T.wall(w).cell1) &&
            !Tissue::isBackground(T.wall(w).cell2))
          bothOnBoundary = false;
      if (bothOnBoundary)
        mobile_.push_back(v);
    }
  }

private:
  std::vector<size_t> mobile_;
};
TISSUE_REGISTER_REACTION(VertexNoUpdateFromList, "VertexNoUpdateFromList")

// Clamps every boundary vertex, either in all directions (no indices) or only
// along the axes listed - a tissue held at its rim while its interior moves.
//
//   VertexNoUpdateBoundary 0 0
//   VertexNoUpdateBoundary 0 1 N   (axes to freeze)
class VertexNoUpdateBoundary : public Reaction {
public:
  VertexNoUpdateBoundary(const ParameterList &p, const IndexLevels &i) {
    if (i.size() > 1)
      throw std::runtime_error(
          "VertexNoUpdateBoundary: either no indices (every direction is "
          "fixed) or the directions to fix, in first level.");
    configure("VertexNoUpdateBoundary", p, i, 0,
              i.empty() ? std::vector<size_t>{}
                        : std::vector<size_t>{kAnyCount},
              {});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dimension = vertexData.cols();
    const bool someAxes = numVariableIndexLevel() > 0;
    if (someAxes)
      for (size_t k = 0; k < numVariableIndex(0); ++k)
        if (variableIndex(0, k) >= dimension)
          throw std::runtime_error(
              "VertexNoUpdateBoundary: direction index outside the tissue's "
              "dimensions.");
    for (size_t v = 0; v < T.numVertex(); ++v) {
      if (!vertexOnBoundary(T, v))
        continue;
      if (someAxes)
        for (size_t k = 0; k < numVariableIndex(0); ++k)
          vertexDerivs[v][variableIndex(0, k)] = 0.0;
      else
        for (size_t d = 0; d < dimension; ++d)
          vertexDerivs[v][d] = 0.0;
    }
  }
};
TISSUE_REGISTER_REACTION(VertexNoUpdateBoundary, "VertexNoUpdateBoundary")

} // namespace
} // namespace tissue
