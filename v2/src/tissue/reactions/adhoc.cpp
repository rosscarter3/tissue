//
// Ad hoc utility reactions. Ported from legacy adhocReaction.cc.
//
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Translates the whole tissue so the (unweighted) vertex centroid sits at the
// origin; applied at initiation and after every step. Keeps plots centered.
class CenterCOM : public Reaction {
public:
  CenterCOM(const ParameterList &p, const IndexLevels &i) {
    configure("CenterCOM", p, i, 0, {}, {});
  }
  void initiate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &, Matrix &, Matrix &) override {
    update(T, cellData, wallData, vertexData, 0.0);
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &, Matrix &, Matrix &, Matrix &vertexData,
              double) override {
    const size_t dimension = vertexData.cols();
    const size_t n = vertexData.rows();
    double com[3] = {0.0, 0.0, 0.0};
    for (size_t v = 0; v < n; ++v)
      for (size_t d = 0; d < dimension; ++d)
        com[d] += vertexData[v][d];
    for (size_t d = 0; d < dimension; ++d)
      com[d] /= static_cast<double>(n);
    for (size_t v = 0; v < n; ++v)
      for (size_t d = 0; d < dimension; ++d)
        vertexData[v][d] -= com[d];
  }
};
TISSUE_REGISTER_REACTION(CenterCOM, "CenterCOM")

// As CenterCOM, but also translates the center-triangulation centers stored
// in cell variables (starting at the given index).
class CenterCOMcenterTriangulation : public Reaction {
public:
  CenterCOMcenterTriangulation(const ParameterList &p, const IndexLevels &i) {
    configure("CenterCOMcenterTriangulation", p, i, 0, {1}, {});
  }
  void initiate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &, Matrix &, Matrix &) override {
    update(T, cellData, wallData, vertexData, 0.0);
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
              double) override {
    const size_t dimension = vertexData.cols();
    const size_t n = vertexData.rows();
    double com[3] = {0.0, 0.0, 0.0};
    for (size_t v = 0; v < n; ++v)
      for (size_t d = 0; d < dimension; ++d)
        com[d] += vertexData[v][d];
    for (size_t d = 0; d < dimension; ++d)
      com[d] /= static_cast<double>(n);
    for (size_t v = 0; v < n; ++v)
      for (size_t d = 0; d < dimension; ++d)
        vertexData[v][d] -= com[d];
    const size_t startIndex = variableIndex(0, 0);
    for (size_t c = 0; c < cellData.rows(); ++c)
      for (size_t d = 0; d < dimension; ++d)
        cellData[c][startIndex + d] -= com[d];
  }
};
TISSUE_REGISTER_REACTION(CenterCOMcenterTriangulation,
                         "CenterCOMcenterTriangulation")

// Sets every wall's resting length (column 0) to factor times the current
// vertex distance, once before the simulation.
class InitiateWallLength : public Reaction {
public:
  InitiateWallLength(const ParameterList &p, const IndexLevels &i) {
    configure("InitiateWallLength", p, i, 1, {}, {"factor"});
  }
  void initiate(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const double factor = parameter(0);
    for (size_t w = 0; w < T.numWall(); ++w)
      wallData[w][0] = factor * T.wallLengthFromVertices(w, vertexData);
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(InitiateWallLength, "InitiateWallLength")

} // namespace
} // namespace tissue
