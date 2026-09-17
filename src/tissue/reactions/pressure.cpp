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
} // namespace tissue
