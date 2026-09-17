//
// Transport reactions between neighboring cells. Ported from legacy
// transport.cc.
//
#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Plain concentration-difference diffusion across shared (non-boundary)
// walls; no geometric factors ("Simple").
class DiffusionSimple : public Reaction {
public:
  DiffusionSimple(const ParameterList &p, const IndexLevels &i) {
    configure("DiffusionSimple", p, i, 1, {1}, {"p_0"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const double p0 = parameter(0);
    // Legacy iterates cell-wall pairs with an i<neighbor guard; iterating
    // interior walls once is equivalent (each interior wall borders exactly
    // the two cells) and parallel-friendly.
    parallelScatter1(T.numWall(), cellDerivs,
                     [&](size_t b, size_t e, Matrix &out) {
                       for (size_t w = b; w < e; ++w) {
                         const Wall &wall = T.wall(w);
                         if (Tissue::isBackground(wall.cell1) ||
                             Tissue::isBackground(wall.cell2))
                           continue;
                         size_t i = std::min(wall.cell1, wall.cell2);
                         size_t neigh = std::max(wall.cell1, wall.cell2);
                         double flux = p0 * (cellData[i][aI] - cellData[neigh][aI]);
                         out[i][aI] -= flux;
                         out[neigh][aI] += flux;
                       }
                     });
  }
  void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &wallData,
                     Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                     Matrix &vertexDerivs, Matrix &, Matrix &,
                     Matrix &) override {
    // Legacy adds zero noise for diffusion (documented deterministic).
    derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs,
           vertexDerivs);
  }
};
TISSUE_REGISTER_REACTION(DiffusionSimple, "DiffusionSimple", "Diffusion::Simple")

} // namespace
} // namespace tissue
