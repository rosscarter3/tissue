//
// Growth-force reactions moving vertices directly. Ported from legacy
// growthForce.cc.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Radial outward movement of all vertices:
//   r_pow=1: dx/dt += k * x        (exponential growth)
//   r_pow=0: dx/dt += k * x / |x|  (constant speed; origin vertices fixed)
class GrowthForceRadial : public Reaction {
public:
  GrowthForceRadial(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 || (p[1] != 0.0 && p[1] != 1.0))
      throw std::runtime_error("MoveVertexRadially: uses two parameters "
                               "(k_growth, r_pow in {0,1}).");
    configure("Force::Radial", p, i, 2, {}, {"k_growth", "r_pow"});
  }
  // Radial expansion is imposed motion, not a force: it does not vanish at
  // force balance, so a force-balance solver must integrate it rather than
  // relax it (see Reaction::prescribesVelocity).
  bool prescribesVelocity() const override { return true; }
  void velocityDerivs(Tissue &T, Matrix &cellData, Matrix &wallData,
                      Matrix &vertexData, Matrix &vertexVel) override {
    Matrix ignoredCell, ignoredWall;
    ignoredCell.reshapeLike(cellData);
    ignoredWall.reshapeLike(wallData);
    derivs(T, cellData, wallData, vertexData, ignoredCell, ignoredWall,
           vertexVel);
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const double k = parameter(0);
    const bool constantSpeed = parameter(1) == 0.0;
    const size_t dimension = vertexData.cols();
    parallelFor(vertexDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t v = b; v < e; ++v) {
        double fac = k;
        if (constantSpeed) {
          double r = 0.0;
          for (size_t d = 0; d < dimension; ++d)
            r += vertexData[v][d] * vertexData[v][d];
          r = std::sqrt(r);
          fac = r > 0.0 ? k / r : 0.0;
        }
        for (size_t d = 0; d < dimension; ++d)
          vertexDerivs[v][d] += fac * vertexData[v][d];
      }
    });
  }
};
TISSUE_REGISTER_REACTION(GrowthForceRadial, "GrowthForce::Radial",
                         "MoveVertexRadially")

} // namespace
} // namespace tissue
