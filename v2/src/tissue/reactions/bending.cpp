//
// Wall-chain bending stiffness: cell walls are stiff plates, so chains of
// collinear walls (e.g. the longitudinal wall files of a tissue) resist
// kinking. At every vertex joined by exactly two chain-flagged walls, a
// straightening force pulls the vertex toward the midpoint of its two chain
// neighbors (equal and opposite halves on the neighbors, momentum
// conserving):
//   F_v = k_bend (0.5 (x_a + x_b) - x_v),  F_a = F_b = -F_v/2
// Short-wavelength wrinkles are strongly penalized while smooth organ-scale
// curvature is barely affected. (New in v2; used by the apical-hook model.)
//
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class WallBendingChain : public Reaction {
public:
  WallBendingChain(const ParameterList &p, const IndexLevels &i) {
    configure("WallMechanics::BendingChain", p, i, 1, {1}, {"k_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t flagIndex = variableIndex(0, 0);
    const double kBend = parameter(0);
    const size_t dimension = vertexData.cols();
    parallelScatter1(
        T.numVertex(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t v = b; v < e; ++v) {
            size_t chainWalls[2];
            size_t found = 0;
            for (size_t w : T.vertex(v).walls) {
              if (wallData[w][flagIndex] != 0.0) {
                if (found < 2)
                  chainWalls[found] = w;
                ++found;
              }
            }
            if (found != 2)
              continue; // chain ends or junctions carry no bending force
            size_t a = T.wall(chainWalls[0]).otherVertex(v);
            size_t c = T.wall(chainWalls[1]).otherVertex(v);
            for (size_t d = 0; d < dimension; ++d) {
              double f = kBend * (0.5 * (vertexData[a][d] + vertexData[c][d]) -
                                  vertexData[v][d]);
              out[v][d] += f;
              out[a][d] -= 0.5 * f;
              out[c][d] -= 0.5 * f;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(WallBendingChain, "WallMechanics::BendingChain")

} // namespace
} // namespace tissue
