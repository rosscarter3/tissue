//
// Wall-chain bending stiffness: cell walls are stiff plates, so chains of
// collinear walls (e.g. the longitudinal wall files of a tissue) resist
// kinking. At every vertex joined by exactly two chain-flagged walls, a
// straightening force pulls the vertex toward the midpoint of its two chain
// neighbors (equal and opposite halves on the neighbors, momentum
// conserving):
//   F_v = k_bend (0.5 (x_a + x_b) - x_v),  F_a = F_b = -F_v/2
// This is a discrete Laplacian stencil, not a curvature operator. Its energy
// is E = sum_v (k/2)|s_v|^2 with s_v = 0.5(x_a + x_b) - x_v (the forces above
// are exactly -dE/dx once each vertex's own and its neighbours' stencils are
// summed, so it is conservative). For a sinusoid of amplitude A and
// wavelength L sampled at spacing h,
//
//   E / length = (k/4) A^2 (1 - cos(2 pi h / L))^2 / h  ->  k pi^4 A^2 h^3 / L^4
//
// verified numerically against the discrete sum (ratio 7.45, 7.86, 7.96, 7.99
// over successive halvings of h, converging on 8).
//
// Two consequences. The wavelength dependence is L^-4, so short-wavelength
// wrinkles really are strongly penalized while smooth organ-scale curvature is
// barely affected - the reaction does the job it claims to. But the effective
// stiffness carries an explicit h^3, so **k_bend is tied to the mesh**:
// halving the wall length weakens the bending eightfold. A model calibrated at
// one resolution is not calibrated at another.
//
// For a mesh-independent form, hold B = k_bend * h^3 / 4 fixed instead and
// recompute k per vertex from the current spacing; B is then a bending modulus
// with energy 4 pi^4 B A^2 / L^4 and no h in it. (Raised by the pavement-cell
// session, which hit this as mesh-dependent lobe spacing. No model in this
// tree uses BendingChain - the apical-hook model has no bending term - so
// nothing here is affected, and it is left as-is rather than changed under
// anything that may be calibrated against it elsewhere.)
//
// (New in v2.)
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
