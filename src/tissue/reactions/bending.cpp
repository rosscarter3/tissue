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
#include <cmath>
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


// ---------------------------------------------------------------------------
// The legacy Bending:: family (legacy/bending.cc). These walk each cell's
// sorted vertex cycle and act on the turn angle at every vertex, so a wall
// shared by two cells carries one angle spring per cell.
//
// Two index bugs in the legacy sources are fixed here (see "Deliberate fixes
// over legacy" in README.md); both made the reaction apply force to the wrong
// vertex, so neither can be reproduced bit-for-bit and neither should be.
//
// Legacy writes pi as the 6-digit literal 3.14159. That is kept, because the
// rest angle written by Bending::AngleInitiate is compared against the angle
// computed by Bending::Angle and Bending::AngleRelax: they must use the same
// constant, and init files carrying rest angles saved from legacy runs are
// expressed in it. It biases a stored rest angle by 2.7e-6 rad.
// ---------------------------------------------------------------------------


constexpr double kLegacyPi = 3.14159;

// Turn geometry at cell-local vertex k: the two edge vectors meeting there,
// their lengths, and their dot product.
struct Turn {
  size_t jm, j, jp; // global vertex indices, in cycle order
  size_t em;        // wall carrying the rest angle (jm--j)
  double Lm, Lp, dot;
};

inline Turn turnAt(const Tissue &T, size_t cell, size_t k,
                   const Matrix &vertexData) {
  const CellTopo &c = T.cell(cell);
  const size_t n = c.numWall();
  const size_t kp = k + 1 < n ? k + 1 : 0;
  const size_t km = k > 0 ? k - 1 : n - 1;
  Turn t{c.vertices[km], c.vertices[k], c.vertices[kp], c.walls[km], 0, 0, 0};
  const size_t dim = vertexData.cols();
  for (size_t d = 0; d < dim; ++d) {
    const double dm = vertexData[t.j][d] - vertexData[t.jm][d];
    const double dp = vertexData[t.jp][d] - vertexData[t.j][d];
    t.dot += dm * dp;
    t.Lp += dp * dp;
    t.Lm += dm * dm;
  }
  t.Lp = std::sqrt(t.Lp);
  t.Lm = std::sqrt(t.Lm);
  return t;
}

// theta = acos(cos) - pi, in [-pi, 0]; -pi is a straight chain.
inline double restAngle(double cosTheta) {
  if (cosTheta > 1.0)
    cosTheta = 1.0; // legacy clamps the upper end only
  return std::acos(cosTheta) - kLegacyPi;
}

// Pulls each vertex toward the point on the line between its two cycle
// neighbours that matches its own rest-length spacing:
//   F_j = -k (x_j - (x_jp L_em + x_jm L_ep)/(L_em + L_ep))
// Unlike WallMechanics::BendingChain this puts no reaction on the neighbours,
// so it is not momentum conserving - it acts as an external restoring field,
// which is what legacy does and what models using it are calibrated against.
//
// FIX: legacy accumulates into `vertexDerivs[k]`, where k is the cell-local
// wall counter, not the global vertex index. Every cell therefore pushed on
// vertices 0..numWall-1 of the whole tissue regardless of which vertices the
// angle belonged to. Written here to the intended vertex j.
class BendingNeighborCenter : public Reaction {
public:
  BendingNeighborCenter(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Bending::NeighborCenter: one index level with "
                               "one index (the wall length).");
    configure("Bending::NeighborCenter", p, i, 1, {1}, {"k_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t Li = variableIndex(0, 0);
    const double kBend = parameter(0);
    const size_t dim = vertexData.cols();
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numWall();
            for (size_t k = 0; k < n; ++k) {
              const size_t kp = k + 1 < n ? k + 1 : 0;
              const size_t km = k > 0 ? k - 1 : n - 1;
              const size_t j = cell.vertices[k];
              const size_t jp = cell.vertices[kp];
              const size_t jm = cell.vertices[km];
              const size_t ep = cell.walls[k];
              const size_t em = cell.walls[km];
              const double sum = wallData[em][Li] + wallData[ep][Li];
              for (size_t d = 0; d < dim; ++d)
                out[j][d] -= kBend * (vertexData[j][d] -
                                      (vertexData[jp][d] * wallData[em][Li] +
                                       vertexData[jm][d] * wallData[ep][Li]) /
                                          sum);
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(BendingNeighborCenter, "Bending::NeighborCenter")

// Angular spring on the turn at each vertex: E = (k/2)(theta - theta_0)^2 with
// theta_0 stored per wall (set by Bending::AngleInitiate). The three vertex
// forces are exactly -dE/dx and sum to zero.
//
// FIX: legacy accumulates the middle (central-vertex) term into
// `vertexDerivs[jm]` a second time instead of into `vertexDerivs[j]`, so the
// centre vertex felt no force, its predecessor felt two, and the triple did
// not sum to zero. Written here to j.
//
// WART kept from legacy: cos(theta) is clamped to +/-0.999 to keep the
// 1/sqrt(1-F^2) factor finite. Near a straight chain the true limit is finite
// anyway (theta - theta_0 vanishes at the same rate), so the clamp leaves a
// residual force: with theta_0 from AngleInitiate on a straight chain,
// theta - theta_0 floors at acos(0.999) = 0.045 rad rather than 0. Models
// using this reaction are calibrated with that floor present, so it is left
// in; a straight chain is not a stress-free state for it.
class BendingAngle : public Reaction {
public:
  BendingAngle(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Bending::Angle: one index level with one "
                               "index (the wall variable holding the rest "
                               "angle).");
    configure("Bending::Angle", p, i, 1, {1}, {"k_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t Ti = variableIndex(0, 0);
    const double kBend = parameter(0);
    const size_t dim = vertexData.cols();
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            for (size_t k = 0; k < T.cell(c).numWall(); ++k) {
              const Turn t = turnAt(T, c, k, vertexData);
              const double gDenom = 1.0 / (t.Lp * t.Lm);
              double F = t.dot * gDenom;
              if (F > 0.999)
                F = 0.999;
              else if (F < -0.999)
                F = -0.999;
              const double theta = std::acos(F) - kLegacyPi;
              const double f0 = kBend * (theta - wallData[t.em][Ti]) * gDenom /
                                std::sqrt(1.0 - F * F);
              const double f1 = t.dot * gDenom * t.Lm / t.Lp;
              const double f2 = t.dot * gDenom * t.Lp / t.Lm;
              for (size_t d = 0; d < dim; ++d) {
                const double dm = vertexData[t.j][d] - vertexData[t.jm][d];
                const double dp = vertexData[t.jp][d] - vertexData[t.j][d];
                out[t.jm][d] += f0 * (dm * f2 - dp);
                out[t.j][d] += f0 * (dp * (1.0 + f1) - dm * (1.0 + f2));
                out[t.jp][d] += f0 * (dm - dp * f1);
              }
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(BendingAngle, "Bending::Angle")

// Stores the current turn angle at each vertex as the rest angle of the wall
// behind it, making the starting shape stress free for Bending::Angle (up to
// the 0.999 clamp noted above).
class BendingAngleInitiate : public Reaction {
public:
  BendingAngleInitiate(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Bending::AngleInitiate: one index level with "
                               "one index (the rest-angle wall variable).");
    configure("Bending::AngleInitiate", p, i, 0, {1}, {});
  }
  void initiate(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const size_t Ti = variableIndex(0, 0);
    for (size_t c = 0; c < T.numCell(); ++c)
      for (size_t k = 0; k < T.cell(c).numWall(); ++k) {
        const Turn t = turnAt(T, c, k, vertexData);
        wallData[t.em][Ti] = restAngle(t.dot / (t.Lp * t.Lm));
      }
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(BendingAngleInitiate, "Bending::AngleInitiate")

// Relaxes the stored rest angle toward the current angle (plastic bending):
// dtheta_0/dt = -k (theta_0 - theta).
class BendingAngleRelax : public Reaction {
public:
  BendingAngleRelax(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Bending::AngleRelax: one index level with one "
                               "index (the rest-angle wall variable).");
    configure("Bending::AngleRelax", p, i, 1, {1}, {"k_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t Ti = variableIndex(0, 0);
    const double kBend = parameter(0);
    parallelScatter1(
        T.numCell(), wallDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c)
            for (size_t k = 0; k < T.cell(c).numWall(); ++k) {
              const Turn t = turnAt(T, c, k, vertexData);
              const double theta = restAngle(t.dot / (t.Lp * t.Lm));
              out[t.em][Ti] -= kBend * (wallData[t.em][Ti] - theta);
            }
        });
  }
};
TISSUE_REGISTER_REACTION(BendingAngleRelax, "Bending::AngleRelax")

} // namespace
} // namespace tissue
