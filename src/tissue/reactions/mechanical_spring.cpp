//
// Spring mechanics. Ported from legacy mechanicalSpring.cc; the shared
// coefficient idiom is coeff = K*(1/L_rest - 1/d), zeroed only when both
// d<=0 and L_rest<=0, scaled by frac_adh when stretched (d > L_rest), applied
// as (x_a - x_b)*coeff with -= on a and += on b.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Asymmetric wall spring moving the wall's two vertices. Optional: save the
// scalar force K_eff*(d-L)/L into a wall variable (level 1), and use K_force2
// for walls whose type variable (level 2) equals 1.
class WallMechanicsSpring : public Reaction {
public:
  WallMechanicsSpring(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error(
          "WallMechanics::Spring: uses two or three parameters "
          "(K_force, frac_adhesion, [K_force2]).");
    bool okLevels =
        i.size() == 1 || i.size() == 2 ||
        (i.size() == 3 && i[2].size() == 1 &&
         (i[1].empty() || i[1].size() == 1));
    if (!okLevels || i.empty() || (i[0].size() != 1 && i[0].size() != 2))
      throw std::runtime_error(
          "WallMechanics::Spring: wall length index in level 0; optional "
          "force-save wall index in level 1; optional wall-type index in "
          "level 2 (with K_force2).");
    std::vector<std::string> ids{"K_force", "frac_adh"};
    if (p.size() == 3)
      ids.push_back("K_force2");
    setId("WallMechanics::Spring");
    setParameter(p);
    setVariableIndex(i);
    setParameterIds(std::move(ids));
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t wallLengthIndex = variableIndex(0, 0);
    const double kForce = parameter(0);
    const double fracAdh = parameter(1);
    const bool typed = numParameter() == 3 && numVariableIndexLevel() == 3;
    const bool saveForce =
        (numVariableIndexLevel() == 2 && numVariableIndex(1) > 0) ||
        (numParameter() == 3 && numVariableIndexLevel() > 1 &&
         numVariableIndex(1) > 0);
    const size_t dimension = vertexData.cols();
    const size_t numWalls = T.numWall();

    auto wallCoeff = [&](size_t w) {
      const size_t v1 = T.wall(w).vertex1;
      const size_t v2 = T.wall(w).vertex2;
      double d = 0.0;
      for (size_t dim = 0; dim < dimension; ++dim) {
        double diff = vertexData[v1][dim] - vertexData[v2][dim];
        d += diff * diff;
      }
      d = std::sqrt(d);
      const double wallLength = wallData[w][wallLengthIndex];
      double coeff = kForce * (1.0 / wallLength - 1.0 / d);
      if (typed && wallData[w][variableIndex(2, 0)] == 1.0)
        coeff = parameter(2) * (1.0 / wallLength - 1.0 / d);
      if (d <= 0.0 && wallLength <= 0.0)
        coeff = 0.0;
      if (d > wallLength)
        coeff *= fracAdh;
      if (saveForce)
        wallData[w][variableIndex(1, 0)] = coeff * d; // overwrite (legacy)
      return coeff;
    };

    if (numWalls < 2048 || ThreadPool::instance().numThreads() == 1) {
      // Small tissues: serial loop, force accumulation in exact legacy order.
      for (size_t w = 0; w < numWalls; ++w) {
        double coeff = wallCoeff(w);
        const size_t v1 = T.wall(w).vertex1;
        const size_t v2 = T.wall(w).vertex2;
        for (size_t dim = 0; dim < dimension; ++dim) {
          double div = (vertexData[v1][dim] - vertexData[v2][dim]) * coeff;
          vertexDerivs[v1][dim] -= div;
          vertexDerivs[v2][dim] += div;
        }
      }
      return;
    }
    // Large tissues: two-pass gather. Pass 1 computes per-wall coefficients;
    // pass 2 lets each vertex sum its own walls' forces — race free and
    // deterministic for any thread count.
    coeff_.resize(numWalls);
    parallelFor(
        numWalls,
        [&](size_t b, size_t e) {
          for (size_t w = b; w < e; ++w)
            coeff_[w] = wallCoeff(w);
        },
        1024);
    parallelFor(
        T.numVertex(),
        [&](size_t b, size_t e) {
          for (size_t v = b; v < e; ++v) {
            for (size_t w : T.vertex(v).walls) {
              const Wall &wall = T.wall(w);
              double sign = wall.vertex1 == v ? -1.0 : 1.0;
              for (size_t dim = 0; dim < dimension; ++dim)
                vertexDerivs[v][dim] +=
                    sign * (vertexData[wall.vertex1][dim] -
                            vertexData[wall.vertex2][dim]) *
                    coeff_[w];
            }
          }
        },
        1024);
  }
  void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &wallData,
                     Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                     Matrix &vertexDerivs, Matrix &, Matrix &,
                     Matrix &) override {
    derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs,
           vertexDerivs); // legacy: no noise contribution
  }

private:
  std::vector<double> coeff_; // per-wall coefficients for the gather path
};
TISSUE_REGISTER_REACTION(WallMechanicsSpring, "WallMechanics::Spring",
                         "VertexFromWallSpring")

// Spring between each cell's center-triangulation center (stored in cell
// variables) and its vertices.
class CTEdgeSpring : public Reaction {
public:
  CTEdgeSpring(const ParameterList &p, const IndexLevels &i) {
    configure("CenterTriangulation::EdgeSpring", p, i, 2, {1},
              {"K_force", "frac_adh"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &vertexDerivs) override {
    const size_t posIndex = variableIndex(0, 0);
    const size_t dimension = vertexData.cols(); // 2D or 3D
    const size_t lengthIndex = posIndex + dimension;
    const double kForce = parameter(0);
    const double fracAdh = parameter(1);
    parallelScatter2(
        T.numCell(), vertexDerivs, cellDerivs,
        [&](size_t b, size_t e, Matrix &vOut, Matrix &cOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            for (size_t k = 0; k < cell.numVertex(); ++k) {
              size_t v = cell.vertices[k];
              double d = 0.0;
              for (size_t dim = 0; dim < dimension; ++dim) {
                double diff = vertexData[v][dim] - cellData[c][posIndex + dim];
                d += diff * diff;
              }
              d = std::sqrt(d);
              double edgeLength = cellData[c][lengthIndex + k];
              double coeff = kForce * (1.0 / edgeLength - 1.0 / d);
              if (d <= 0.0 && edgeLength <= 0.0)
                coeff = 0.0;
              if (d > edgeLength)
                coeff *= fracAdh;
              for (size_t dim = 0; dim < dimension; ++dim) {
                double div =
                    (vertexData[v][dim] - cellData[c][posIndex + dim]) * coeff;
                vOut[v][dim] -= div;
                cOut[c][posIndex + dim] += div;
              }
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(CTEdgeSpring, "CenterTriangulation::EdgeSpring",
                         "CenterTriangulation::Spring")

} // namespace
} // namespace tissue
