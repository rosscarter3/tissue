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

// The shared wall-spring loop. Every reaction in legacy's mechanicalSpring.cc
// is this loop with two things changed: which walls act (`active`) and how
// stiff each one is (`stiffness`). The coefficient idiom itself -
// K (1/L - 1/d), zeroed only when both d and L are non-positive, scaled by
// frac_adh once stretched, applied as (x1 - x2) coeff with -= on v1 and += on
// v2 - is identical throughout.
//
// Two execution paths, as before: a serial loop that accumulates in legacy's
// exact order for small tissues, and a two-pass gather for large ones where
// each vertex sums its own walls, which is race free and deterministic for any
// thread count.
template <class Active, class Stiffness>
void wallSpringLoop(Tissue &T, Matrix &wallData, Matrix &vertexData,
                    Matrix &vertexDerivs, std::vector<double> &scratch,
                    size_t wallLengthIndex, double fracAdh, bool saveForce,
                    size_t saveIndex, bool zeroInactiveForce, Active active,
                    Stiffness stiffness) {
  const size_t dimension = vertexData.cols();
  const size_t numWalls = T.numWall();

  auto wallCoeff = [&](size_t w) {
    if (!active(w)) {
      if (saveForce && zeroInactiveForce)
        wallData[w][saveIndex] = 0.0;
      return 0.0;
    }
    const size_t v1 = T.wall(w).vertex1;
    const size_t v2 = T.wall(w).vertex2;
    double d = 0.0;
    for (size_t dim = 0; dim < dimension; ++dim) {
      const double diff = vertexData[v1][dim] - vertexData[v2][dim];
      d += diff * diff;
    }
    d = std::sqrt(d);
    const double wallLength = wallData[w][wallLengthIndex];
    double coeff = stiffness(w) * (1.0 / wallLength - 1.0 / d);
    if (d <= 0.0 && wallLength <= 0.0)
      coeff = 0.0;
    if (d > wallLength)
      coeff *= fracAdh;
    if (saveForce)
      wallData[w][saveIndex] = coeff * d; // overwrite, as legacy does
    return coeff;
  };

  if (numWalls < 2048 || ThreadPool::instance().numThreads() == 1) {
    for (size_t w = 0; w < numWalls; ++w) {
      const double coeff = wallCoeff(w);
      const size_t v1 = T.wall(w).vertex1;
      const size_t v2 = T.wall(w).vertex2;
      for (size_t dim = 0; dim < dimension; ++dim) {
        const double div = (vertexData[v1][dim] - vertexData[v2][dim]) * coeff;
        vertexDerivs[v1][dim] -= div;
        vertexDerivs[v2][dim] += div;
      }
    }
    return;
  }
  scratch.resize(numWalls);
  parallelFor(
      numWalls,
      [&](size_t b, size_t e) {
        for (size_t w = b; w < e; ++w)
          scratch[w] = wallCoeff(w);
      },
      1024);
  parallelFor(
      T.numVertex(),
      [&](size_t b, size_t e) {
        for (size_t v = b; v < e; ++v)
          for (size_t w : T.vertex(v).walls) {
            const Wall &wall = T.wall(w);
            const double sign = wall.vertex1 == v ? -1.0 : 1.0;
            for (size_t dim = 0; dim < dimension; ++dim)
              vertexDerivs[v][dim] += sign *
                                      (vertexData[wall.vertex1][dim] -
                                       vertexData[wall.vertex2][dim]) *
                                      scratch[w];
          }
      },
      1024);
}

// True when the wall lies on the tissue boundary.
inline bool onBoundaryWall(const Tissue &T, size_t w) {
  return Tissue::isBackground(T.wall(w).cell1) ||
         Tissue::isBackground(T.wall(w).cell2);
}

// True when a cell touches the background (so its walls are "epidermal").
inline bool cellTouchesBackground(const Tissue &T, size_t c) {
  if (Tissue::isBackground(c))
    return false;
  for (size_t w : T.cell(c).walls)
    if (Tissue::isBackground(T.wall(w).otherCell(c)))
      return true;
  return false;
}

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
    const double kForce = parameter(0);
    const bool typed = numParameter() == 3 && numVariableIndexLevel() == 3;
    const size_t typeIndex = typed ? variableIndex(2, 0) : 0;
    const bool saveForce =
        (numVariableIndexLevel() == 2 && numVariableIndex(1) > 0) ||
        (numParameter() == 3 && numVariableIndexLevel() > 1 &&
         numVariableIndex(1) > 0);
    wallSpringLoop(
        T, wallData, vertexData, vertexDerivs, coeff_, variableIndex(0, 0),
        parameter(1), saveForce, saveForce ? variableIndex(1, 0) : 0,
        /*zeroInactiveForce=*/false, [](size_t) { return true; },
        [&](size_t w) {
          return (typed && wallData[w][typeIndex] == 1.0) ? parameter(2)
                                                          : kForce;
        });
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

// Spring acting only on walls that lie on the tissue boundary - a stiffer (or
// softer) epidermal layer without giving those walls their own variables.
class WallMechanicsSpringEpidermal : public Reaction {
public:
  WallMechanicsSpringEpidermal(const ParameterList &p, const IndexLevels &i) {
    if (i.empty() || i.size() > 2 || i[0].size() != 1 ||
        (i.size() == 2 && i[1].size() != 1))
      throw std::runtime_error(
          "WallMechanics::SpringEpidermal: wall length index in level 0; "
          "optional wall index to save the force in level 1.");
    configure("WallMechanics::SpringEpidermal", p, i, 2,
              std::vector<size_t>(i.size(), 1), {"K_force", "frac_adh"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const bool save = numVariableIndexLevel() > 1;
    wallSpringLoop(
        T, wallData, vertexData, vertexDerivs, coeff_, variableIndex(0, 0),
        parameter(1), save, save ? variableIndex(1, 0) : 0,
        /*zeroInactiveForce=*/true,
        [&](size_t w) { return onBoundaryWall(T, w); },
        [&](size_t) { return parameter(0); });
  }

private:
  std::vector<double> coeff_;
};
TISSUE_REGISTER_REACTION(WallMechanicsSpringEpidermal,
                         "WallMechanics::SpringEpidermal",
                         "VertexFromEpidermalWallSpring")

// As SpringEpidermal, but the whole *cell* counts as epidermal: every wall of
// a cell that touches the background acts, not just the boundary wall itself.
class WallMechanicsSpringEpidermalCell : public Reaction {
public:
  WallMechanicsSpringEpidermalCell(const ParameterList &p,
                                   const IndexLevels &i) {
    if (i.empty() || i.size() > 2 || i[0].size() != 1 ||
        (i.size() == 2 && i[1].size() != 1))
      throw std::runtime_error(
          "WallMechanics::SpringEpidermalCell: wall length index in level 0; "
          "optional wall index to save the force in level 1.");
    configure("WallMechanics::SpringEpidermalCell", p, i, 2,
              std::vector<size_t>(i.size(), 1), {"K_force", "frac_adh"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const bool save = numVariableIndexLevel() > 1;
    wallSpringLoop(
        T, wallData, vertexData, vertexDerivs, coeff_, variableIndex(0, 0),
        parameter(1), save, save ? variableIndex(1, 0) : 0,
        /*zeroInactiveForce=*/true,
        [&](size_t w) {
          const Wall &wall = T.wall(w);
          return Tissue::isBackground(wall.cell1) ||
                 cellTouchesBackground(T, wall.cell1) ||
                 Tissue::isBackground(wall.cell2) ||
                 cellTouchesBackground(T, wall.cell2);
        },
        [&](size_t) { return parameter(0); });
  }

private:
  std::vector<double> coeff_;
};
TISSUE_REGISTER_REACTION(WallMechanicsSpringEpidermalCell,
                         "WallMechanics::SpringEpidermalCell",
                         "VertexFromEpidermalCellWallSpring")

// Stiffness set by an *inhibitory* Hill function of a cell concentration,
// summed over the wall's two cells: K = K_min + K_max (f(c1) + f(c2)) with
// f(c) = K_Hill^n/(K_Hill^n + c^n). A wall on the boundary gets only one
// contribution, so its stiffness runs to K_min + K_max rather than
// K_min + 2 K_max - the same boundary asymmetry legacy has, and preserved.
class WallMechanicsSpringConcentrationHill : public Reaction {
public:
  WallMechanicsSpringConcentrationHill(const ParameterList &p,
                                       const IndexLevels &i) {
    if (i.empty() || i.size() > 2 || i[0].size() != 2 ||
        (i.size() == 2 && i[1].size() != 1))
      throw std::runtime_error(
          "WallMechanics::SpringConcentrationHill: wall length and cell "
          "concentration indices in level 0; optional wall index to save the "
          "force in level 1.");
    configure("WallMechanics::SpringConcentrationHill", p, i, 5,
              i.size() == 2 ? std::vector<size_t>{2, 1}
                            : std::vector<size_t>{2},
              {"K_min", "K_max", "K_Hill", "n_Hill", "frac_adh"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &, Matrix &,
              Matrix &vertexDerivs) override {
    const size_t concIndex = variableIndex(0, 1);
    const double kPow = std::pow(parameter(2), parameter(3));
    const bool save = numVariableIndexLevel() > 1;
    auto fac = [&](size_t c) {
      if (Tissue::isBackground(c))
        return 0.0;
      return kPow / (kPow + std::pow(cellData[c][concIndex], parameter(3)));
    };
    wallSpringLoop(
        T, wallData, vertexData, vertexDerivs, coeff_, variableIndex(0, 0),
        parameter(4), save, save ? variableIndex(1, 0) : 0,
        /*zeroInactiveForce=*/false, [](size_t) { return true; },
        [&](size_t w) {
          return parameter(0) +
                 parameter(1) * (fac(T.wall(w).cell1) + fac(T.wall(w).cell2));
        });
  }

private:
  std::vector<double> coeff_;
};
TISSUE_REGISTER_REACTION(WallMechanicsSpringConcentrationHill,
                         "WallMechanics::SpringConcentrationHill",
                         "VertexFromWallSpringConcentrationHill")

// WallMechanics::Spring restricted to boundary walls, with the same optional
// second stiffness selected by a wall type variable. Unlike SpringEpidermal
// this leaves the saved force on interior walls untouched rather than zeroing
// it, which is legacy behaviour.
class VertexFromWallBoundarySpring : public Reaction {
public:
  VertexFromWallBoundarySpring(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error("VertexFromWallBoundarySpring: uses two or "
                               "three parameters (K_force, frac_adhesion, "
                               "[K_force2]).");
    const bool ok =
        !i.empty() && i[0].size() == 1 &&
        (i.size() == 1 || i.size() == 2 ||
         (i.size() == 3 && i[2].size() == 1 &&
          (i[1].empty() || i[1].size() == 1)));
    if (!ok)
      throw std::runtime_error(
          "VertexFromWallBoundarySpring: wall length index in level 0; "
          "optional force-save wall index in level 1; optional wall-type "
          "index in level 2 (with K_force2).");
    std::vector<std::string> ids{"K_force", "frac_adh"};
    if (p.size() == 3)
      ids.push_back("K_force2");
    setId("VertexFromWallBoundarySpring");
    setParameter(p);
    setVariableIndex(i);
    setParameterIds(std::move(ids));
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const bool typed = numParameter() == 3 && numVariableIndexLevel() == 3;
    const size_t typeIndex = typed ? variableIndex(2, 0) : 0;
    const bool save = numVariableIndexLevel() == 2 ||
                      (numParameter() == 3 && numVariableIndexLevel() > 1 &&
                       numVariableIndex(1) == 1);
    wallSpringLoop(
        T, wallData, vertexData, vertexDerivs, coeff_, variableIndex(0, 0),
        parameter(1), save, save ? variableIndex(1, 0) : 0,
        /*zeroInactiveForce=*/false,
        [&](size_t w) { return onBoundaryWall(T, w); },
        [&](size_t w) {
          return (typed && wallData[w][typeIndex] == 1.0) ? parameter(2)
                                                          : parameter(0);
        });
  }

private:
  std::vector<double> coeff_;
};
TISSUE_REGISTER_REACTION(VertexFromWallBoundarySpring,
                         "VertexFromWallBoundarySpring")

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
