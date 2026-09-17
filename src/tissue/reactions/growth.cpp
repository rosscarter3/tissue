//
// Wall growth reactions: irreversible resting-length growth of walls (edges).
// Ported from legacy growth.cc.
//
#include <cmath>
#include <stdexcept>

#include <vector>
#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// dL/dt += k_growth [* L] [* (1 - L/L_trunc)]
class WallGrowthConstant : public Reaction {
public:
  WallGrowthConstant(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error(
          "WallGrowth::Constant: uses two or three parameters "
          "(k_growth, linearFlag, [L_trunc]).");
    configure("WallGrowth::Constant", p, i, p.size(), {1},
              p.size() == 3
                  ? std::vector<std::string>{"k_growth", "linearFlag", "L_trunc"}
                  : std::vector<std::string>{"k_growth", "linearFlag"});
  }
  void derivs(Tissue &, Matrix &, Matrix &wallData, Matrix &,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const double k = parameter(0);
    const bool linear = parameter(1) == 1.0;
    const bool truncated = numParameter() > 2;
    const double lTrunc = truncated ? parameter(2) : 0.0;
    parallelFor(wallDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        double arg = k;
        if (linear)
          arg *= wallData[w][lengthIndex];
        if (truncated)
          arg *= (1.0 - wallData[w][lengthIndex] / lTrunc);
        wallDerivs[w][lengthIndex] += arg;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthConstant, "WallGrowth::Constant",
                         "WallGrowthConstant")

// Stress/strain-gated wall growth:
//   stress = sum of wall stress variables (strain_flag=0) or geometric strain
//   dL/dt += k_growth * (stress - threshold) [* L] [* (1 - L/L_trunc)]
// applied when threshold==0 (always, shrinkage allowed) or stress>threshold.
class WallGrowthStress : public Reaction {
public:
  WallGrowthStress(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 4 && p.size() != 5)
      throw std::runtime_error(
          "WallGrowth::Stress: uses four or five parameters "
          "(k_growth, stress_threshold, strain_flag, linear_flag, [L_trunc]).");
    if (p[2] != 0.0 && p[2] != 1.0)
      throw std::runtime_error(
          "WallGrowth::Stress: strain_flag (parameter 2) must be 0 or 1.");
    if (p[3] != 0.0 && p[3] != 1.0)
      throw std::runtime_error(
          "WallGrowth::Stress: linear_flag (parameter 3) must be 0 or 1.");
    if (!((i.size() == 1 || i.size() == 2) && i[0].size() == 1) ||
        (p[2] == 0.0 && (i.size() < 2 || i[1].empty())))
      throw std::runtime_error(
          "WallGrowth::Stress: wall length index in level 0; stress variable "
          "indices in level 1 (required when strain_flag=0).");
    std::vector<std::string> ids{"k_growth", "stress_threshold", "strain_flag",
                                 "linear_flag"};
    if (p.size() == 5)
      ids.push_back("L_trunc");
    setId("WallGrowth::Stress");
    setParameter(p);
    setVariableIndex(i);
    setParameterIds(std::move(ids));
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const double k = parameter(0);
    const double threshold = parameter(1);
    const bool strainMode = parameter(2) == 1.0;
    const bool linear = parameter(3) == 1.0;
    const bool truncated = numParameter() > 4;
    const double lTrunc = truncated ? parameter(4) : 0.0;
    parallelFor(wallDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        double stress = 0.0;
        if (!strainMode) {
          for (size_t k2 = 0; k2 < numVariableIndex(1); ++k2)
            stress += wallData[w][variableIndex(1, k2)];
        } else {
          double d = T.wallLengthFromVertices(w, vertexData);
          stress = (d - wallData[w][lengthIndex]) / wallData[w][lengthIndex];
        }
        if (threshold == 0.0 || stress > threshold) {
          double growthRate = k * (stress - threshold);
          if (linear)
            growthRate *= wallData[w][lengthIndex];
          if (truncated)
            growthRate *= (1.0 - wallData[w][lengthIndex] / lTrunc);
          wallDerivs[w][lengthIndex] += growthRate;
        }
      }
    });
  }
  void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &wallData,
                     Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                     Matrix &vertexDerivs, Matrix &, Matrix &,
                     Matrix &) override {
    // Legacy duplicates derivs and never touches the noise matrices.
    derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs,
           vertexDerivs);
  }
};
TISSUE_REGISTER_REACTION(WallGrowthStress, "WallGrowth::Stress",
                         "WallGrowthStress")

// Strain-gated growth applied as a direct per-step state update (legacy
// derivs is empty; growth uses explicit Euler in update()). Almansi strain:
// eps = (1 - (L0/d)^2)/2. With 5 parameters, growth only happens when every
// cell's stored velocity is below velocity_threshold (mechanical equilibrium)
// and uses a 2-term Taylor form. Note: legacy reads parameter(4) and
// variableIndex(0,1) even in the 4-parameter form (out-of-range, undefined);
// v2 skips the equilibrium gate in that case.
class WallGrowthStrain : public Reaction {
public:
  WallGrowthStrain(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 4 && p.size() != 5)
      throw std::runtime_error(
          "WallGrowth::Strain: uses four or five parameters "
          "(k_growth, s_threshold, strain_flag, linear_flag, "
          "[velocity_threshold]).");
    if (p[2] != 0.0)
      throw std::runtime_error(
          "WallGrowth::Strain: only strain_flag=0 (Almansi) implemented.");
    if (p[3] != 0.0 && p[3] != 1.0)
      throw std::runtime_error(
          "WallGrowth::Strain: linear_flag (parameter 3) must be 0 or 1.");
    if (i.size() != 1 || (i[0].size() != 1 && i[0].size() != 2))
      throw std::runtime_error(
          "WallGrowth::Strain: level 0 holds wall length index "
          "[+ cell velocity-store index for the 5-parameter form].");
    if (p.size() == 5 && i[0].size() != 2)
      throw std::runtime_error(
          "WallGrowth::Strain: 5-parameter form needs the cell velocity-store "
          "index as second index.");
    std::vector<std::string> ids{"k_growth", "s_threshold", "strain_flag",
                                 "linear_flag"};
    if (p.size() == 5)
      ids.push_back("velocity_threshold");
    setId("WallGrowth::Strain");
    setParameter(p);
    setVariableIndex(i);
    setParameterIds(std::move(ids));
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, double h) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const double k = parameter(0);
    const double threshold = parameter(1);
    bool equil = true;
    if (numParameter() == 5) {
      const double velocityThreshold = parameter(4);
      const size_t velocityStoreIndex = variableIndex(0, 1);
      for (size_t c = 0; c < cellData.rows(); ++c)
        if (cellData[c][velocityStoreIndex] > velocityThreshold) {
          equil = false;
          break;
        }
    }
    if (!equil)
      return;
    for (size_t w = 0; w < T.numWall(); ++w) {
      double restingL = wallData[w][lengthIndex];
      double d = T.wallLengthFromVertices(w, vertexData);
      double strain = 0.5 * (1.0 - (restingL / d) * (restingL / d));
      if (strain > threshold) {
        if (numParameter() != 5) {
          wallData[w][lengthIndex] += h * restingL * k * (strain - threshold);
        } else {
          double factor1 = k * h * (strain - threshold) / (1.0 - 2.0 * strain);
          factor1 = factor1 + factor1 * factor1;
          wallData[w][lengthIndex] += factor1 * restingL;
        }
      }
    }
  }
};
TISSUE_REGISTER_REACTION(WallGrowthStrain, "WallGrowth::Strain",
                         "WallGrowthStrain")

// Constant growth of the internal (center-triangulation) edges stored in
// cell variables: dL_k/dt += k_growth [* L_k] [* (1 - L_k/L_trunc)].
class CTWallGrowthConstant : public Reaction {
public:
  CTWallGrowthConstant(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error(
          "CenterTriangulation::WallGrowth::Constant: uses two or three "
          "parameters (k_growth, linearFlag, [L_trunc]).");
    configure("WallGrowth::CenterTriangulation::Constant", p, i, p.size(), {1},
              p.size() == 3
                  ? std::vector<std::string>{"k_growth", "linearFlag", "L_trunc"}
                  : std::vector<std::string>{"k_growth", "linearFlag"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t lengthStartIndex = variableIndex(0, 0) + 3; // 3D center first
    const double k = parameter(0);
    const bool linear = parameter(1) == 1.0;
    const bool truncated = numParameter() > 2;
    const double lTrunc = truncated ? parameter(2) : 0.0;
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t c = b; c < e; ++c) {
        for (size_t k2 = 0; k2 < T.cell(c).numVertex(); ++k2) {
          double arg = k;
          if (linear)
            arg *= cellData[c][k2 + lengthStartIndex];
          if (truncated)
            arg *= (1.0 - cellData[c][k2 + lengthStartIndex] / lTrunc);
          cellDerivs[c][k2 + lengthStartIndex] += arg;
        }
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CTWallGrowthConstant,
                         "WallGrowth::CenterTriangulation::Constant",
                         "CenterTriangulation::WallGrowth::Constant")

// Wall stress, either summed from listed wall variables or computed as the
// Almansi-style extension (d - L)/L. Shared by the spatial and Hill variants
// below; `useStrain` corresponds to legacy's stress/strain flag.
static inline double wallStress(const Tissue &T, size_t w, Matrix &wallData,
                                Matrix &vertexData, size_t lengthIndex,
                                const Reaction &r, bool useStrain) {
  if (!useStrain) {
    double s = 0.0;
    for (size_t k = 0; k < r.numVariableIndex(1); ++k)
      s += wallData[w][r.variableIndex(1, k)];
    return s;
  }
  const double d = T.wallLengthFromVertices(w, vertexData);
  return (d - wallData[w][lengthIndex]) / wallData[w][lengthIndex];
}

// Growth above a stress threshold, damped by distance from the tissue point
// with the largest value of a chosen vertex coordinate.
class WallGrowthStressSpatial : public Reaction {
public:
  WallGrowthStressSpatial(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2)
      throw std::runtime_error(
          "WallGrowth::StressSpatial: level 0 = wall length index and the "
          "vertex coordinate defining the spatial maximum; level 1 = the wall "
          "variables summed as stress.");
    configure("WallGrowth::StressSpatial", p, i, 6, {2, i[1].size()},
              {"k_growth", "stress_threshold", "K_hill", "n_hill",
               "strain_flag", "linear_flag"});
    kPow_ = std::pow(p[2], p[3]);
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0), sI = variableIndex(0, 1);
    const size_t dim = vertexData.cols();
    size_t maxI = 0;
    for (size_t v = 1; v < vertexData.rows(); ++v)
      if (vertexData[v][sI] > vertexData[maxI][sI])
        maxI = v;
    std::vector<double> maxPos(dim);
    for (size_t d = 0; d < dim; ++d)
      maxPos[d] = vertexData[maxI][d];
    for (size_t w = 0; w < T.numWall(); ++w) {
      const double stress = wallStress(T, w, wallData, vertexData, lengthIndex,
                                       *this, parameter(4) != 0.0);
      if (stress <= parameter(1))
        continue;
      const size_t v1 = T.wall(w).vertex1, v2 = T.wall(w).vertex2;
      double maxDistance = 0.0;
      for (size_t d = 0; d < dim; ++d) {
        const double pos = 0.5 * (vertexData[v1][d] + vertexData[v2][d]);
        maxDistance += (maxPos[d] - pos) * (maxPos[d] - pos);
      }
      maxDistance = std::sqrt(maxDistance);
      const double spatial =
          kPow_ / (kPow_ + std::pow(maxDistance, parameter(3)));
      double rate = parameter(0) * (stress - parameter(1)) * spatial;
      if (parameter(5))
        rate *= wallData[w][lengthIndex];
      wallDerivs[w][lengthIndex] += rate;
    }
  }

private:
  double kPow_ = 1.0;
};
TISSUE_REGISTER_REACTION(WallGrowthStressSpatial, "WallGrowth::StressSpatial",
                         "WallGrowthStressSpatial")

// As above, but the spatial distance is measured along the chosen coordinate
// alone rather than in full space.
class WallGrowthStressSpatialSingle : public Reaction {
public:
  WallGrowthStressSpatialSingle(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2)
      throw std::runtime_error(
          "WallGrowth::StressSpatialSingle: level 0 = wall length index and "
          "the vertex coordinate; level 1 = wall variables summed as stress.");
    configure("WallGrowth::StressSpatialSingle", p, i, 6, {2, i[1].size()},
              {"k_growth", "stress_threshold", "K_hill", "n_hill",
               "strain_flag", "linear_flag"});
    kPow_ = std::pow(p[2], p[3]);
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0), sI = variableIndex(0, 1);
    double sMax = vertexData[0][sI];
    for (size_t v = 1; v < vertexData.rows(); ++v)
      sMax = std::max(sMax, vertexData[v][sI]);
    for (size_t w = 0; w < T.numWall(); ++w) {
      const double stress = wallStress(T, w, wallData, vertexData, lengthIndex,
                                       *this, parameter(4) != 0.0);
      if (stress <= parameter(1))
        continue;
      const size_t v1 = T.wall(w).vertex1, v2 = T.wall(w).vertex2;
      const double maxDistance =
          sMax - 0.5 * (vertexData[v1][sI] + vertexData[v2][sI]);
      const double spatial =
          kPow_ / (kPow_ + std::pow(maxDistance, parameter(3)));
      double rate = parameter(0) * (stress - parameter(1)) * spatial;
      if (parameter(5))
        rate *= wallData[w][lengthIndex];
      wallDerivs[w][lengthIndex] += rate;
    }
  }

private:
  double kPow_ = 1.0;
};
TISSUE_REGISTER_REACTION(WallGrowthStressSpatialSingle,
                         "WallGrowth::StressSpatialSingle",
                         "WallGrowthStressSpatialSingle")

// Growth above a stress threshold, with the rate raised by a Hill function of
// a concentration in the adjacent cells. Both neighbours contribute, so the
// Hill factor runs to 2 for an interior wall and 1 at the boundary.
class WallGrowthStressConcentrationHill : public Reaction {
public:
  WallGrowthStressConcentrationHill(const ParameterList &p,
                                    const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2)
      throw std::runtime_error(
          "WallGrowth::StressConcentrationHill: level 0 = wall length index "
          "and the cell concentration index; level 1 = wall variables summed "
          "as stress.");
    configure("WallGrowth::StressConcentrationHill", p, i, 7,
              {2, i[1].size()},
              {"k_growthConst", "k_growthHill", "K_hill", "n_hill",
               "stress_threshold", "strain_flag", "linear_flag"});
    kPow_ = std::pow(p[2], p[3]);
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &, Matrix &wallDerivs,
              Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0),
                 concIndex = variableIndex(0, 1);
    for (size_t w = 0; w < T.numWall(); ++w) {
      const double stress = wallStress(T, w, wallData, vertexData, lengthIndex,
                                       *this, parameter(5) != 0.0);
      if (stress <= parameter(4))
        continue;
      double hill = 0.0;
      for (size_t c : {T.wall(w).cell1, T.wall(w).cell2})
        if (!Tissue::isBackground(c)) {
          const double cp = std::pow(cellData[c][concIndex], parameter(3));
          hill += cp / (kPow_ + cp);
        }
      double rate =
          (parameter(0) + hill * parameter(1)) * (stress - parameter(4));
      if (parameter(6))
        rate *= wallData[w][lengthIndex];
      wallDerivs[w][lengthIndex] += rate;
    }
  }

private:
  double kPow_ = 1.0;
};
TISSUE_REGISTER_REACTION(WallGrowthStressConcentrationHill,
                         "WallGrowth::StressConcentrationHill",
                         "WallGrowthStressConcentrationHill")

// Constant stress-driven growth with a separate rate factor for walls facing
// the background, i.e. the epidermis.
class WallGrowthConstantStressEpidermalAsymmetric : public Reaction {
public:
  WallGrowthConstantStressEpidermalAsymmetric(const ParameterList &p,
                                              const IndexLevels &i) {
    configure("WallGrowth::ConstantStressEpidermalAsymmetric", p, i, 2, {1},
              {"k_growth", "f_epidermis"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        double k = parameter(0);
        if (Tissue::isBackground(T.wall(w).cell1) ||
            Tissue::isBackground(T.wall(w).cell2))
          k *= parameter(1);
        const double d = T.wallLengthFromVertices(w, vertexData);
        if (d > wallData[w][lengthIndex])
          wallDerivs[w][lengthIndex] += k * (d - wallData[w][lengthIndex]);
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthConstantStressEpidermalAsymmetric,
                         "WallGrowth::ConstantStressEpidermalAsymmetric",
                         "WallGrowthConstantStressEpidermalAsymmetric")

// dL/dt += k_L (F - phi) for the summed wall force F above a threshold.
class WallGrowthForce : public Reaction {
public:
  WallGrowthForce(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1)
      throw std::runtime_error(
          "WallGrowth::Force: level 0 = wall length index; level 1 = the wall "
          "variables summed as force.");
    configure("WallGrowth::Force", p, i, 2, {1, i[1].size()}, {"k_L", "phi"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        double F = 0.0;
        for (size_t j = 0; j < numVariableIndex(1); ++j)
          F += wallData[w][variableIndex(1, j)];
        const double arg = F - parameter(1);
        if (arg > 0.0)
          wallDerivs[w][lengthIndex] += parameter(0) * arg;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthForce, "WallGrowth::Force")

} // namespace

} // namespace tissue
