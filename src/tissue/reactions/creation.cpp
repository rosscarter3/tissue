//
// Creation reactions: constant and state-dependent production terms for cell
// variables. Ported from legacy creation.cc; identical math, parallel over
// cells (each iteration writes only its own row).
//
#include <cmath>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// dc/dt += k_c
class CreationZero : public Reaction {
public:
  CreationZero(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::Zero", p, i, 1, {1}, {"k_c"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &cellDerivs,
              Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kc = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] += kc;
    });
  }
  void derivsWithAbs(Tissue &T, Matrix &cd, Matrix &wd, Matrix &vd,
                     Matrix &cellDerivs, Matrix &wallDerivs,
                     Matrix &vertexDerivs, Matrix &sdydtCell, Matrix &,
                     Matrix &) override {
    derivs(T, cd, wd, vd, cellDerivs, wallDerivs, vertexDerivs);
    const size_t cIndex = variableIndex(0, 0);
    const double kc = parameter(0);
    for (size_t cell = 0; cell < sdydtCell.rows(); ++cell)
      sdydtCell[cell][cIndex] += kc;
  }
};
TISSUE_REGISTER_REACTION(CreationZero, "Creation::Zero", "CreationZero")

// dc/dt += k_c * X
class CreationOne : public Reaction {
public:
  CreationOne(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::One", p, i, 1, {1, 1}, {"k_c"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double kc = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] += kc * cellData[cell][xIndex];
    });
  }
  void derivsWithAbs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
                     Matrix &cellDerivs, Matrix &, Matrix &, Matrix &sdydtCell,
                     Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double kc = parameter(0);
    for (size_t cell = 0; cell < cellDerivs.rows(); ++cell) {
      double value = kc * cellData[cell][xIndex];
      cellDerivs[cell][cIndex] += value;
      sdydtCell[cell][cIndex] += value;
    }
  }
};
TISSUE_REGISTER_REACTION(CreationOne, "Creation::One", "CreationOne")

// dc/dt += k_c * X * Y
class CreationTwo : public Reaction {
public:
  CreationTwo(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::Two", p, i, 1, {1, 2}, {"k_c"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const size_t yIndex = variableIndex(1, 1);
    const double kc = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] +=
            kc * cellData[cell][xIndex] * cellData[cell][yIndex];
    });
  }
};
TISSUE_REGISTER_REACTION(CreationTwo, "Creation::Two", "CreationTwo")

// dc/dt += k_c * X * Y * Z
class CreationThree : public Reaction {
public:
  CreationThree(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::Three", p, i, 1, {1, 3}, {"k_c"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const size_t yIndex = variableIndex(1, 1);
    const size_t zIndex = variableIndex(1, 2);
    const double kc = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] += kc * cellData[cell][xIndex] *
                                    cellData[cell][yIndex] *
                                    cellData[cell][zIndex];
    });
  }
};
TISSUE_REGISTER_REACTION(CreationThree, "Creation::Three", "CreationThree")

// Hill-function production dependent on cell distance from the origin:
// outside (R_sign>0) or inside (R_sign<0) a sphere of radius R.
class CreationSpatialSphere : public Reaction {
public:
  CreationSpatialSphere(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::SpatialSphere", p, i, 4, {1},
              {"V_max", "R_Hill", "n_Hill", "R_sign"});
    if (p[3] != -1.0 && p[3] != 1.0)
      throw std::runtime_error(
          "Creation::SpatialSphere: R_sign (parameter 3) must be +/-1.");
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double vMax = parameter(0);
    const double nHill = parameter(2);
    const double powK = std::pow(parameter(1), nHill);
    const bool outside = parameter(3) > 0.0;
    const size_t dimension = vertexData.cols();
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        Vec3 center = T.cellPosition(cell, vertexData);
        double r2 = 0.0;
        for (size_t d = 0; d < dimension; ++d)
          r2 += center[d] * center[d];
        double powR = std::pow(std::sqrt(r2), nHill);
        cellDerivs[cell][cIndex] +=
            outside ? vMax * powR / (powK + powR) : vMax * powK / (powK + powR);
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CreationSpatialSphere, "Creation::SpatialSphere",
                         "CreationSpatialSphere")

} // namespace
} // namespace tissue
