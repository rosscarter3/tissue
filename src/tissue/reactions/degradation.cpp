//
// Degradation reactions: first order and state-dependent decay of cell
// variables. Ported from legacy degradation.cc.
//
#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// dc/dt -= k_d * c
class DegradationOne : public Reaction {
public:
  DegradationOne(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::One", p, i, 1, {1}, {"k_d"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kd = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] -= kd * cellData[cell][cIndex];
    });
  }
  void derivsWithAbs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
                     Matrix &cellDerivs, Matrix &, Matrix &, Matrix &sdydtCell,
                     Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kd = parameter(0);
    for (size_t cell = 0; cell < cellDerivs.rows(); ++cell) {
      double value = kd * cellData[cell][cIndex];
      cellDerivs[cell][cIndex] -= value;
      sdydtCell[cell][cIndex] += value;
    }
  }
};
TISSUE_REGISTER_REACTION(DegradationOne, "Degradation::One", "DegradationOne")

// dc/dt -= k_d * X * c
class DegradationTwo : public Reaction {
public:
  DegradationTwo(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::Two", p, i, 1, {1, 1}, {"k_d"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double kd = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] -=
            kd * cellData[cell][xIndex] * cellData[cell][cIndex];
    });
  }
};
TISSUE_REGISTER_REACTION(DegradationTwo, "Degradation::Two", "DegradationTwo")

// dc/dt -= k_d * (prod_i X_i) * c
class DegradationN : public Reaction {
public:
  DegradationN(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::N", p, i, 1, {1, kAnyCount}, {"k_d"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kd = parameter(0);
    const size_t numX = numVariableIndex(1);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        double product = kd;
        for (size_t k = 0; k < numX; ++k)
          product *= cellData[cell][variableIndex(1, k)];
        cellDerivs[cell][cIndex] -= product * cellData[cell][cIndex];
      }
    });
  }
};
TISSUE_REGISTER_REACTION(DegradationN, "Degradation::N", "DegradationN")

} // namespace
} // namespace tissue
