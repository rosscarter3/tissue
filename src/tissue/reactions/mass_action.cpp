//
// Mass action reactions. Ported from legacy massAction.cc.
//
#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// r1 + r2 -> P at rate k1*[r1]*[r2]. Legacy quirk preserved: no positivity
// guard, and r1==r2 consumes twice.
class MassActionTwoToOne : public Reaction {
public:
  MassActionTwoToOne(const ParameterList &p, const IndexLevels &i) {
    // Level 0: {r1, r2, P}; level 1 must be present and empty ("1 2 3 0").
    configure("TwoToOne", p, i, 1, {3, 0}, {"k1"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t r1I = variableIndex(0, 0);
    const size_t r2I = variableIndex(0, 1);
    const size_t pI = variableIndex(0, 2);
    const double k1 = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        double fac = k1 * cellData[cell][r1I] * cellData[cell][r2I];
        cellDerivs[cell][r1I] -= fac;
        cellDerivs[cell][r2I] -= fac;
        cellDerivs[cell][pI] += fac;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionTwoToOne, "MassAction::TwoToOne",
                         "MassActionTwoToOne")

// General reactants -> products; rate k_f * prod(reactants); only applied
// when rate > 0 (legacy guard).
class MassActionGeneral : public Reaction {
public:
  MassActionGeneral(const ParameterList &p, const IndexLevels &i) {
    configure("MassAction::General", p, i, 1, {kAnyCount, kAnyCount}, {"k_f"});
    if (numVariableIndex(0) < 1)
      throw std::runtime_error(
          "MassAction::General: at least one reactant index required.");
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const double kf = parameter(0);
    const size_t numReactant = numVariableIndex(0);
    const size_t numProduct = numVariableIndex(1);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        double rate = kf;
        for (size_t k = 0; k < numReactant; ++k)
          rate *= cellData[cell][variableIndex(0, k)];
        if (rate > 0.0) {
          for (size_t k = 0; k < numReactant; ++k)
            cellDerivs[cell][variableIndex(0, k)] -= rate;
          for (size_t k = 0; k < numProduct; ++k)
            cellDerivs[cell][variableIndex(1, k)] += rate;
        }
      }
    });
  }
  void derivsWithAbs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
                     Matrix &cellDerivs, Matrix &, Matrix &, Matrix &sdydtCell,
                     Matrix &, Matrix &) override {
    const double kf = parameter(0);
    const size_t numReactant = numVariableIndex(0);
    const size_t numProduct = numVariableIndex(1);
    for (size_t cell = 0; cell < cellDerivs.rows(); ++cell) {
      double rate = kf;
      for (size_t k = 0; k < numReactant; ++k)
        rate *= cellData[cell][variableIndex(0, k)];
      if (rate > 0.0) {
        for (size_t k = 0; k < numReactant; ++k) {
          cellDerivs[cell][variableIndex(0, k)] -= rate;
          sdydtCell[cell][variableIndex(0, k)] += rate;
        }
        for (size_t k = 0; k < numProduct; ++k) {
          cellDerivs[cell][variableIndex(1, k)] += rate;
          sdydtCell[cell][variableIndex(1, k)] += rate;
        }
      }
    }
  }
};
TISSUE_REGISTER_REACTION(MassActionGeneral, "MassAction::General",
                         "MassActionGeneral")

} // namespace
} // namespace tissue
