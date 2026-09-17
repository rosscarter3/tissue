//
// Mass action reactions. Ported from legacy massAction.cc.
//
#include <cmath>
#include <stdexcept>
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

// A -> B + C, first order in A
class MassActionOneToTwo : public Reaction {
public:
  MassActionOneToTwo(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 3 || !i[1].empty())
      throw std::runtime_error(
          "MassAction::OneToTwo: level 0 = reactant, product 1, product 2; "
          "level 1 empty.");
    configure("MassAction::OneToTwo", p, i, 1, {3, 0}, {"k"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t r = variableIndex(0, 0), p1 = variableIndex(0, 1),
                 p2 = variableIndex(0, 2);
    const double k = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t c = b; c < e; ++c) {
        const double fac = k * cellData[c][r];
        cellDerivs[c][p1] += fac;
        cellDerivs[c][p2] += fac;
        cellDerivs[c][r] -= fac;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionOneToTwo, "MassAction::OneToTwo")

// A -> B + C on walls. Wall variables come in pairs (one per side of the
// wall), so each index also drives index+1 with its own rate.
class MassActionOneToTwoWall : public Reaction {
public:
  MassActionOneToTwoWall(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || !i[0].empty() || i[1].size() != 3)
      throw std::runtime_error(
          "MassAction::OneToTwoWall: level 0 empty; level 1 = reactant, "
          "product 1, product 2 (wall variables, each paired with index+1).");
    configure("MassAction::OneToTwoWall", p, i, 1, {0, 3}, {"k"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t r = variableIndex(1, 0), p1 = variableIndex(1, 1),
                 p2 = variableIndex(1, 2);
    const double k = parameter(0);
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        const double fac = k * wallData[w][r];
        const double fac2 = k * wallData[w][r + 1];
        wallDerivs[w][p1] += fac;
        wallDerivs[w][p2] += fac;
        wallDerivs[w][r] -= fac;
        wallDerivs[w][p1 + 1] += fac2;
        wallDerivs[w][p2 + 1] += fac2;
        wallDerivs[w][r + 1] -= fac2;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionOneToTwoWall, "MassAction::OneToTwoWall")

// A + B -> C on walls, paired per side as above.
class MassActionTwoToOneWall : public Reaction {
public:
  MassActionTwoToOneWall(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || !i[0].empty() || i[1].size() != 3)
      throw std::runtime_error(
          "MassAction::TwoToOneWall: level 0 empty; level 1 = reactant 1, "
          "reactant 2, product (wall variables, each paired with index+1).");
    configure("MassAction::TwoToOneWall", p, i, 1, {0, 3}, {"k"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t r1 = variableIndex(1, 0), r2 = variableIndex(1, 1),
                 pr = variableIndex(1, 2);
    const double k = parameter(0);
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        const double fac = k * wallData[w][r1] * wallData[w][r2];
        const double fac2 = k * wallData[w][r1 + 1] * wallData[w][r2 + 1];
        wallDerivs[w][r1] -= fac;
        wallDerivs[w][r2] -= fac;
        wallDerivs[w][pr] += fac;
        wallDerivs[w][r1 + 1] -= fac2;
        wallDerivs[w][r2 + 1] -= fac2;
        wallDerivs[w][pr + 1] += fac2;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionTwoToOneWall, "MassAction::TwoToOneWall")

// General stoichiometry on wall variables: rate = k * prod(reactants), applied
// negatively to every reactant and positively to every product. Legacy gates
// on rate > 0, so a zero or negative rate leaves the wall untouched.
class MassActionGeneralWall : public Reaction {
public:
  MassActionGeneralWall(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2)
      throw std::runtime_error(
          "MassAction::GeneralWall: level 0 = reactants, level 1 = products "
          "(wall variables).");
    configure("MassAction::GeneralWall", p, i, 1, {i[0].size(), i[1].size()},
              {"k"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    if (!numVariableIndex(0))
      return;
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        double rate = parameter(0);
        for (size_t j = 0; j < numVariableIndex(0); ++j)
          rate *= wallData[w][variableIndex(0, j)];
        if (rate <= 0.0)
          continue;
        for (size_t j = 0; j < numVariableIndex(0); ++j)
          wallDerivs[w][variableIndex(0, j)] -= rate;
        for (size_t j = 0; j < numVariableIndex(1); ++j)
          wallDerivs[w][variableIndex(1, j)] += rate;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionGeneralWall, "MassAction::GeneralWall")

// A -> B with a Hill in the product: rate = k A P^n/(K^n + P^n)
class MassActionHillSimple : public Reaction {
public:
  MassActionHillSimple(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "MassAction::HillSimple: level 0 = reactant, level 1 = product.");
    configure("MassAction::HillSimple", p, i, 3, {1, 1},
              {"k", "K_hill", "n_hill"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t r = variableIndex(0, 0), pr = variableIndex(1, 0);
    const double kPow = std::pow(parameter(1), parameter(2));
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t c = b; c < e; ++c) {
        const double power = std::pow(cellData[c][pr], parameter(2));
        double rate = parameter(0) * power / (kPow + power);
        rate *= cellData[c][r];
        if (rate > 0.0) {
          cellDerivs[c][r] -= rate;
          cellDerivs[c][pr] += rate;
        }
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionHillSimple, "MassAction::HillSimple")

// General stoichiometry with enzymes: rate = k * prod(reactants) *
// prod(enzymes). Enzymes (level 2) multiply the rate but are not consumed.
class MassActionGeneralEnzymatic : public Reaction {
public:
  MassActionGeneralEnzymatic(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 3)
      throw std::runtime_error(
          "MassAction::GeneralEnzymatic: level 0 = reactants, level 1 = "
          "products, level 2 = enzymes (not consumed).");
    configure("MassAction::GeneralEnzymatic", p, i, 1,
              {i[0].size(), i[1].size(), i[2].size()}, {"k"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    if (!numVariableIndex(0))
      return;
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t c = b; c < e; ++c) {
        double rate = parameter(0);
        for (size_t j = 0; j < numVariableIndex(0); ++j)
          rate *= cellData[c][variableIndex(0, j)];
        for (size_t j = 0; j < numVariableIndex(2); ++j)
          rate *= cellData[c][variableIndex(2, j)];
        if (rate <= 0.0)
          continue;
        for (size_t j = 0; j < numVariableIndex(0); ++j)
          cellDerivs[c][variableIndex(0, j)] -= rate;
        for (size_t j = 0; j < numVariableIndex(1); ++j)
          cellDerivs[c][variableIndex(1, j)] += rate;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(MassActionGeneralEnzymatic,
                         "MassAction::GeneralEnzymatic")

} // namespace

} // namespace tissue
