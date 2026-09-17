//
// Degradation reactions: first order and state-dependent decay of cell
// variables. Ported from legacy degradation.cc.
//
#include <cmath>
#include <stdexcept>

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

// dc/dt -= k_deg * c * y^n / (K^n + y^n)   (degradation activated by y)
class DegradationHill : public Reaction {
public:
  DegradationHill(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::Hill", p, i, 3, {1, 1}, {"k_deg", "K_hill", "n_hill"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t degraded = variableIndex(0, 0);
    const size_t degrader = variableIndex(1, 0);
    const double kPow = std::pow(parameter(1), parameter(2));
    const double n = parameter(2), k = parameter(0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        const double yPow = std::pow(cellData[cell][degrader], n);
        cellDerivs[cell][degraded] -=
            k * cellData[cell][degraded] * yPow / (kPow + yPow);
      }
    });
  }
};
TISSUE_REGISTER_REACTION(DegradationHill, "Degradation::Hill", "DegradationHill")

// dc/dt -= k_deg * c * prod(activator Hills) * prod(repressor Hills).
// Parameters after the first come in (K, n) pairs, activators then repressors,
// matching the index levels 1 and 2.
class DegradationHillN : public Reaction {
public:
  DegradationHillN(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 3 || i[0].size() != 1)
      throw std::runtime_error(
          "Degradation::HillN: level 0 = variable degraded (one index); "
          "level 1 = activators; level 2 = repressors.");
    if (p.size() != 1 + 2 * (i[1].size() + i[2].size()))
      throw std::runtime_error(
          "Degradation::HillN: uses k_deg followed by a (K, n) pair for each "
          "activator and then each repressor.");
    std::vector<std::string> ids{"k_deg"};
    for (size_t k = 0; k < i[1].size() + i[2].size(); ++k) {
      ids.push_back("K_hill");
      ids.push_back("n_hill");
    }
    configure("Degradation::HillN", p, i, p.size(),
              {1, i[1].size(), i[2].size()}, std::move(ids));
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    parallelFor(cellDerivs.rows(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        double contribution = parameter(0);
        size_t pi = 1;
        for (size_t j = 0; j < numVariableIndex(1); ++j) {
          const double c = std::pow(cellData[cell][variableIndex(1, j)],
                                    parameter(pi + 1));
          contribution *= c / (std::pow(parameter(pi), parameter(pi + 1)) + c);
          pi += 2;
        }
        for (size_t j = 0; j < numVariableIndex(2); ++j) {
          const double c = std::pow(parameter(pi), parameter(pi + 1));
          contribution *= c / (c + std::pow(cellData[cell][variableIndex(2, j)],
                                            parameter(pi + 1)));
          pi += 2;
        }
        cellDerivs[cell][cIndex] -= contribution * cellData[cell][cIndex];
      }
    });
  }
};
TISSUE_REGISTER_REACTION(DegradationHillN, "Degradation::HillN",
                         "DegradationHillN")

// dc/dt -= V * k_d * X * c   (second order, scaled by cell volume)
class DegradationTwoGeometric : public Reaction {
public:
  DegradationTwoGeometric(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::TwoGeometric", p, i, 1, {1, 1}, {"k_d"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double kd = parameter(0);
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] -= T.cellVolume(cell, vertexData) * kd *
                                    cellData[cell][xIndex] *
                                    cellData[cell][cIndex];
    });
  }
};
TISSUE_REGISTER_REACTION(DegradationTwoGeometric, "Degradation::TwoGeometric",
                         "DegradationTwoGeometric")

// dw/dt -= k_cw * w, on wall variables rather than cell variables
class DegradationOneWall : public Reaction {
public:
  DegradationOneWall(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::OneWall", p, i, 1, {1}, {"k_cw"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t wIndex = variableIndex(0, 0);
    const double k = parameter(0);
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w)
        wallDerivs[w][wIndex] -= k * wallData[w][wIndex];
    });
  }
  void derivsWithAbs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
                     Matrix &wallDerivs, Matrix &, Matrix &, Matrix &sdydtWall,
                     Matrix &) override {
    const size_t wIndex = variableIndex(0, 0);
    const double k = parameter(0);
    for (size_t w = 0; w < T.numWall(); ++w) {
      const double value = k * wallData[w][wIndex];
      wallDerivs[w][wIndex] -= value;
      sdydtWall[w][wIndex] += value;
    }
  }
};
TISSUE_REGISTER_REACTION(DegradationOneWall, "Degradation::OneWall",
                         "DegradationOneWall")

// dc/dt -= k_c * c, but only in cells touching the background
class DegradationOneBoundary : public Reaction {
public:
  DegradationOneBoundary(const ParameterList &p, const IndexLevels &i) {
    configure("Degradation::OneBoundary", p, i, 1, {1}, {"k_c"});
  }
  static bool onBoundary(const Tissue &T, size_t cell) {
    for (size_t w : T.cell(cell).walls)
      if (Tissue::isBackground(T.wall(w).otherCell(cell)))
        return true;
    return false;
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double k = parameter(0);
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        if (onBoundary(T, cell))
          cellDerivs[cell][cIndex] -= k * cellData[cell][cIndex];
    });
  }
  void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
                     Matrix &cellDerivs, Matrix &, Matrix &, Matrix &sdydtCell,
                     Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double k = parameter(0);
    for (size_t cell = 0; cell < T.numCell(); ++cell)
      if (onBoundary(T, cell)) {
        const double value = k * cellData[cell][cIndex];
        cellDerivs[cell][cIndex] -= value;
        sdydtCell[cell][cIndex] += value;
      }
  }
};
TISSUE_REGISTER_REACTION(DegradationOneBoundary, "Degradation::OneBoundary",
                         "DegradationOneBoundary")

// dc/dt -= k_d * c, in a listed set of cells only. The second index level is
// a list of cell indices rather than variable indices.
class DegradationOneFromList : public Reaction {
public:
  DegradationOneFromList(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].empty())
      throw std::runtime_error(
          "Degradation::OneFromList: level 0 = variable degraded (one index); "
          "level 1 = the cell indices to degrade in.");
    configure("Degradation::OneFromList", p, i, 1, {1, i[1].size()}, {"k_d"});
  }
  void derivs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kd = parameter(0);
    for (size_t j = 0; j < numVariableIndex(1); ++j) {
      const size_t cell = variableIndex(1, j);
      if (cell < cellDerivs.rows())
        cellDerivs[cell][cIndex] -= kd * cellData[cell][cIndex];
    }
  }
  void derivsWithAbs(Tissue &, Matrix &cellData, Matrix &, Matrix &,
                     Matrix &cellDerivs, Matrix &, Matrix &, Matrix &sdydtCell,
                     Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kd = parameter(0);
    for (size_t j = 0; j < numVariableIndex(1); ++j) {
      const size_t cell = variableIndex(1, j);
      if (cell < cellDerivs.rows()) {
        const double value = kd * cellData[cell][cIndex];
        cellDerivs[cell][cIndex] -= value;
        sdydtCell[cell][cIndex] += value;
      }
    }
  }
};
TISSUE_REGISTER_REACTION(DegradationOneFromList, "Degradation::OneFromList",
                         "DegradationOneFromList")

} // namespace

} // namespace tissue
