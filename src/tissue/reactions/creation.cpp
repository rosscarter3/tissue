//
// Creation reactions: constant and state-dependent production terms for cell
// variables. Ported from legacy creation.cc; identical math, parallel over
// cells (each iteration writes only its own row).
//
#include <cmath>

#include <stdexcept>
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

// Spatial creation reactions share a Hill in some scalar derived from the
// cell centre: dc/dt += k * H, where H is r^n/(K^n+r^n) when the sign
// parameter is positive and K^n/(K^n+r^n) when it is not.
static inline double spatialHill(double value, double kPow, double n,
                                 double sign) {
  const double vPow = std::pow(value, n);
  return sign > 0.0 ? vPow / (kPow + vPow) : kPow / (kPow + vPow);
}

// Distance from the z axis (first two coordinates only).
class CreationSpatialCylinder : public Reaction {
public:
  CreationSpatialCylinder(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::SpatialCylinder", p, i, 4, {1},
              {"k_c", "K_hill", "n_hill", "R_sign"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kPow = std::pow(parameter(1), parameter(2));
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        const auto c = T.cellPosition(cell, vertexData);
        const double r = std::sqrt(c[0] * c[0] + c[1] * c[1]);
        cellDerivs[cell][cIndex] +=
            parameter(0) * spatialHill(r, kPow, parameter(2), parameter(3));
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CreationSpatialCylinder, "Creation::SpatialCylinder",
                         "CreationSpatialCylinder")

// Distance to a ring of radius R about the origin.
class CreationSpatialRing : public Reaction {
public:
  CreationSpatialRing(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::SpatialRing", p, i, 5, {1},
              {"k_c", "K_hill", "R", "n_hill", "R_sign"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double kPow = std::pow(parameter(1), parameter(3));
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        const auto c = T.cellPosition(cell, vertexData);
        double r = 0.0;
        for (double x : c)
          r += x * x;
        r = std::fabs(std::sqrt(r) - parameter(2));
        cellDerivs[cell][cIndex] +=
            parameter(0) * spatialHill(r, kPow, parameter(3), parameter(4));
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CreationSpatialRing, "Creation::SpatialRing",
                         "CreationSpatialRing")

// Hill in one Cartesian coordinate of the cell centre.
class CreationSpatialCoordinate : public Reaction {
public:
  CreationSpatialCoordinate(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::SpatialCoordinate", p, i, 4, {1, 1},
              {"k_c", "K_hill", "n_hill", "X_sign"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double kPow = std::pow(parameter(1), parameter(2));
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        const auto c = T.cellPosition(cell, vertexData);
        cellDerivs[cell][cIndex] +=
            parameter(0) *
            spatialHill(c[xIndex], kPow, parameter(2), parameter(3));
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CreationSpatialCoordinate,
                         "Creation::SpatialCoordinate",
                         "CreationSpatialCoordinate")

// Step in one coordinate: production on one side of the plane x = X only.
class CreationSpatialPlane : public Reaction {
public:
  CreationSpatialPlane(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::SpatialPlane", p, i, 3, {1, 1},
              {"k_c", "X", "X_sign"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double X = parameter(1), k = parameter(0);
    const bool below = parameter(2) < 0.0;
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell) {
        const double x = T.cellPosition(cell, vertexData)[xIndex];
        if (below ? (x <= X) : (x >= X))
          cellDerivs[cell][cIndex] += k;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CreationSpatialPlane, "Creation::SpatialPlane",
                         "CreationSpatialPlane")

// Constant production in a listed set of cells. With a second parameter set,
// the rate is a number of molecules and is divided by cell volume.
//
// NOTE the volume index: legacy divides by the volume of cell `k`, the loop
// counter, rather than of the listed cell `variableIndex(1,k)` it is adding
// to. That is almost certainly a defect, but it is reproduced here so models
// written against legacy behave identically; the two coincide whenever the
// list is 0,1,2,... See tools/port/NOTES.md.
class CreationFromList : public Reaction {
public:
  CreationFromList(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1 && p.size() != 2)
      throw std::runtime_error(
          "Creation::FromList: uses k_c, and optionally a number flag that "
          "divides the rate by cell volume.");
    if (i.size() != 2 || i[0].size() != 1 || i[1].empty())
      throw std::runtime_error(
          "Creation::FromList: level 0 = variable produced; level 1 = the "
          "cell indices to produce in.");
    std::vector<std::string> ids{"k_c"};
    if (p.size() == 2)
      ids.push_back("number_flag");
    configure("Creation::FromList", p, i, p.size(), {1, i[1].size()},
              std::move(ids));
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const bool perVolume = numParameter() == 2 && parameter(1) != 0.0;
    for (size_t k = 0; k < numVariableIndex(1); ++k) {
      const size_t cell = variableIndex(1, k);
      if (cell >= cellDerivs.rows())
        continue;
      cellDerivs[cell][cIndex] +=
          perVolume ? parameter(0) / T.cellVolume(k, vertexData)
                    : parameter(0);
    }
  }
};
TISSUE_REGISTER_REACTION(CreationFromList, "Creation::FromList",
                         "CreationFromList")

// dc/dt += V * k_c * X   (first order in X, scaled by cell volume)
class CreationOneGeometric : public Reaction {
public:
  CreationOneGeometric(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::OneGeometric", p, i, 1, {1, 1}, {"k_c"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t xIndex = variableIndex(1, 0);
    const double k = parameter(0);
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t cell = b; cell < e; ++cell)
        cellDerivs[cell][cIndex] +=
            T.cellVolume(cell, vertexData) * k * cellData[cell][xIndex];
    });
  }
};
TISSUE_REGISTER_REACTION(CreationOneGeometric, "Creation::OneGeometric",
                         "CreationOneGeometric")

// dc/dt += A (1 + sin(2 pi (t/T + phase))), with t accumulated between steps.
// Legacy uses 6.28 rather than 2 pi; kept, since the period it implies is what
// existing models were tuned against.
class CreationSinus : public Reaction {
public:
  CreationSinus(const ParameterList &p, const IndexLevels &i) {
    configure("Creation::Sinus", p, i, 3, {1},
              {"amplitude", "period", "phase"});
  }
  void initiate(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
                Matrix &) override {
    time_ = 0.0;
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &, Matrix &cellDerivs,
              Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    const double v =
        parameter(0) *
        (1.0 + std::sin(6.28 * (time_ / parameter(1) + parameter(2))));
    for (size_t cell = 0; cell < T.numCell(); ++cell)
      cellDerivs[cell][cIndex] += v;
  }
  void update(Tissue &, Matrix &, Matrix &, Matrix &, double h) override {
    time_ += h;
  }

private:
  double time_ = 0.0;
};
TISSUE_REGISTER_REACTION(CreationSinus, "Creation::Sinus", "CreationSinus",
                         "creationSinus")

} // namespace

} // namespace tissue
