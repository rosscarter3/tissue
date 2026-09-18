//
// Auxin/PIN network models. Ported from legacy network.cc, which holds 41 of
// them in 5700 lines - families of near-identical reactions differing in one
// production term, one Hill factor, or whether a flux is weighted by geometry.
//
// Nearly all share one shape: a per-cell block of production and degradation
// ODEs, then a polarised transport block that walks the cell's interior walls,
// builds a normalising sum over the neighbours, and moves auxin out of the
// cell in proportion to the PIN it allocates to each membrane. The traversal
// and the sum are here as helpers; each reaction supplies its own chemistry.
//
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// A wall with a real cell on both sides.
inline bool interiorWall(const Tissue &T, size_t w) {
  return !Tissue::isBackground(T.wall(w).cell1) &&
         !Tissue::isBackground(T.wall(w).cell2);
}

// Visits every interior wall of cell c in the cell's own wall order, which is
// the order legacy accumulates in. `side` is 0 when c is the wall's cell1 and
// 1 otherwise, so `paired + side` addresses this cell's face of a paired wall
// variable.
template <class Fn>
void forEachInteriorWall(const Tissue &T, size_t c, Fn fn) {
  const CellTopo &cell = T.cell(c);
  for (size_t n = 0; n < cell.numWall(); ++n) {
    const size_t w = cell.walls[n];
    if (!interiorWall(T, w))
      continue;
    const Wall &wall = T.wall(w);
    const size_t side = wall.cell1 == c ? 0 : 1;
    fn(n, w, wall.otherCell(c), side);
  }
}

// True when the cell has at least one wall on the tissue boundary. Legacy
// spells this `T.cell(i).isNeighbor(T.background())` and uses it to mark the
// outer cell layer (L1).
inline bool onBoundary(const Tissue &T, size_t c) {
  for (size_t w : T.cell(c).walls)
    if (Tissue::isBackground(T.wall(w).otherCell(c)))
      return true;
  return false;
}

// Sum of a cell variable over the neighbours across interior walls.
inline double neighbourSum(const Tissue &T, const Matrix &cellData, size_t c,
                           size_t index) {
  double sum = 0.0;
  forEachInteriorWall(T, c, [&](size_t, size_t, size_t neigh, size_t) {
    sum += cellData[neigh][index];
  });
  return sum;
}

// --- the AuxinModelSimple1 family -------------------------------------------
//
// Auxin is produced at a rate set by a membrane-bound marker M (which is
// itself produced only in cells touching the boundary, so it marks the outer
// layer), degraded, and exported to neighbours in proportion to the PIN each
// membrane carries. PIN allocation to a membrane is the cell's PIN times the
// neighbour's X, normalised by the sum of X over all neighbours plus p_pol -
// so a cell polarises toward whichever neighbour carries the most X, and X is
// itself driven by auxin.
//
// kGeometric weights the flux by wall length / cell volume (AuxinModel1);
// without it the flux is per-membrane and geometry-free (AuxinModelSimple1).
template <bool kGeometric>
class AuxinModelSimple1Base : public Reaction {
public:
  AuxinModelSimple1Base(const ParameterList &p, const IndexLevels &i) {
    if ((i.size() != 1 && i.size() != 2) || i[0].size() != 4 ||
        (i.size() == 2 && i[1].size() != 1))
      throw std::runtime_error(
          "AuxinModelSimple1: four cell variable indices (auxin, PIN, X, M); "
          "an optional second level names a paired wall variable to store the "
          "membrane PIN in.");
    std::vector<size_t> shape{4};
    if (i.size() == 2)
      shape.push_back(1);
    configure(kGeometric ? "AuxinModel1" : "AuxinModelSimple1", p, i, 12, shape,
              {"p_auxin(M)", "p_auxin", "d_auxin", "p_pol", "T_auxin",
               "D_auxin", "p_pin", "d_pin", "p_X", "d_X", "p_M", "d_M"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &,
              Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t xI = variableIndex(0, 2), mI = variableIndex(0, 3);
    const bool storePin = numVariableIndexLevel() == 2;
    const size_t pinWallIndex = storePin ? variableIndex(1, 0) : 0;

    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) * cellData[c][mI] + parameter(1) -
                           parameter(2) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) - parameter(7) * cellData[c][pI];
      cellDerivs[c][xI] +=
          parameter(8) * cellData[c][aI] - parameter(9) * cellData[c][xI];
      cellDerivs[c][mI] -= parameter(11) * cellData[c][mI];
      if (onBoundary(T, c))
        cellDerivs[c][mI] += parameter(10);

      const double sum = neighbourSum(T, cellData, c, xI) + parameter(3);
      forEachInteriorWall(
          T, c, [&](size_t, size_t w, size_t neigh, size_t side) {
            // A cell with no X anywhere around it spreads PIN evenly, which
            // legacy writes as a rate of exactly 1 rather than 1/n.
            const double polRate =
                sum != 0.0 ? cellData[c][pI] * cellData[neigh][xI] / sum : 1.0;
            if (storePin)
              wallData[w][pinWallIndex + side] = polRate;
            const double rate = parameter(4) * polRate + parameter(5);
            if (kGeometric) {
              const double flux = T.wallLengthFromVertices(w, vertexData) *
                                  rate * cellData[c][aI];
              cellDerivs[c][aI] -= flux / T.cellVolume(c, vertexData);
              cellDerivs[neigh][aI] += flux / T.cellVolume(neigh, vertexData);
            } else {
              cellDerivs[c][aI] -= rate * cellData[c][aI];
              cellDerivs[neigh][aI] += rate * cellData[c][aI];
            }
          });
    }
  }
};
using AuxinModelSimple1 = AuxinModelSimple1Base<false>;
using AuxinModel1 = AuxinModelSimple1Base<true>;
TISSUE_REGISTER_REACTION(AuxinModelSimple1, "AuxinModelSimple1")
TISSUE_REGISTER_REACTION(AuxinModel1, "AuxinModel1")

// As AuxinModelSimple1, but the quantity the cell polarises on lives in the
// *wall* rather than in the neighbouring cell, and every wall contributes to
// the normalising sum - boundary walls included, unlike the transport itself.
//
// Legacy shifts the whole set up by the most negative value before summing,
// so a negative wall signal cannot make the denominator smaller than the
// numerator. That correction is applied inside the accumulation loop rather
// than after it, so it only sees the walls visited so far; reproduced as
// written.
class AuxinModelSimple1Wall : public Reaction {
public:
  AuxinModelSimple1Wall(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 3)
      throw std::runtime_error("AuxinModelSimple1Wall: three indices (cell "
                               "auxin, cell PIN, wall signal).");
    configure("AuxinModelSimple1Wall", p, i, 7, {3},
              {"p_auxin", "d_auxin", "p_pol", "T_auxin", "D_auxin", "p_pin",
               "d_pin"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t xI = variableIndex(0, 2);
    std::vector<double> pin;
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(5) - parameter(6) * cellData[c][pI];
      const CellTopo &cell = T.cell(c);
      const size_t numWalls = cell.numWall();
      pin.assign(numWalls, 0.0);
      double sum = 0.0, minPin = 0.0;
      for (size_t n = 0; n < numWalls; ++n) {
        sum += pin[n] = wallData[cell.walls[n]][xI];
        if (pin[n] < minPin)
          minPin = pin[n];
        if (minPin < 0.0) {
          sum += numWalls * minPin;
          for (size_t k = 0; k < numWalls; ++k)
            pin[k] += minPin;
        }
      }
      sum += parameter(2);
      forEachInteriorWall(T, c, [&](size_t n, size_t, size_t neigh, size_t) {
        double polRate = 0.0;
        if (sum > 0.0)
          polRate = cellData[c][pI] * pin[n] / sum;
        const double rate = parameter(3) * polRate + parameter(4);
        cellDerivs[c][aI] -= rate * cellData[c][aI];
        cellDerivs[neigh][aI] += rate * cellData[c][aI];
      });
    }
  }
};
TISSUE_REGISTER_REACTION(AuxinModelSimple1Wall, "AuxinModelSimple1Wall")

// Transport only, with no production or degradation: passive diffusion plus a
// PIN-driven active term that saturates in both the donor's auxin and the two
// cells' AUX. The polarisation signal is a Hill function of the neighbour's X
// rather than X itself, and boundary walls still contribute the baseline k1 to
// the normalising sum even though nothing moves across them.
class AuxinTransportCellCellNoGeometry : public Reaction {
public:
  AuxinTransportCellCellNoGeometry(const ParameterList &p,
                                   const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 4)
      throw std::runtime_error("AuxinTransportCellCellNoGeometry: four cell "
                               "variable indices (auxin, PIN, AUX, X).");
    configure("AuxinTransportCellCellNoGeometry", p, i, 7, {4},
              {"d", "T", "k1", "k2", "K_H", "n", "K_M"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0), PI = variableIndex(0, 1);
    const size_t AI = variableIndex(0, 2), xI = variableIndex(0, 3);
    const double kPow = std::pow(parameter(4), parameter(5));
    std::vector<double> pin;
    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      pin.assign(cell.numWall(), 0.0);
      double sum = 1.0;
      for (size_t n = 0; n < cell.numWall(); ++n) {
        const size_t w = cell.walls[n];
        if (interiorWall(T, w)) {
          const double x =
              std::pow(cellData[T.wall(w).otherCell(c)][xI], parameter(5));
          sum += pin[n] = parameter(2) + parameter(3) * x / (kPow + x);
        } else {
          sum += pin[n] = parameter(2);
        }
      }
      forEachInteriorWall(T, c, [&](size_t n, size_t, size_t j, size_t) {
        const double passive = parameter(0) * cellData[c][aI];
        const double Pij = cellData[c][PI] * pin[n] / sum;
        const double active = Pij * parameter(1) * cellData[j][AI] *
                              cellData[c][aI] /
                              ((parameter(6) + cellData[c][aI]) *
                               (cellData[c][AI] + cellData[j][AI]));
        cellDerivs[c][aI] -= passive + active;
        cellDerivs[j][aI] += passive + active;
      });
    }
  }
};
TISSUE_REGISTER_REACTION(AuxinTransportCellCellNoGeometry,
                         "AuxinTransportCellCellNoGeometry")


// --- the SimpleROPModel family ----------------------------------------------
//
// Seven variants of one ROP/PIN polarisation model. PIN cycles between the
// cell's cytoplasmic pool and each membrane; what makes a membrane accumulate
// PIN is a Hill function of the PIN on the *opposite* face of the same wall,
// which is the positive feedback that polarises neighbouring cells against
// each other. The first three carry auxin in a wall compartment; the last
// four move it cell to cell directly and make PIN production saturate in
// auxin.
//
// They share the traversal and differ only in a term or two each, so the
// per-membrane chemistry is written out in each rather than hidden behind
// flags - the differences are the point.
//
// Common index layout: level 0 = (cell auxin, cell PIN); level 1 = (paired
// wall auxin, paired membrane PIN), with the wall auxin unused in models 4-7.

// Hill function of the opposite membrane's PIN.
inline double ropHill(double pOther, double K, double n) {
  const double x = std::pow(pOther, n);
  return x / (std::pow(K, n) + x);
}

#define ROP_CHECK(NAME, NPARAM, NLEVEL1)                                       \
  if (i.size() != 2 || i[0].size() != 2 || i[1].size() != NLEVEL1)             \
    throw std::runtime_error(NAME ": level 0 = (cell auxin, cell PIN); "       \
                                  "level 1 = paired wall variables.");         \
  configure(NAME, p, i, NPARAM, {2, NLEVEL1},                                  \
            {"c_IAA", "d_IAA", "p_IAAH(in)", "p_IAAH(out)", "p_IAA-",          \
             "D_IAA", "c_PIN", "d_PIN", "endo_PIN", "exo_PIN", "K_hill",       \
             "n_hill", "endo_PIN_back", "p_13", "p_14", "p_15"});

// Auxin lives in the wall: it is secreted into the wall compartment, diffuses
// across to the opposite face, and PIN is removed from the membrane in
// proportion to the cell's PIN pool times the wall auxin on that face.
class SimpleROPModel : public Reaction {
public:
  SimpleROPModel(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel", 13, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t awI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) - parameter(7) * cellData[c][pI];
      forEachInteriorWall(T, c, [&](size_t, size_t w, size_t, size_t side) {
        const size_t own = side, other = 1 - side;
        double fac = parameter(3) * cellData[c][aI] -
                     parameter(2) * wallData[w][awI + own] +
                     parameter(4) * cellData[c][aI] * wallData[w][pwI + own];
        wallDerivs[w][awI + own] += fac;
        cellDerivs[c][aI] -= fac;
        wallDerivs[w][awI + own] -= parameter(5) * wallData[w][awI + own];
        wallDerivs[w][awI + other] += parameter(5) * wallData[w][awI + own];
        fac = parameter(12) * wallData[w][pwI + own] +
              parameter(8) * wallData[w][pwI + own] *
                  ropHill(wallData[w][pwI + other], parameter(10),
                          parameter(11)) -
              parameter(9) * cellData[c][pI] * wallData[w][awI + own];
        wallDerivs[w][pwI + own] -= fac;
        cellDerivs[c][pI] += fac;
      });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel, "SimpleROPModel")

// SimpleROPModel with the PIN removal no longer proportional to the wall
// auxin on that face - the only difference between the two.
class SimpleROPModel2 : public Reaction {
public:
  SimpleROPModel2(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel2", 13, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t awI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) - parameter(7) * cellData[c][pI];
      forEachInteriorWall(T, c, [&](size_t, size_t w, size_t, size_t side) {
        const size_t own = side, other = 1 - side;
        double fac = parameter(3) * cellData[c][aI] -
                     parameter(2) * wallData[w][awI + own] +
                     parameter(4) * cellData[c][aI] * wallData[w][pwI + own];
        wallDerivs[w][awI + own] += fac;
        cellDerivs[c][aI] -= fac;
        wallDerivs[w][awI + own] -= parameter(5) * wallData[w][awI + own];
        wallDerivs[w][awI + other] += parameter(5) * wallData[w][awI + own];
        fac = parameter(12) * wallData[w][pwI + own] +
              parameter(8) * wallData[w][pwI + own] *
                  ropHill(wallData[w][pwI + other], parameter(10),
                          parameter(11)) -
              parameter(9) * cellData[c][pI];
        wallDerivs[w][pwI + own] -= fac;
        cellDerivs[c][pI] += fac;
      });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel2, "SimpleROPModel2")

// As SimpleROPModel2, with the constant PIN return term written last rather
// than first. Algebraically the same reaction; kept separate because legacy
// registers both names.
class SimpleROPModel3 : public Reaction {
public:
  SimpleROPModel3(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel3", 13, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t awI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) - parameter(7) * cellData[c][pI];
      forEachInteriorWall(T, c, [&](size_t, size_t w, size_t, size_t side) {
        const size_t own = side, other = 1 - side;
        double fac = parameter(3) * cellData[c][aI] -
                     parameter(2) * wallData[w][awI + own] +
                     parameter(4) * cellData[c][aI] * wallData[w][pwI + own];
        wallDerivs[w][awI + own] += fac;
        cellDerivs[c][aI] -= fac;
        wallDerivs[w][awI + own] -= parameter(5) * wallData[w][awI + own];
        wallDerivs[w][awI + other] += parameter(5) * wallData[w][awI + own];
        fac = parameter(8) * wallData[w][pwI + own] *
                  ropHill(wallData[w][pwI + other], parameter(10),
                          parameter(11)) -
              parameter(9) * cellData[c][pI] +
              parameter(12) * wallData[w][pwI + own];
        wallDerivs[w][pwI + own] -= fac;
        cellDerivs[c][pI] += fac;
      });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel3, "SimpleROPModel3")

// From here the wall auxin compartment is dropped: auxin moves straight from
// cell to cell, carried by the membrane PIN, and PIN production saturates in
// the cell's own auxin. Level 1 still takes two indices, but the first (the
// wall auxin) is no longer read.
class SimpleROPModel4 : public Reaction {
public:
  SimpleROPModel4(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel4", 13, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) * cellData[c][aI] /
                               (parameter(4) + cellData[c][aI]) -
                           parameter(7) * cellData[c][pI];
      forEachInteriorWall(
          T, c, [&](size_t, size_t w, size_t neigh, size_t side) {
            const size_t own = side, other = 1 - side;
            double fac =
                parameter(2) * cellData[c][aI] +
                parameter(3) * cellData[c][aI] * wallData[w][pwI + own];
            cellDerivs[c][aI] -= fac;
            cellDerivs[neigh][aI] += fac;
            fac = parameter(8) * wallData[w][pwI + own] *
                      ropHill(wallData[w][pwI + other], parameter(10),
                              parameter(11)) -
                  parameter(9) * cellData[c][pI] +
                  parameter(12) * wallData[w][pwI + own];
            wallDerivs[w][pwI + own] -= fac;
            cellDerivs[c][pI] += fac;
          });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel4, "SimpleROPModel4")

// SimpleROPModel4 with a second, self-limiting removal term: PIN is taken off
// a membrane faster the more PIN that membrane already carries.
class SimpleROPModel5 : public Reaction {
public:
  SimpleROPModel5(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel5", 16, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) * cellData[c][aI] /
                               (parameter(4) + cellData[c][aI]) -
                           parameter(7) * cellData[c][pI];
      forEachInteriorWall(
          T, c, [&](size_t, size_t w, size_t neigh, size_t side) {
            const size_t own = side, other = 1 - side;
            double fac =
                parameter(2) * cellData[c][aI] +
                parameter(3) * cellData[c][aI] * wallData[w][pwI + own];
            cellDerivs[c][aI] -= fac;
            cellDerivs[neigh][aI] += fac;
            fac = parameter(8) * wallData[w][pwI + own] *
                      ropHill(wallData[w][pwI + other], parameter(10),
                              parameter(11)) -
                  parameter(9) * cellData[c][pI] +
                  parameter(12) * wallData[w][pwI + own] -
                  parameter(13) * cellData[c][pI] *
                      ropHill(wallData[w][pwI + own], parameter(15),
                              parameter(14));
            wallDerivs[w][pwI + own] -= fac;
            cellDerivs[c][pI] += fac;
          });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel5, "SimpleROPModel5")

// SimpleROPModel4 that also cycles PIN on the cell's *boundary* membranes,
// where there is no neighbour to transport to and no opposite face to feed
// back from. Legacy reads the cell1 face there whichever side the cell is on,
// which is reproduced.
class SimpleROPModel6 : public Reaction {
public:
  SimpleROPModel6(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel6", 13, 2)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 1);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) * cellData[c][aI] /
                               (parameter(4) + cellData[c][aI]) -
                           parameter(7) * cellData[c][pI];
      const CellTopo &cell = T.cell(c);
      for (size_t n = 0; n < cell.numWall(); ++n) {
        const size_t w = cell.walls[n];
        if (!interiorWall(T, w)) { // boundary membrane: cycling only
          const double fac = -parameter(9) * cellData[c][pI] +
                             parameter(12) * wallData[w][pwI];
          wallDerivs[w][pwI] -= fac;
          cellDerivs[c][pI] += fac;
          continue;
        }
        const size_t neigh = T.wall(w).otherCell(c);
        const size_t own = T.wall(w).cell1 == c ? 0 : 1, other = 1 - own;
        double fac = parameter(2) * cellData[c][aI] +
                     parameter(3) * cellData[c][aI] * wallData[w][pwI + own];
        cellDerivs[c][aI] -= fac;
        cellDerivs[neigh][aI] += fac;
        fac = parameter(8) * wallData[w][pwI + own] *
                  ropHill(wallData[w][pwI + other], parameter(10),
                          parameter(11)) -
              parameter(9) * cellData[c][pI] +
              parameter(12) * wallData[w][pwI + own];
        wallDerivs[w][pwI + own] -= fac;
        cellDerivs[c][pI] += fac;
      }
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel6, "SimpleROPModel6")

// SimpleROPModel4 gated by a wall marker: walls flagged 1 run the polarising
// model, walls flagged 0 instead leak auxin and PIN at fixed rates.
//
// The leak branch has a legacy slip that is harmless and reproduced: on the
// cell1 side it sets the neighbour to `cell1`, which is the cell itself, so
// the two derivative contributions cancel exactly and nothing moves. Only the
// cell2 side of a flagged-0 wall actually leaks.
class SimpleROPModel7 : public Reaction {
public:
  SimpleROPModel7(const ParameterList &p, const IndexLevels &i) {
    ROP_CHECK("SimpleROPModel7", 15, 3)
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 1), mwI = variableIndex(1, 2);
    for (size_t c = 0; c < T.numCell(); ++c) {
      cellDerivs[c][aI] += parameter(0) - parameter(1) * cellData[c][aI];
      cellDerivs[c][pI] += parameter(6) * cellData[c][aI] /
                               (parameter(4) + cellData[c][aI]) -
                           parameter(7) * cellData[c][pI];
      forEachInteriorWall(
          T, c, [&](size_t, size_t w, size_t neigh, size_t side) {
            const size_t own = side, other = 1 - side;
            if (wallData[w][mwI] == 1) {
              double fac =
                  parameter(2) * cellData[c][aI] +
                  parameter(3) * cellData[c][aI] * wallData[w][pwI + own];
              cellDerivs[c][aI] -= fac;
              cellDerivs[neigh][aI] += fac;
              fac = parameter(8) * wallData[w][pwI + own] *
                        ropHill(wallData[w][pwI + other], parameter(10),
                                parameter(11)) -
                    parameter(9) * cellData[c][pI] +
                    parameter(12) * wallData[w][pwI + own];
              wallDerivs[w][pwI + own] -= fac;
              cellDerivs[c][pI] += fac;
            } else if (wallData[w][mwI] == 0) {
              // Legacy sends the cell1 side's leak to cell1 - itself - so it
              // cancels; only the cell2 side moves anything.
              const size_t target = side == 1 ? neigh : c;
              const double fa = parameter(13) * cellData[c][aI];
              const double fp = parameter(14) * cellData[c][pI];
              cellDerivs[c][aI] -= fa;
              cellDerivs[c][pI] -= fp;
              cellDerivs[target][aI] += fa;
              cellDerivs[target][pI] += fp;
            }
          });
    }
  }
};
TISSUE_REGISTER_REACTION(SimpleROPModel7, "SimpleROPModel7")

} // namespace
} // namespace tissue
