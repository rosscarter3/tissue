//
// Transport reactions between neighboring cells. Ported from legacy
// transport.cc.
//
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Plain concentration-difference diffusion across shared (non-boundary)
// walls; no geometric factors ("Simple").
class DiffusionSimple : public Reaction {
public:
  DiffusionSimple(const ParameterList &p, const IndexLevels &i) {
    configure("DiffusionSimple", p, i, 1, {1}, {"p_0"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const double p0 = parameter(0);
    // Legacy iterates cell-wall pairs with an i<neighbor guard; iterating
    // interior walls once is equivalent (each interior wall borders exactly
    // the two cells) and parallel-friendly.
    parallelScatter1(T.numWall(), cellDerivs,
                     [&](size_t b, size_t e, Matrix &out) {
                       for (size_t w = b; w < e; ++w) {
                         const Wall &wall = T.wall(w);
                         if (Tissue::isBackground(wall.cell1) ||
                             Tissue::isBackground(wall.cell2))
                           continue;
                         size_t i = std::min(wall.cell1, wall.cell2);
                         size_t neigh = std::max(wall.cell1, wall.cell2);
                         double flux = p0 * (cellData[i][aI] - cellData[neigh][aI]);
                         out[i][aI] -= flux;
                         out[neigh][aI] += flux;
                       }
                     });
  }
  void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &wallData,
                     Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                     Matrix &vertexDerivs, Matrix &, Matrix &,
                     Matrix &) override {
    // Legacy adds zero noise for diffusion (documented deterministic).
    derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs,
           vertexDerivs);
  }
};
TISSUE_REGISTER_REACTION(DiffusionSimple, "DiffusionSimple", "Diffusion::Simple")

// Walls carry *paired* variables throughout this file: index `k` is the face
// belonging to wall.cell1 and `k+1` the face belonging to wall.cell2, which is
// how legacy represents the two membranes of one wall. Legacy reaches them by
// walking each cell's walls and asking which side it is on; iterating walls
// directly and reading cell1/cell2 is the same traversal without the lookup.

// Diffusion of a membrane-bound species around the inside of each cell: each
// wall face loses to, and gains from, the two faces adjacent to it in the
// cell's wall cycle.
class MembraneDiffusionSimple : public Reaction {
public:
  MembraneDiffusionSimple(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error("Diffusion::MembraneSimple: one index (the "
                               "paired wall variable that diffuses).");
    configure("Diffusion::MembraneSimple", p, i, 1, {1}, {"D_mem"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t pwI = variableIndex(0, 0);
    const double D = parameter(0);
    // Which face of wall w belongs to cell c.
    auto face = [&](size_t w, size_t c) {
      return T.wall(w).cell1 == c ? pwI : pwI + 1;
    };
    parallelScatter1(
        T.numCell(), wallDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numWall();
            for (size_t k = 0; k < n; ++k) {
              const size_t w = cell.walls[k];
              const double fac = D * wallData[w][face(w, c)];
              out[w][face(w, c)] -= 2.0 * fac;
              const size_t wNext = cell.walls[(k + 1) % n];
              out[wNext][face(wNext, c)] += fac;
              const size_t wPrev = cell.walls[(k + n - 1) % n];
              out[wPrev][face(wPrev, c)] += fac;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(MembraneDiffusionSimple, "Diffusion::MembraneSimple",
                         "MembraneDiffusionSimple")

// Diffusion of a cell variable whose mobility is gated by a second cell
// variable: the flux is set by the difference in the *product* of the two.
class DiffusionSimpleOne : public Reaction {
public:
  DiffusionSimpleOne(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("Diffusion::SimpleOne: two levels of one index "
                               "each (the diffusing variable, then the gate).");
    configure("Diffusion::SimpleOne", p, i, 1, {1, 1}, {"p_0"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t hI = variableIndex(1, 0);
    const double p0 = parameter(0);
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const size_t i = std::min(wall.cell1, wall.cell2);
            const size_t neigh = std::max(wall.cell1, wall.cell2);
            const double fac = p0 * (cellData[i][aI] * cellData[i][hI] -
                                     cellData[neigh][aI] * cellData[neigh][hI]);
            out[i][aI] -= fac;
            out[neigh][aI] += fac;
          }
        });
  }
};
TISSUE_REGISTER_REACTION(DiffusionSimpleOne, "Diffusion::SimpleOne",
                         "DiffusionSimpleOne")

// Diffusion whose conductance is a wall variable that adapts to the flux it
// carries - the standard flux-reinforcement rule for vein patterning:
//   dC/dt = p_1 (|flux|^p_2 / C^(p_3+1) - p_4) C
// A wall at zero conductance stays there (legacy gates the update on C > 0),
// so the conductance field has to be seeded non-zero to pattern at all.
class DiffusionConductiveSimple : public Reaction {
public:
  DiffusionConductiveSimple(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "Diffusion::ConductiveSimple: two levels of one index each (the "
          "diffusing cell variable, then the wall conductivity).");
    configure("Diffusion::ConductiveSimple", p, i, 5, {1, 1},
              {"p_0", "p_1", "p_2", "p_3", "p_4"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t cI = variableIndex(0, 0);
    const size_t CI = variableIndex(1, 0);
    parallelScatter2(
        T.numWall(), cellDerivs, wallDerivs,
        [&](size_t b, size_t e, Matrix &outCell, Matrix &outWall) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const size_t i = std::min(wall.cell1, wall.cell2);
            const size_t neigh = std::max(wall.cell1, wall.cell2);
            const double conductance = wallData[w][CI];
            double flux =
                parameter(0) * conductance * (cellData[i][cI] - cellData[neigh][cI]);
            outCell[i][cI] -= flux;
            outCell[neigh][cI] += flux;
            if (conductance > 0.0) {
              flux = std::fabs(flux);
              outWall[w][CI] +=
                  parameter(1) *
                  ((std::pow(flux, parameter(2)) /
                    std::pow(conductance, parameter(3) + 1)) -
                   parameter(4)) *
                  conductance;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(DiffusionConductiveSimple,
                         "Diffusion::ConductiveSimple",
                         "DiffusionConductiveSimple")

// Diffusion with the geometric factors of a finite-volume scheme: flux scales
// with the shared wall length and falls with the distance between cell
// centres, and is divided by the receiving cell's area to give a concentration
// rate.
//
// Two legacy behaviours to be aware of:
//
// FIX: legacy takes the two cell centres from `Cell::positionFromVertex()`,
// the no-argument overload, which reads each Vertex's *cached* position rather
// than the live `vertexData` passed into derivs. That cache is only refreshed
// by BaseSolver::setTissueVariables(), i.e. at print points, so on a moving
// mesh the centre-to-centre distance lagged the geometry by up to a whole
// print interval while the wall length and cell area in the same expression
// were current. v2 reads all three from `vertexData`. On a static mesh the two
// agree exactly.
//
// KEPT: legacy visits each interior wall once from each side without an
// ordering guard, so every flux is applied twice and the effective diffusion
// constant is 2*p_0. Models are calibrated against that, so the factor of two
// is preserved here rather than silently halving existing parameters.
class Diffusion2d : public Reaction {
public:
  Diffusion2d(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "Diffusion::2D: one index (the diffusing cell variable).");
    configure("Diffusion::2D", p, i, 1, {1}, {"p_0"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const double p0 = parameter(0);
    const size_t dim = vertexData.cols();
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const size_t i = wall.cell1;
            const size_t neigh = wall.cell2;
            const Vec3 pi = T.cellPosition(i, vertexData);
            const Vec3 pn = T.cellPosition(neigh, vertexData);
            double distance = 0.0;
            for (size_t d = 0; d < dim; ++d)
              distance += (pn[d] - pi[d]) * (pn[d] - pi[d]);
            distance = std::sqrt(distance);
            const double contactLength =
                T.wallLengthFromVertices(w, vertexData);
            const double drive =
                p0 * contactLength * (cellData[i][aI] - cellData[neigh][aI]) /
                distance;
            // Twice, as legacy does - see the note above.
            out[i][aI] -= 2.0 * drive / T.cellVolume(i, vertexData);
            out[neigh][aI] += 2.0 * drive / T.cellVolume(neigh, vertexData);
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Diffusion2d, "Diffusion::2D", "Diffusion2D",
                         "Diffusion::2d", "Diffusion2d")

// PIN-driven efflux: each cell loses auxin across each membrane at a rate set
// by that membrane's own PIN level. Purely directional - there is no return
// term, so the two faces of a wall act independently.
class ActiveTransportCellEfflux : public Reaction {
public:
  ActiveTransportCellEfflux(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("ActiveTransportCellEfflux: one cell variable "
                               "(auxin) and one paired wall variable (PIN).");
    configure("ActiveTransportCellEfflux", p, i, 1, {1, 1}, {"T"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t pwI = variableIndex(1, 0);
    const double rate = parameter(0);
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const double out1 = rate * cellData[wall.cell1][aI] * wallData[w][pwI];
            out[wall.cell1][aI] -= out1;
            out[wall.cell2][aI] += out1;
            const double out2 =
                rate * cellData[wall.cell2][aI] * wallData[w][pwI + 1];
            out[wall.cell2][aI] -= out2;
            out[wall.cell1][aI] += out2;
          }
        });
  }
};
TISSUE_REGISTER_REACTION(ActiveTransportCellEfflux, "ActiveTransportCellEfflux")

// ActiveTransportCellEfflux with the carrier saturating in the donor cell's
// auxin: T A/(A + K_T) P instead of T A P.
class ActiveTransportCellEffluxMM : public Reaction {
public:
  ActiveTransportCellEffluxMM(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("ActiveTransportCellEffluxMM: one cell variable "
                               "(auxin) and one paired wall variable (PIN).");
    configure("ActiveTransportCellEffluxMM", p, i, 2, {1, 1}, {"T", "K_T"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t pwI = variableIndex(1, 0);
    const double rate = parameter(0), K = parameter(1);
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const double a1 = cellData[wall.cell1][aI];
            const double out1 = rate * (a1 / (a1 + K)) * wallData[w][pwI];
            out[wall.cell1][aI] -= out1;
            out[wall.cell2][aI] += out1;
            const double a2 = cellData[wall.cell2][aI];
            const double out2 = rate * (a2 / (a2 + K)) * wallData[w][pwI + 1];
            out[wall.cell2][aI] -= out2;
            out[wall.cell1][aI] += out2;
          }
        });
  }
};
TISSUE_REGISTER_REACTION(ActiveTransportCellEffluxMM,
                         "ActiveTransportCellEffluxMM")

// Passive diffusion plus PIN-driven efflux in one net cell-to-cell flux, with
// no explicit wall compartment. An optional third index writes the flux into a
// paired wall variable for visualisation: the positive part goes to the face
// of the cell it leaves.
//
// That write goes into `wallData`, not `wallDerivs` - it is a report, not a
// state derivative, and it happens once per derivative evaluation, so what
// ends up in the output depends on how many stages the solver takes and which
// one ran last. Legacy does the same; it is a diagnostic, not a variable to
// integrate or feed back on.
class DiffusionActiveTransportCell : public Reaction {
public:
  DiffusionActiveTransportCell(const ParameterList &p, const IndexLevels &i) {
    const bool ok = (i.size() == 2 || i.size() == 3) && i[0].size() == 1 &&
                    i[1].size() == 1 && (i.size() == 2 || i[2].size() == 1);
    if (!ok)
      throw std::runtime_error(
          "DiffusionActiveTransportCell: one cell variable (auxin) and one "
          "paired wall variable (PIN); an optional third level names a paired "
          "wall variable to record the flux in.");
    std::vector<size_t> shape{1, 1};
    if (i.size() == 3)
      shape.push_back(1);
    configure("DiffusionActiveTransportCell", p, i, 2, shape, {"D", "T"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t pwI = variableIndex(1, 0);
    const bool saveFlux = numVariableIndexLevel() == 3;
    const size_t fI = saveFlux ? variableIndex(2, 0) : 0;
    const double D = parameter(0), rate = parameter(1);
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const size_t c1 = wall.cell1, c2 = wall.cell2;
            const double fac = (D + rate * wallData[w][pwI]) * cellData[c1][aI] -
                               (D + rate * wallData[w][pwI + 1]) * cellData[c2][aI];
            out[c1][aI] -= fac;
            out[c2][aI] += fac;
            if (saveFlux) {
              wallData[w][fI] = fac >= 0.0 ? fac : 0.0;
              wallData[w][fI + 1] = fac >= 0.0 ? 0.0 : -fac;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(DiffusionActiveTransportCell,
                         "DiffusionActiveTransportCell")

// The full two-compartment picture: auxin lives in the wall as well as in the
// cell. Each membrane exchanges with its own cell (protonated influx p_0,
// AUX-mediated influx p_1, passive efflux p_2, PIN-mediated efflux p_3) and
// the two wall faces exchange with each other by diffusion p_4.
class ActiveTransportWall : public Reaction {
public:
  ActiveTransportWall(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 2)
      throw std::runtime_error(
          "ActiveTransportWall: two cell variables (auxin, AUX) and two paired "
          "wall variables (auxin, PIN).");
    configure("ActiveTransportWall", p, i, 5, {2, 2},
              {"p_IAAH(in)", "p_IAA-(in,AUX)", "p_IAAH(out)",
               "p_IAA-(out,PIN)", "D_IAA"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t auxI = variableIndex(0, 1);
    const size_t awI = variableIndex(1, 0);
    const size_t pwI = variableIndex(1, 1);
    // Cell order, as legacy: each cell settles its own membrane before the
    // neighbour touches the shared wall's other face.
    parallelScatter2(
        T.numCell(), cellDerivs, wallDerivs,
        [&](size_t b, size_t e, Matrix &outCell, Matrix &outWall) {
          for (size_t c = b; c < e; ++c) {
            for (size_t w : T.cell(c).walls) {
              const Wall &wall = T.wall(w);
              if (Tissue::isBackground(wall.cell1) ||
                  Tissue::isBackground(wall.cell2))
                continue;
              const bool first = wall.cell1 == c;
              const size_t mine = first ? awI : awI + 1;
              const size_t other = first ? awI + 1 : awI;
              const size_t pin = first ? pwI : pwI + 1;
              const double fac =
                  (parameter(2) + parameter(3) * wallData[w][pin]) *
                      cellData[c][aI] -
                  (parameter(0) + parameter(1) * cellData[c][auxI]) *
                      wallData[w][mine];
              outWall[w][mine] += fac;
              outCell[c][aI] -= fac;
              const double diff = parameter(4) * wallData[w][mine];
              outWall[w][mine] -= diff;
              outWall[w][other] += diff;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(ActiveTransportWall, "ActiveTransportWall")

// AUX/LAX-driven influx: each membrane's carrier pulls auxin *in* from the
// neighbour, and the two faces compete for the same shared pool through the
// normalisation k0 + AUX_1 + AUX_2.
//
// The optional flux record differs from DiffusionActiveTransportCell's in
// being accumulated (`+=`) rather than assigned, so it integrates over every
// derivative evaluation the solver makes rather than reporting the current
// flux. That is legacy behaviour and is kept; read it as a running total whose
// scale depends on the solver, not as a rate.
class InfluxActiveTransportCell : public Reaction {
public:
  InfluxActiveTransportCell(const ParameterList &p, const IndexLevels &i) {
    const bool ok = (i.size() == 2 || i.size() == 3) && i[0].size() == 1 &&
                    i[1].size() == 1 && (i.size() == 2 || i[2].size() == 1);
    if (!ok)
      throw std::runtime_error(
          "InfluxActiveTransportCell: one cell variable (auxin) and one paired "
          "wall variable (AUX/LAX); an optional third level names a paired "
          "wall variable to accumulate the flux in.");
    std::vector<size_t> shape{1, 1};
    if (i.size() == 3)
      shape.push_back(1);
    configure("InfluxActiveTransportCell", p, i, 2, shape, {"T", "k0"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t aI = variableIndex(0, 0);
    const size_t awI = variableIndex(1, 0);
    const bool saveFlux = numVariableIndexLevel() == 3;
    const size_t fI = saveFlux ? variableIndex(2, 0) : 0;
    const double rate = parameter(0), k0 = parameter(1);
    parallelScatter1(
        T.numWall(), cellDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            if (Tissue::isBackground(wall.cell1) ||
                Tissue::isBackground(wall.cell2))
              continue;
            const size_t c1 = wall.cell1, c2 = wall.cell2;
            const double facnorm = k0 + wallData[w][awI] + wallData[w][awI + 1];
            const double fac =
                rate * (wallData[w][awI] * cellData[c2][aI] / facnorm -
                        wallData[w][awI + 1] * cellData[c1][aI] / facnorm);
            out[c1][aI] += fac;
            out[c2][aI] -= fac;
            if (saveFlux) {
              if (fac >= 0.0)
                wallData[w][fI] += fac;
              else
                wallData[w][fI + 1] += -fac;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(InfluxActiveTransportCell, "InfluxActiveTransportCell")

} // namespace
} // namespace tissue
