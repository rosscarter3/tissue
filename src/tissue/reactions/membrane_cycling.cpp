//
// PIN (and other carrier) cycling between a cell's interior pool and its
// individual membranes. Ported from legacy membraneCycling.cc and
// membraneCyclingAll.cc.
//
// Every reaction here has the same shape: walk each cell's walls, and for the
// membrane facing that cell move carrier on at rate k_on and off at rate
// k_off, each modulated by something - the opposite membrane, a wall signal,
// the neighbour's auxin, the cell's own auxin, the membrane's own load, or the
// auxin flux it carries. Wall variables are paired, index k being the face of
// wall.cell1 and k+1 that of wall.cell2.
//
// The MembraneCyclingAll:: variants are the same rules with the "the wall must
// separate two real cells" test dropped, so cycling also happens on the
// tissue boundary.
//
// Four legacy typos are fixed here; see "Deliberate fixes over legacy" in
// README.md. Each one breaks a symmetry the rule cannot sensibly break: which
// of a wall's two cells is stored as cell1 is an arbitrary artefact of the
// init file, so the two branches must be mirror images, and in each of these
// four they are not.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// One cell's view of one of its walls.
struct Membrane {
  size_t cell;
  size_t wall;
  size_t own;      // 0 if this cell is the wall's cell1, else 1
  size_t neighbor; // kBackground on a boundary wall
};

// Visits cells [b, e) and their walls in legacy's order.
template <class Fn>
void forEachMembrane(Tissue &T, size_t b, size_t e, bool interiorOnly, Fn fn) {
  for (size_t c = b; c < e; ++c)
    for (size_t w : T.cell(c).walls) {
      const Wall &wall = T.wall(w);
      const size_t other = wall.otherCell(c);
      if (interiorOnly && Tissue::isBackground(other))
        continue;
      fn(Membrane{c, w, wall.cell1 == c ? size_t{0} : size_t{1}, other});
    }
}

inline double hill(double x, double k, double n) {
  const double xn = std::pow(x, n);
  return xn / (xn + std::pow(k, n));
}

// Shared skeleton: `rate` returns the net on-rate for one membrane, which is
// added to that membrane and taken from the cell pool.
template <class Rate>
void cycle(Tissue &T, Matrix &cellDerivs, Matrix &wallDerivs, size_t pI,
           size_t pwI, bool interiorOnly, Rate rate) {
  parallelScatter2(T.numCell(), cellDerivs, wallDerivs,
                   [&](size_t b, size_t e, Matrix &oc, Matrix &ow) {
                     forEachMembrane(T, b, e, interiorOnly, [&](Membrane m) {
                       const double fac = rate(m);
                       ow[m.wall][pwI + m.own] += fac;
                       oc[m.cell][pI] -= fac;
                     });
                   });
}

// --- base rates: no modulation ---------------------------------------------

// dP_mem/dt = k_on P_cell - k_off P_mem.
template <bool kInteriorOnly>
class ConstantCycling : public Reaction {
public:
  ConstantCycling(const ParameterList &p, const IndexLevels &i)
      : Reaction() {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("MembraneCycling::Constant: one cell variable "
                               "(the carrier pool) and one paired wall "
                               "variable (the membrane carrier).");
    configure(kInteriorOnly ? "MembraneCycling::Constant"
                            : "MembraneCyclingAll::Constant",
              p, i, 2, {1, 1}, {"k_on", "k_off"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0), pwI = variableIndex(1, 0);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, kInteriorOnly, [&](Membrane m) {
      return parameter(0) * cellData[m.cell][pI] -
             parameter(1) * wallData[m.wall][pwI + m.own];
    });
  }
};
using ConstantCyclingInterior = ConstantCycling<true>;
using ConstantCyclingAll = ConstantCycling<false>;
TISSUE_REGISTER_REACTION(ConstantCyclingInterior, "MembraneCycling::Constant")
TISSUE_REGISTER_REACTION(ConstantCyclingAll, "MembraneCyclingAll::Constant")

// --- modulated by the *opposite* membrane -----------------------------------

// Cycling on each face is gated by how loaded the face across the wall is, as
// a Hill function - the two membranes of one wall talking to each other.
class CrossMembraneNonLinear : public Reaction {
public:
  CrossMembraneNonLinear(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("MembraneCycling::CrossMembraneNonLinear: one "
                               "cell variable and one paired wall variable.");
    configure("MembraneCycling::CrossMembraneNonLinear", p, i, 4, {1, 1},
              {"k_on", "k_off", "k", "n"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0), pwI = variableIndex(1, 0);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, true, [&](Membrane m) {
      const double h = hill(wallData[m.wall][pwI + 1 - m.own], parameter(2),
                            parameter(3));
      return (parameter(0) * cellData[m.cell][pI] -
              parameter(1) * wallData[m.wall][pwI + m.own]) *
             h;
    });
  }
};
TISSUE_REGISTER_REACTION(CrossMembraneNonLinear,
                         "MembraneCycling::CrossMembraneNonLinear")

// --- modulated by a wall signal ---------------------------------------------

// Cycling gated by a second paired wall variable - a local signal in the wall
// recruiting carrier to the membrane beside it.
template <bool kInteriorOnly>
class LocalWallFeedbackNonLinear : public Reaction {
public:
  LocalWallFeedbackNonLinear(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 2)
      throw std::runtime_error(
          "MembraneCycling::LocalWallFeedbackNonLinear: one cell variable, "
          "then two paired wall variables (the signal, then the carrier).");
    configure(kInteriorOnly ? "MembraneCycling::LocalWallFeedbackNonLinear"
                            : "MembraneCyclingAll::LocalWallFeedbackNonLinear",
              p, i, 4, {1, 2}, {"k_on", "k_off", "k", "n"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0);
    const size_t xwI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, kInteriorOnly, [&](Membrane m) {
      const double h =
          hill(wallData[m.wall][xwI + m.own], parameter(2), parameter(3));
      return (parameter(0) * cellData[m.cell][pI] -
              parameter(1) * wallData[m.wall][pwI + m.own]) *
             h;
    });
  }
};
using LocalWallFeedbackNonLinearInterior = LocalWallFeedbackNonLinear<true>;
using LocalWallFeedbackNonLinearAll = LocalWallFeedbackNonLinear<false>;
TISSUE_REGISTER_REACTION(LocalWallFeedbackNonLinearInterior,
                         "MembraneCycling::LocalWallFeedbackNonLinear")
TISSUE_REGISTER_REACTION(LocalWallFeedbackNonLinearAll,
                         "MembraneCyclingAll::LocalWallFeedbackNonLinear")

// LocalWallFeedbackNonLinear with the wall signal *inhibiting* cycling rather
// than promoting it. The two terms are not simply the same factor: the on-rate
// is divided by (x^n + k^n) while the off-rate is divided by (1 + x^n/k^n),
// so they differ by a factor k^n. That is legacy's algebra and is preserved -
// it rescales k_on relative to k_off but leaves both monotonically inhibited.
//
// FIX: legacy's second branch reads the signal from the *first* face
// (`xwI`) in its off-rate term while using `xwI+1` in the on-rate term, so a
// wall's two membranes saw different signals depending only on which cell the
// init file listed first. Both terms read the membrane's own face here.
class LocalWallFeedbackNonLinearInhibition : public Reaction {
public:
  LocalWallFeedbackNonLinearInhibition(const ParameterList &p,
                                       const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 2)
      throw std::runtime_error(
          "MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition: one cell "
          "variable, then two paired wall variables (signal, carrier).");
    configure("MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition", p, i,
              4, {1, 2}, {"k_on", "k_off", "k", "n"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0);
    const size_t xwI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, false, [&](Membrane m) {
      const double xn = std::pow(wallData[m.wall][xwI + m.own], parameter(3));
      const double kn = std::pow(parameter(2), parameter(3));
      return parameter(0) * cellData[m.cell][pI] / (xn + kn) -
             parameter(1) * wallData[m.wall][pwI + m.own] / (1.0 + xn / kn);
    });
  }
};
TISSUE_REGISTER_REACTION(
    LocalWallFeedbackNonLinearInhibition,
    "MembraneCyclingAll::LocalWallFeedbackNonLinearInhibition")

// Linear form of the same rule.
//
// FIX: legacy's first branch reads `wallData[i][pwI]`, indexing the wall table
// with the *cell* index - so the off-rate was taken from an unrelated wall, or
// out of range wherever a tissue has more cells than walls. Reads the
// membrane's own value here, as the second branch already did.
class LocalWallFeedbackLinear : public Reaction {
public:
  LocalWallFeedbackLinear(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 2)
      throw std::runtime_error(
          "MembraneCycling::LocalWallFeedbackLinear: one cell variable, then "
          "two paired wall variables (the signal, then the carrier).");
    configure("MembraneCycling::LocalWallFeedbackLinear", p, i, 2, {1, 2},
              {"k_on", "k_off"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0);
    const size_t xwI = variableIndex(1, 0), pwI = variableIndex(1, 1);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, true, [&](Membrane m) {
      return (parameter(0) * cellData[m.cell][pI] -
              parameter(1) * wallData[m.wall][pwI + m.own]) *
             wallData[m.wall][xwI + m.own];
    });
  }
};
TISSUE_REGISTER_REACTION(LocalWallFeedbackLinear,
                         "MembraneCycling::LocalWallFeedbackLinear")

// --- modulated by the membrane's own load -----------------------------------

// Self-reinforcing recruitment: a membrane already carrying PIN recruits more.
//
// FIX: legacy's first branch raises the membrane load to `parameter(2)` - the
// half-max constant - instead of `parameter(3)`, the Hill exponent, while the
// second branch uses the exponent correctly. One face of every wall was
// running a different power law from the other.
class PINFeedbackNonLinear : public Reaction {
public:
  PINFeedbackNonLinear(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("MembraneCycling::PINFeedbackNonLinear: one "
                               "cell variable and one paired wall variable.");
    configure("MembraneCycling::PINFeedbackNonLinear", p, i, 4, {1, 1},
              {"k_on", "k_off", "k", "n"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0), pwI = variableIndex(1, 0);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, true, [&](Membrane m) {
      const double h =
          hill(wallData[m.wall][pwI + m.own], parameter(2), parameter(3));
      return (parameter(0) * cellData[m.cell][pI] -
              parameter(1) * wallData[m.wall][pwI + m.own]) *
             h;
    });
  }
};
TISSUE_REGISTER_REACTION(PINFeedbackNonLinear,
                         "MembraneCycling::PINFeedbackNonLinear")

// Linear self-reinforcement.
//
// FIX: legacy's first branch reads `wallData[i][pwI]` in its off-rate term,
// indexing the wall table with the cell index - the same typo as
// LocalWallFeedbackLinear. The second branch squares the membrane's own value,
// which is what both branches do here.
class PINFeedbackLinear : public Reaction {
public:
  PINFeedbackLinear(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error("MembraneCycling::PINFeedbackLinear: one cell "
                               "variable and one paired wall variable.");
    configure("MembraneCycling::PINFeedbackLinear", p, i, 2, {1, 1},
              {"k_on", "k_off"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t pI = variableIndex(0, 0), pwI = variableIndex(1, 0);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, true, [&](Membrane m) {
      const double pw = wallData[m.wall][pwI + m.own];
      return (parameter(0) * cellData[m.cell][pI] - parameter(1) * pw) * pw;
    });
  }
};
TISSUE_REGISTER_REACTION(PINFeedbackLinear,
                         "MembraneCycling::PINFeedbackLinear")

// --- modulated by auxin -----------------------------------------------------

// "Up the gradient": recruitment to a membrane scales with the auxin of the
// cell on the far side of it, so carrier accumulates facing richer neighbours.
// kOwnAuxin instead reads the cell's own auxin (InternalCell), which has no
// directional information and cycles all of a cell's membranes together.
template <bool kNonLinear, bool kOwnAuxin>
class AuxinModulatedCycling : public Reaction {
public:
  AuxinModulatedCycling(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
      throw std::runtime_error(
          "MembraneCycling: two cell variables (auxin, carrier pool) and one "
          "paired wall variable (the membrane carrier).");
    configure(name(), p, i, kNonLinear ? 4 : 2, {2, 1},
              kNonLinear ? std::vector<std::string>{"k_on", "k_off", "k", "n"}
                         : std::vector<std::string>{"k_on", "k_off"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 0);
    cycle(T, cellDerivs, wallDerivs, pI, pwI, true, [&](Membrane m) {
      const double auxin = cellData[kOwnAuxin ? m.cell : m.neighbor][aI];
      const double mod =
          kNonLinear ? hill(auxin, parameter(2), parameter(3)) : auxin;
      return (parameter(0) * cellData[m.cell][pI] -
              parameter(1) * wallData[m.wall][pwI + m.own]) *
             mod;
    });
  }

private:
  static const char *name() {
    if (kOwnAuxin)
      return kNonLinear ? "MembraneCycling::InternalCellNonLinear"
                        : "MembraneCycling::InternalCellLinear";
    return kNonLinear ? "MembraneCycling::CellUpTheGradientNonLinear"
                      : "MembraneCycling::CellUpTheGradientLinear";
  }
};
using CellUpTheGradientNonLinear = AuxinModulatedCycling<true, false>;
using CellUpTheGradientLinear = AuxinModulatedCycling<false, false>;
using InternalCellNonLinear = AuxinModulatedCycling<true, true>;
using InternalCellLinear = AuxinModulatedCycling<false, true>;
TISSUE_REGISTER_REACTION(CellUpTheGradientNonLinear,
                         "MembraneCycling::CellUpTheGradientNonLinear")
TISSUE_REGISTER_REACTION(CellUpTheGradientLinear,
                         "MembraneCycling::CellUpTheGradientLinear")
TISSUE_REGISTER_REACTION(InternalCellNonLinear,
                         "MembraneCycling::InternalCellNonLinear")
TISSUE_REGISTER_REACTION(InternalCellLinear,
                         "MembraneCycling::InternalCellLinear")

// Flux-based feedback: auxin is effluxed across each membrane (passively at
// p_0, PIN-driven at p_1), and carrier is *removed* from a membrane in
// proportion to the square of the net outward flux it carries, but only while
// that flux is outward. An inward-facing membrane is left alone rather than
// gaining carrier, which is legacy's asymmetry and is kept.
class CellFluxExocytosis : public Reaction {
public:
  CellFluxExocytosis(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
      throw std::runtime_error("MembraneCycling::CellFluxExocytosis: two cell "
                               "variables (auxin, carrier pool) and one paired "
                               "wall variable (membrane carrier).");
    configure("MembraneCycling::CellFluxExocytosis", p, i, 3, {2, 1},
              {"D", "T", "k_flux"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData, Matrix &,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &) override {
    const size_t aI = variableIndex(0, 0), pI = variableIndex(0, 1);
    const size_t pwI = variableIndex(1, 0);
    parallelScatter2(
        T.numCell(), cellDerivs, wallDerivs,
        [&](size_t b, size_t e, Matrix &oc, Matrix &ow) {
          forEachMembrane(T, b, e, true, [&](Membrane m) {
            const double out =
                parameter(0) * cellData[m.cell][aI] +
                parameter(1) * cellData[m.cell][aI] * wallData[m.wall][pwI + m.own];
            oc[m.cell][aI] -= out;
            oc[m.neighbor][aI] += out;
            const double back =
                parameter(0) * cellData[m.neighbor][aI] +
                parameter(1) * cellData[m.neighbor][aI] *
                    wallData[m.wall][pwI + 1 - m.own];
            const double net = out - back;
            if (net < 0.0)
              return; // inward-facing membranes are not cycled at all
            const double fac = -parameter(2) * cellData[m.cell][pI] * net * net;
            ow[m.wall][pwI + m.own] -= fac;
            oc[m.cell][pI] += fac;
          });
        });
  }
};
TISSUE_REGISTER_REACTION(CellFluxExocytosis,
                         "MembraneCycling::CellFluxExocytosis")

} // namespace
} // namespace tissue
