//
// CMT::SisterInhibition: the two faces of a wall talk to each other.
//
// Every patterning rule in this model is cell-autonomous -- a wall reads its
// own curvature, its own stress, its own reinforcement -- and measured on
// real leaves that is not enough. After removing a quadratic surface in
// position, which takes out the leaf's maturation gradient, and the cell's
// own shape and size, a cell's lobe sharpness still correlates with its
// neighbours' at 0.40, 0.22 and 0.15 in three images, and it survives at
// graph distance two where no wall is shared. Something coordinates lobes
// across cells over about two cell diameters, and nothing here can produce
// it.
//
// The form is Abley, Sauret-Gueto, Maree and Coen (2016), Formation of
// polarity convergences underlying shoot outgrowths, whose indirect coupling
// model represents the wall explicitly with several compartments per cell
// edge and has neighbours coordinate through it rather than each reading its
// own geometry. Two ingredients, both of which this model lacks:
//
//   across the wall   the two faces of one wall inhibit each other, so they
//                     break symmetry: the face that gets ahead suppresses its
//                     partner. That is what makes one cell's lobe the other's
//                     indentation rather than both trying to bulge.
//   along the wall    lateral diffusion between neighbouring wall elements of
//                     the same cell, which sets how wide a domain gets and so
//                     how far apart lobes sit.
//
//     dp_w/dt = k_on * S(src_w) / (1 + (h * p_sister(w))^n)
//               - k_off * p_w
//               + D * sum over outline neighbours (p_n - p_w)
//
//     S(kappa) = 1 / (1 + exp(kappa / kappa_half))
//
// which is exactly CMT::SignedCurvatureRecruitment's cue, so with h = 0 and
// D = 0 this reaction *is* that one and the arms of the comparison differ in
// one thing. It is bounded in (0, 1), so p is bounded by k_on/k_off as the
// reinforcement variable m is and the same beta_stiff spans the same
// stiffness range; and it is positive everywhere, so nothing has to be
// clamped.
//
// Not clamping is the point. Concave recruits -- measured, cortical
// microtubule density is 1.191 of the cell mean on concave cortex against
// 0.670 on convex -- but that is a factor of 1.8, not a switch. An earlier
// version took max(-kappa, 0), which makes every convex wall produce exactly
// nothing, and that hard zero does the coupling's job before the coupling
// runs: on the starting mesh it already left 90.2% of wall pairs with one
// face holding three times the other, and sweeping h from 0 to 1000 moved
// the relative face difference only 1.856 to 1.959. A cue that overstates
// the asymmetry it is supposed to be the seed of cannot test whether mutual
// inhibition amplifies it.
//
// The exponent is not decoration. With n = 1 the symmetric steady state
// k_off p* = k_on s / (1 + h p*) is a quadratic with a single positive root,
// so the system has exactly one symmetric solution and cannot break symmetry
// however hard the two faces inhibit each other -- measured, at h = 5 the
// face difference came out 0.76 times the uncoupled one, slightly damped
// rather than amplified. Mutual inhibition is a switch only when it is
// cooperative, which is why every such motif in the literature carries a
// Hill exponent. n >= 2 gives the bistable pair of asymmetric states that
// makes one face win.
//
// `p` is an ordinary wall variable, so it can drive whatever the model
// already reads: WallMechanics::SpringModulated takes an arbitrary wall
// variable as its reinforcement index, and the growth rules take one as their
// inhibition index.
//
// This needs a per-face mesh (make_leaf_init.py --per-face). On a shared-wall
// mesh each wall is its own sister and the inhibition term degenerates to
// simple self-limitation, which the code detects and reports rather than
// silently doing nothing interesting.
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class CMTSisterInhibition : public Reaction {
public:
  CMTSisterInhibition(const ParameterList &p, const IndexLevels &i) {
    configure("CMT::SisterInhibition", p, i, 6, {2},
              {"k_on", "k_off", "kappa_half", "h_inhibit", "D_lateral",
               "n_inhibit"});
    if (parameter(5) < 1.0)
      throw std::runtime_error("CMT::SisterInhibition: n_inhibit must be at "
                               "least 1 (and at least 2 to break symmetry).");
    if (parameter(2) <= 0.0)
      throw std::runtime_error("CMT::SisterInhibition: kappa_half must be "
                               "positive; it sets the cue scale.");
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
                Matrix &) override {
    // vertex -> its sister, from the pairs the init file declared
    std::vector<std::size_t> sis(T.numVertex(), kBackground);
    for (std::size_t i = 0; i < T.numSisterVertex(); ++i) {
      sis[T.sisterVertex(i, 0)] = T.sisterVertex(i, 1);
      sis[T.sisterVertex(i, 1)] = T.sisterVertex(i, 0);
    }

    // wall -> the wall whose endpoints are this one's sisters
    sister_.assign(T.numWall(), kBackground);
    std::size_t paired = 0;
    for (std::size_t w = 0; w < T.numWall(); ++w) {
      const std::size_t a = sis[T.wall(w).vertex1], b = sis[T.wall(w).vertex2];
      if (a == kBackground || b == kBackground)
        continue;
      for (std::size_t x : T.vertex(a).walls)
        if (x != w && T.wall(x).hasVertex(b)) {
          sister_[w] = x;
          ++paired;
          break;
        }
    }

    // walls sharing a vertex within the same cell, for the lateral term
    lateral_.assign(T.numWall(), {});
    for (std::size_t c = 0; c < T.numCell(); ++c) {
      const auto &ws = T.cell(c).walls;
      for (std::size_t k = 0; k < ws.size(); ++k)
        for (std::size_t j = k + 1; j < ws.size(); ++j) {
          const Wall &A = T.wall(ws[k]);
          const Wall &B = T.wall(ws[j]);
          if (A.hasVertex(B.vertex1) || A.hasVertex(B.vertex2)) {
            lateral_[ws[k]].push_back(ws[j]);
            lateral_[ws[j]].push_back(ws[k]);
          }
        }
    }
    std::cerr << "CMT::SisterInhibition: " << paired << " of " << T.numWall()
              << " walls paired with a sister" << std::endl;
    if (paired == 0)
      std::cerr << "  no sisters found -- this needs a per-face mesh, and "
                   "without one the cross-wall term does nothing" << std::endl;
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const std::size_t pIdx = variableIndex(0, 0);
    const std::size_t srcIdx = variableIndex(0, 1);
    const double kOn = parameter(0), kOff = parameter(1);
    const double kHalf = parameter(2);
    const double h = parameter(3), D = parameter(4), nExp = parameter(5);
    if (sister_.size() != T.numWall())
      return;

    parallelFor(T.numWall(), [&](std::size_t begin, std::size_t end) {
      for (std::size_t w = begin; w < end; ++w) {
        const double p = wallData[w][pIdx];
        // Signed curvature, positive where the wall bulges out of its
        // own cell, through the same sigmoid the cell-autonomous version
        // uses: 1 at a deep neck, 1/2 on a straight wall, towards 0 on a
        // lobe. The two faces of one wall read opposite signs, which is the
        // asymmetry the coupling acts on, and it is graded rather than
        // switched so that the coupling is what decides a marginal wall.
        const double src =
            1.0 / (1.0 + std::exp(wallData[w][srcIdx] / kHalf));
        const std::size_t s = sister_[w];
        // A wall with no sister sits on the tissue boundary. Inhibiting it by
        // nothing would let it run away relative to every interior wall, so
        // it is inhibited by itself instead: the same steady state, no edge
        // artefact.
        const double pOther = (s == kBackground) ? p : wallData[s][pIdx];
        const double x = h * (pOther > 0.0 ? pOther : 0.0);
        double d = kOn * src / (1.0 + std::pow(x, nExp)) - kOff * p;
        for (std::size_t n : lateral_[w])
          d += D * (wallData[n][pIdx] - p);
        wallDerivs[w][pIdx] += d;
      }
    });
  }

private:
  std::vector<std::size_t> sister_;
  std::vector<std::vector<std::size_t>> lateral_;
};


// Diagnostic::WallCurvatureSigned
//
// The same measurement as Diagnostic::WallCurvature but keeping the sign:
// positive where the wall bulges out of the cell it belongs to, negative
// where it is indented. Pure instrumentation, run between solver steps.
//
// It exists because the unsigned version cannot express the one asymmetry the
// two faces of a wall actually have. The same physical wall is convex for one
// cell and concave for the other, and this project has measured what that
// does -- cortical microtubule density is 0.670 of the cell mean on convex
// cortex and 1.191 on concave, a factor of 1.8. An unsigned curvature
// discards exactly that, which is why the existing recruitment had to be a
// function of |kappa|: a consequence of the shared-wall mesh rather than a
// modelling choice. On a per-face mesh the two faces get opposite signs here,
// and that is the seed a mutual-inhibition coupling needs -- without some
// asymmetry it sits in its symmetric steady state forever.
class DiagnosticWallCurvatureSigned : public Reaction {
public:
  DiagnosticWallCurvatureSigned(const ParameterList &p, const IndexLevels &i) {
    configure("Diagnostic::WallCurvatureSigned", p, i, 0, {1}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              double) override {
    const std::size_t out = variableIndex(0, 0);
    std::vector<double> k(T.numVertex(), 0.0);
    std::vector<char> ok(T.numVertex(), 0);

    for (std::size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const std::size_t n = cell.numVertex();
      if (n < 3)
        continue;
      // orientation of this cell's outline, so "out of the cell" is well
      // defined whichever way its vertices happen to be listed
      double two = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        const auto a = vertexData[cell.vertices[j]];
        const auto b = vertexData[cell.vertices[(j + 1) % n]];
        two += a[0] * b[1] - b[0] * a[1];
      }
      const double orient = two >= 0.0 ? 1.0 : -1.0;
      for (std::size_t j = 0; j < n; ++j) {
        const std::size_t v = cell.vertices[j];
        const auto p0 = vertexData[cell.vertices[(j + n - 1) % n]];
        const auto p1 = vertexData[v];
        const auto p2 = vertexData[cell.vertices[(j + 1) % n]];
        const double ax = p1[0] - p0[0], ay = p1[1] - p0[1];
        const double bx = p2[0] - p1[0], by = p2[1] - p1[1];
        const double la = std::hypot(ax, ay), lb = std::hypot(bx, by);
        if (la <= 0.0 || lb <= 0.0)
          continue;
        // signed turning angle per unit arc: positive bulging out of the cell
        const double ang = std::atan2(ax * by - ay * bx, ax * bx + ay * by);
        k[v] = orient * ang / (0.5 * (la + lb));
        ok[v] = 1;
      }
    }
    for (std::size_t w = 0; w < T.numWall(); ++w) {
      double sum = 0.0;
      std::size_t n = 0;
      for (std::size_t v : {T.wall(w).vertex1, T.wall(w).vertex2})
        if (v < ok.size() && ok[v]) {
          sum += k[v];
          ++n;
        }
      wallData[w][out] = n ? sum / double(n) : 0.0;
    }
  }
};

} // namespace

TISSUE_REGISTER_REACTION(CMTSisterInhibition, "CMT::SisterInhibition")
TISSUE_REGISTER_REACTION(DiagnosticWallCurvatureSigned,
                         "Diagnostic::WallCurvatureSigned")

} // namespace tissue
