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
//     dp_w/dt = k_on * src_w / (1 + h * p_sister(w))
//               - k_off * p_w
//               + D * sum over outline neighbours (p_n - p_w)
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
    configure("CMT::SisterInhibition", p, i, 4, {2},
              {"k_on", "k_off", "h_inhibit", "D_lateral"});
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
    const double h = parameter(2), D = parameter(3);
    if (sister_.size() != T.numWall())
      return;

    parallelFor(T.numWall(), [&](std::size_t begin, std::size_t end) {
      for (std::size_t w = begin; w < end; ++w) {
        const double p = wallData[w][pIdx];
        const double src = wallData[w][srcIdx];
        const std::size_t s = sister_[w];
        // A wall with no sister sits on the tissue boundary. Inhibiting it by
        // nothing would let it run away relative to every interior wall, so
        // it is inhibited by itself instead: the same steady state, no edge
        // artefact.
        const double pOther = (s == kBackground) ? p : wallData[s][pIdx];
        double d = kOn * src / (1.0 + h * pOther) - kOff * p;
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
