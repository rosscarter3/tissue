//
// Pavement-cell morphogenesis: the cortical-microtubule / cellulose feedback
// that turns a smooth polygonal epidermal cell into an interdigitated
// "jigsaw puzzle" cell.
//
// The mechanism implemented here is the synthesis of three converging lines
// of evidence (see ../../../../lit/README.md):
//
//  - Sampathkumar et al. (2014) eLife 3:e01967: cell shape sets the stress
//    pattern, tensile stress concentrates in the indenting neck regions,
//    cortical microtubules (CMTs) localise there, and AFM finds stiffer wall
//    exactly there.
//  - Majda et al. (2017) Dev Cell and Bidhendi et al. (2019) Plant Physiol:
//    CMTs, aligned cellulose and de-esterified pectin all co-localise at the
//    necks; lobe amplification is driven by the stiffening they produce.
//  - Sapala et al. (2018) eLife 7:e32794: lobes arise from *restriction of
//    growth in the indentations*, not promotion of growth in the
//    protrusions, and need no mobile signal to coordinate a lobe with the
//    facing indentation.
//
// The last point is what makes a shared-wall vertex model the natural
// discretisation: one wall, two cells, so a lobe of one cell IS an
// indentation of its neighbour by construction.
//
// Why local curvature is the recruitment cue. A microtubule is stiff on the
// scale of a cell (persistence length ~mm). Against a concave cortex it lies
// comfortably along the surface; against a convex one it must bend against
// its own rigidity and detaches or undergoes catastrophe. That is precisely
// the physics the tubulaton microtubule simulator implements at the membrane
// (its Angle_mb_limite forces a catastrophe when the required turn is too
// large), so the curvature-to-CMT-density response used here is *measurable*
// rather than assumed -- see ../../../../model/README.md for the calibration
// against tubulaton.
//
// A wall is shared, and a neck for one cell is a lobe for the other, so
// exactly one of the two adjacent cells sees any given wall as concave. The
// reinforcement that cell deposits stiffens the shared wall. The recruitment
// cue is therefore the *unsigned* local curvature, which makes the mechanism
// automatically symmetric: which cell ends up lobing into which is broken by
// noise and by cell size, not prescribed.
//
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Discrete unsigned curvature of the wall chain at every vertex.
//
// At a vertex v joined by exactly two walls (a chain-interior vertex; a
// tri-cellular junction has three and is skipped) with chain neighbours a and
// b, the sagitta vector s = (x_a + x_b)/2 - x_v has |s| ~ kappa h^2 / 2 for a
// smooth curve of curvature kappa sampled at spacing h, so
//
//     kappa_v = 2 |s| / h^2,     h = (|x_a - x_v| + |x_b - x_v|) / 2
//
// which is mesh-spacing independent to leading order -- important here, since
// wall segments lengthen as the tissue grows. Junction vertices get kappa = 0
// and are excluded from the per-wall average rather than counted as flat.
void computeVertexCurvature(const Tissue &T, const Matrix &vertexData,
                            std::vector<double> &kappa,
                            std::vector<char> &valid) {
  const size_t dim = vertexData.cols();
  kappa.assign(T.numVertex(), 0.0);
  valid.assign(T.numVertex(), 0);
  parallelFor(T.numVertex(), [&](size_t begin, size_t end) {
    for (size_t v = begin; v < end; ++v) {
      const VertexTopo &vt = T.vertex(v);
      if (vt.walls.size() != 2)
        continue; // junction or chain end: curvature undefined
      const size_t a = T.wall(vt.walls[0]).otherVertex(v);
      const size_t b = T.wall(vt.walls[1]).otherVertex(v);
      double sag = 0.0, da = 0.0, db = 0.0;
      for (size_t d = 0; d < dim; ++d) {
        const double xa = vertexData[a][d];
        const double xb = vertexData[b][d];
        const double xv = vertexData[v][d];
        const double s = 0.5 * (xa + xb) - xv;
        sag += s * s;
        da += (xa - xv) * (xa - xv);
        db += (xb - xv) * (xb - xv);
      }
      const double h = 0.5 * (std::sqrt(da) + std::sqrt(db));
      if (h <= 0.0)
        continue;
      kappa[v] = 2.0 * std::sqrt(sag) / (h * h);
      valid[v] = 1;
    }
  });
}

// Mean of the valid endpoint curvatures of wall w (0 if neither endpoint is a
// chain-interior vertex).
inline double wallCurvature(const Tissue &T, size_t w,
                            const std::vector<double> &kappa,
                            const std::vector<char> &valid) {
  const Wall &wall = T.wall(w);
  double sum = 0.0;
  size_t n = 0;
  if (valid[wall.vertex1]) {
    sum += kappa[wall.vertex1];
    ++n;
  }
  if (valid[wall.vertex2]) {
    sum += kappa[wall.vertex2];
    ++n;
  }
  return n ? sum / static_cast<double>(n) : 0.0;
}

// ---------------------------------------------------------------------------
// CMT::CurvatureRecruitment
//
//   dm/dt = k_on * S * (1 - m) - k_off * m,   m in [0, 1]
//
// with m the per-wall cortical-microtubule / cellulose reinforcement density
// and S the recruitment cue, a Hill function of the local unsigned curvature:
//
//   S = kappa^n / (kappa_half^n + kappa^n)
//
// In the six-parameter form the cue is blended with the wall's elastic
// strain, which is the Sampathkumar et al. (2014) stress-recruitment term:
//
//   S = (1 - w_strain) * Hill(kappa) + w_strain * Hill(eps),  eps = (d-L)/L
//
// An optional cue_mode inverts the curvature term to S = 1 - Hill(kappa),
// i.e. reinforcement is highest on FLAT wall and falls as curvature grows.
// That is not an alternative guess, it is what the filament physics actually
// gives. Running tubulaton inside these simulated cell geometries (see
// ../../../../model/README.md; 144 495 pooled cortical tubulin elements over
// five runs in four cells) measures a single cell's cortical density as
//
//   convex (lobe tip)  0.670 +/- 0.109     (relative to that cell's mean)
//   flat               1.429 +/- 0.183
//   concave (neck)     1.191 +/- 0.220
//
// so a stiff filament is depleted ~2x on convex cortex, while concave and
// flat are indistinguishable at this precision. A wall is shared and a neck
// of one cell is a lobe of the other, so its total reinforcement is
// g(kappa) + g(-kappa): 1.86 for curved wall against 2.86 for flat, i.e.
// flat wall carries ~1.5x more. Hence the inverted, |kappa|-decreasing cue.
//
// The saturating (1 - m) factor bounds m in [0, 1] for any m(0) in [0, 1], so
// downstream reactions can treat m as a fraction of maximal reinforcement.
class CmtCurvatureRecruitment : public Reaction {
public:
  CmtCurvatureRecruitment(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 4 && p.size() != 5 && p.size() != 6 && p.size() != 7)
      throw std::runtime_error(
          "CMT::CurvatureRecruitment: uses four to seven parameters "
          "(k_on, k_off, kappa_half, n_hill, [w_strain, strain_half], "
          "[cue_mode: 0 = curvature-enriched (default), 1 = "
          "curvature-depleted]).");
    if (p.size() == 6 && (p[4] < 0.0 || p[4] > 1.0))
      throw std::runtime_error(
          "CMT::CurvatureRecruitment: w_strain (parameter 5) must lie in "
          "[0, 1].");
    if (p[2] <= 0.0)
      throw std::runtime_error(
          "CMT::CurvatureRecruitment: kappa_half (parameter 3) must be "
          "positive.");
    const bool withStrain = p.size() >= 6;
    if (p.size() == 5 || p.size() == 7) {
      const double mode = p.back();
      if (mode != 0.0 && mode != 1.0)
        throw std::runtime_error(
            "CMT::CurvatureRecruitment: cue_mode must be 0 "
            "(curvature-enriched) or 1 (curvature-depleted).");
    }
    if (i.size() != (withStrain ? 2u : 1u) || i[0].size() != 1 ||
        (withStrain && i[1].size() != 1))
      throw std::runtime_error(
          "CMT::CurvatureRecruitment: level 0 = wall CMT variable index; "
          "level 1 = wall resting-length index (six-parameter form only).");
    std::vector<std::string> ids{"k_on", "k_off", "kappa_half", "n_hill"};
    if (withStrain) {
      ids.push_back("w_strain");
      ids.push_back("strain_half");
    }
    if (p.size() == 5 || p.size() == 7)
      ids.push_back("cue_mode");
    configure("CMT::CurvatureRecruitment", p, i, p.size(),
              std::vector<size_t>(i.size(), 1), std::move(ids));
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t cmtIndex = variableIndex(0, 0);
    const double kOn = parameter(0);
    const double kOff = parameter(1);
    const double n = parameter(3);
    const double kappaHalfN = std::pow(parameter(2), n);
    const bool withStrain = numParameter() >= 6;
    const double wStrain = withStrain ? parameter(4) : 0.0;
    const double strainHalfN = withStrain ? std::pow(parameter(5), n) : 1.0;
    const size_t lengthIndex = withStrain ? variableIndex(1, 0) : 0;
    const bool depleted =
        (numParameter() == 5 || numParameter() == 7) &&
        parameter(numParameter() - 1) == 1.0;

    computeVertexCurvature(T, vertexData, kappa_, valid_);

    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        const double kappa = wallCurvature(T, w, kappa_, valid_);
        const double kN = std::pow(kappa, n);
        double cue = kN / (kappaHalfN + kN);
        if (depleted)
          cue = 1.0 - cue;
        if (withStrain) {
          const double L = wallData[w][lengthIndex];
          const double d = T.wallLengthFromVertices(w, vertexData);
          const double eps = L > 0.0 ? (d - L) / L : 0.0;
          const double e = eps > 0.0 ? eps : 0.0;
          const double eN = std::pow(e, n);
          cue = (1.0 - wStrain) * cue + wStrain * (eN / (strainHalfN + eN));
        }
        const double m = wallData[w][cmtIndex];
        wallDerivs[w][cmtIndex] += kOn * cue * (1.0 - m) - kOff * m;
      }
    });
  }

private:
  std::vector<double> kappa_;
  std::vector<char> valid_;
};
TISSUE_REGISTER_REACTION(CmtCurvatureRecruitment, "CMT::CurvatureRecruitment")

// ---------------------------------------------------------------------------
// WallMechanics::SpringModulated
//
// WallMechanics::Spring with the spring constant scaled continuously by a
// per-wall reinforcement variable m:
//
//   K_eff = K_force * (1 + beta * m)
//
// so beta is the fold stiffening at full reinforcement, minus one. The legacy
// Spring reaction can only switch between two stiffnesses on a binary wall
// flag, which cannot represent a graded cellulose density. Sampathkumar et
// al. (2014) measure a 5x modulus ratio between cellulose microfibrils and
// the wall matrix, so beta ~ 4 is the literature-scale value.
//
// Everything else follows the legacy spring idiom exactly: the coefficient
// K_eff (1/L_rest - 1/d), zeroed only when both d and L_rest vanish, scaled by
// frac_adh when the wall is stretched, applied as (x_a - x_b) * coeff.
class WallMechanicsSpringModulated : public Reaction {
public:
  WallMechanicsSpringModulated(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3)
      throw std::runtime_error(
          "WallMechanics::SpringModulated: uses three parameters "
          "(K_force, frac_adhesion, beta_stiffening).");
    if (i.size() != 2 && i.size() != 3)
      throw std::runtime_error(
          "WallMechanics::SpringModulated: level 0 = wall resting-length "
          "index, level 1 = wall reinforcement index, optional level 2 = "
          "wall index to save the spring force into.");
    configure("WallMechanics::SpringModulated", p, i, 3,
              std::vector<size_t>(i.size(), 1),
              {"K_force", "frac_adh", "beta_stiffening"});
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const size_t cmtIndex = variableIndex(1, 0);
    const bool saveForce = numVariableIndexLevel() == 3;
    const size_t forceIndex = saveForce ? variableIndex(2, 0) : 0;
    const double kForce = parameter(0);
    const double fracAdh = parameter(1);
    const double beta = parameter(2);
    const size_t dim = vertexData.cols();

    parallelScatter1(
        T.numWall(), vertexDerivs, [&](size_t begin, size_t end, Matrix &out) {
          for (size_t w = begin; w < end; ++w) {
            const size_t v1 = T.wall(w).vertex1;
            const size_t v2 = T.wall(w).vertex2;
            double d = 0.0;
            for (size_t dd = 0; dd < dim; ++dd) {
              const double diff = vertexData[v1][dd] - vertexData[v2][dd];
              d += diff * diff;
            }
            d = std::sqrt(d);
            const double L = wallData[w][lengthIndex];
            const double m = wallData[w][cmtIndex];
            const double kEff = kForce * (1.0 + beta * m);
            double coeff = kEff * (1.0 / L - 1.0 / d);
            if (d <= 0.0 && L <= 0.0)
              coeff = 0.0;
            if (d > L)
              coeff *= fracAdh;
            if (saveForce)
              wallData[w][forceIndex] = coeff * d;
            for (size_t dd = 0; dd < dim; ++dd) {
              const double div =
                  (vertexData[v1][dd] - vertexData[v2][dd]) * coeff;
              out[v1][dd] -= div;
              out[v2][dd] += div;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(WallMechanicsSpringModulated,
                         "WallMechanics::SpringModulated")

// ---------------------------------------------------------------------------
// WallGrowth::StrainWallInhibited
//
// Irreversible wall extension whose rate is damped by the per-wall
// reinforcement variable. This is the "restriction of growth in the
// indentations" that Sapala et al. (2018) identify as the actual driver of
// lobing -- lobes are not pushed out, the necks are held back. gamma is the
// fold reduction in extensibility at full reinforcement, minus one.
//
// Two drives, selected by the optional growth_mode parameter:
//
//  1 (strain-gated, the default) -- Lockhart yielding, as in Bozorg et al.
//    (2016) Phys Biol 13:065002:
//
//      dL/dt = k (eps - eps_th)+ L (1 - L/L_max)+ / (1 + gamma m)
//
//  0 (active) -- wall synthesis proceeding at its own rate, with the
//    threshold acting only as an on/off gate:
//
//      dL/dt = k L (1 - L/L_max)+ / (1 + gamma m),   for eps > eps_th
//
// Mode 0 is what lobing needs, and the reason is worth stating. Under pure
// strain gating the outline's extension is slaved to the cell's areal growth:
// if the perimeter runs ahead, the elastic strain falls, growth slows, and
// the cell simply relaxes to a rounder shape -- measured here as circularity
// rising from 0.892 to 0.915 over 24 h with no lobing whatsoever. Lobes
// require the anticlinal outline to extend faster than sqrt(area), which is
// precisely the periclinal/anticlinal growth-rate mismatch documented by
// Armour et al. (2015) Plant Physiol 167:1039. With an independent synthesis
// rate the outline carries more length than a smooth curve of that area can
// hold; the excess goes slack and, against the wall's bending stiffness,
// buckles at a selected wavelength.
//
// Set eps_th slightly negative in mode 0 so a wall keeps extending as it goes
// slack but stops before being driven deep into compression.
//
// Note the asymmetry with the stiffening term: stiffening alone only changes
// the elastic response and would not by itself select a wavelength, whereas
// gating the irreversible extension is what makes arc length accumulate in
// the low-curvature stretches and amplify the undulation.
class WallGrowthStrainWallInhibited : public Reaction {
public:
  WallGrowthStrainWallInhibited(const ParameterList &p, const IndexLevels &i) {
    if (p.size() < 3 || p.size() > 5)
      throw std::runtime_error(
          "WallGrowth::StrainWallInhibited: uses three to five parameters "
          "(k_growth, strain_threshold, gamma_inhibition, [L_max, 0=off], "
          "[growth_mode: 1=strain-gated (default), 0=active]).");
    if (p[2] < 0.0)
      throw std::runtime_error(
          "WallGrowth::StrainWallInhibited: gamma_inhibition (parameter 3) "
          "must be non-negative.");
    if (p.size() == 5 && p[4] != 0.0 && p[4] != 1.0)
      throw std::runtime_error(
          "WallGrowth::StrainWallInhibited: growth_mode (parameter 5) must "
          "be 0 (active) or 1 (strain-gated).");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "WallGrowth::StrainWallInhibited: level 0 = wall resting-length "
          "index, level 1 = wall reinforcement index.");
    std::vector<std::string> ids{"k_growth", "strain_threshold",
                                 "gamma_inhibition"};
    if (p.size() >= 4)
      ids.push_back("L_max");
    if (p.size() == 5)
      ids.push_back("growth_mode");
    configure("WallGrowth::StrainWallInhibited", p, i, p.size(), {1, 1},
              std::move(ids));
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const size_t cmtIndex = variableIndex(1, 0);
    const double k = parameter(0);
    const double threshold = parameter(1);
    const double gamma = parameter(2);
    const double lMax = numParameter() >= 4 ? parameter(3) : 0.0;
    const bool strainGated = numParameter() < 5 || parameter(4) == 1.0;

    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        const double L = wallData[w][lengthIndex];
        if (L <= 0.0)
          continue;
        const double d = T.wallLengthFromVertices(w, vertexData);
        const double eps = (d - L) / L;
        if (eps <= threshold)
          continue; // below yield, or too deep into compression
        double sat = 1.0;
        if (lMax > 0.0) {
          sat = 1.0 - L / lMax;
          if (sat < 0.0)
            sat = 0.0;
        }
        const double drive = strainGated ? (eps - threshold) : 1.0;
        const double m = wallData[w][cmtIndex];
        wallDerivs[w][lengthIndex] +=
            k * drive * L * sat / (1.0 + gamma * m);
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthStrainWallInhibited,
                         "WallGrowth::StrainWallInhibited")

// ---------------------------------------------------------------------------
// Diagnostic::WallCurvature
//
// Writes the discrete unsigned wall curvature into a wall variable so it can
// be written out and plotted. Pure instrumentation: it contributes nothing to
// any derivative and runs in update(), between solver steps.
class DiagnosticWallCurvature : public Reaction {
public:
  DiagnosticWallCurvature(const ParameterList &p, const IndexLevels &i) {
    configure("Diagnostic::WallCurvature", p, i, 0, {1}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              double) override {
    const size_t outIndex = variableIndex(0, 0);
    computeVertexCurvature(T, vertexData, kappa_, valid_);
    for (size_t w = 0; w < T.numWall(); ++w)
      wallData[w][outIndex] = wallCurvature(T, w, kappa_, valid_);
  }

private:
  std::vector<double> kappa_;
  std::vector<char> valid_;
};
TISSUE_REGISTER_REACTION(DiagnosticWallCurvature, "Diagnostic::WallCurvature")


// ---------------------------------------------------------------------------
// CMT::TransverseReinforcement
//
// The periclinal (outer) wall, which a purely outline-based 2D model
// otherwise omits entirely -- and without which the cell simply inflates
// smoothly, because a wall chain held in tension cannot wrinkle.
//
// Cortical microtubule bundles run ACROSS the periclinal face of a pavement
// cell, anchored at an indentation on one side and at the indentation facing
// it on the other (Sampathkumar et al. 2014; Belteton et al. 2018 report the
// morphologically potent CMTs spanning the periclinal and anticlinal faces).
// The cellulose they template is therefore a band spanning the cell, and it
// resists the cell getting wider -- it pinches the cell in. This is the
// element Sapala et al. (2018) introduce as transverse springs, and it is
// what actually creates indentations; the outline stiffening on its own only
// changes the elastic response.
//
// Implementation, following their rules:
//
//  - Each outline vertex is paired with the nearest vertex facing it across
//    the cell interior (chord within cos_tol of both inward normals).
//  - A pair pulls only when its span exceeds L_target, and only inwards: it
//    is a tether, not a strut. L_target is the cell width the reinforcement
//    holds the cell to, i.e. Sapala's target largest empty circle.
//  - Attachment is inhibited where the wall is MORE convex than a circle of
//    the same perimeter, so indentations accumulate connections and deepen
//    while lobe tips accumulate none. This is the geometric asymmetry that
//    makes the pattern amplify. The comparison has to be against the cell's
//    own mean curvature 2 pi / P rather than against zero: a young pavement
//    cell is convex everywhere, so an absolute test finds no attachment site
//    anywhere and the mechanism can never start. (It was written that way
//    first, and produced results bit-identical to having the reaction
//    switched off entirely.) convex_cut scales that threshold.
//  - Tether strength is scaled by (1 + beta_m * m) with m the local
//    reinforcement from CMT::CurvatureRecruitment, so the curvature-recruited
//    microtubule density -- the quantity tubulaton predicts -- sets how
//    strongly a neck is held in.
//
// Onset is therefore purely geometric (a cell must first grow wider than
// L_target), which reproduces the observation that lobing begins as cells
// enlarge, while amplification is microtubule-dependent, which reproduces the
// loss of lobing when microtubules are depolymerised.
//
// The pairing is an O(n^2) search per cell and is refreshed in update(),
// between solver steps, not inside every derivative evaluation -- cortical
// array reorganisation is an hours-timescale process, and holding the
// connection set fixed within a step keeps the adaptive solver from chasing
// its own feedback between Runge-Kutta stages.
class CmtTransverseReinforcement : public Reaction {
public:
  CmtTransverseReinforcement(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 5 && p.size() != 6)
      throw std::runtime_error(
          "CMT::TransverseReinforcement: uses five or six parameters "
          "(K_tether, L_target, beta_m, cos_tol, refresh_interval, "
          "[convex_cut, default 1]).");
    if (p[1] <= 0.0)
      throw std::runtime_error(
          "CMT::TransverseReinforcement: L_target (parameter 2) must be "
          "positive.");
    if (p[3] <= 0.0 || p[3] > 1.0)
      throw std::runtime_error(
          "CMT::TransverseReinforcement: cos_tol (parameter 4) must lie in "
          "(0, 1].");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "CMT::TransverseReinforcement: level 0 = wall reinforcement "
          "index.");
    std::vector<std::string> ids{"K_tether", "L_target", "beta_m", "cos_tol",
                                 "refresh_interval"};
    if (p.size() == 6)
      ids.push_back("convex_cut");
    configure("CMT::TransverseReinforcement", p, i, p.size(), {1},
              std::move(ids));
  }

  void initiate(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    rebuild(T, wallData, vertexData);
    nextRefresh_ = parameter(4);
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "CMT::TransverseReinforcement requires a 2D tissue.");
    const double lTarget = parameter(1);
    for (const Tether &t : tethers_) {
      const double dx = vertexData[t.v2][0] - vertexData[t.v1][0];
      const double dy = vertexData[t.v2][1] - vertexData[t.v1][1];
      const double d = std::sqrt(dx * dx + dy * dy);
      if (d <= lTarget || d <= 0.0)
        continue; // slack: a tether, not a strut
      const double f = t.strength * (d - lTarget) / d;
      vertexDerivs[t.v1][0] += f * dx;
      vertexDerivs[t.v1][1] += f * dy;
      vertexDerivs[t.v2][0] -= f * dx;
      vertexDerivs[t.v2][1] -= f * dy;
    }
  }

  void update(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              double h) override {
    elapsed_ += h;
    if (elapsed_ < nextRefresh_)
      return;
    elapsed_ = 0.0;
    nextRefresh_ = parameter(4);
    rebuild(T, wallData, vertexData);
  }

private:
  struct Tether {
    size_t v1;
    size_t v2;
    double strength;
  };

  // Mean reinforcement of the walls meeting at a vertex.
  static double vertexReinforcement(const Tissue &T, const Matrix &wallData,
                                    size_t v, size_t mIndex) {
    double sum = 0.0;
    size_t n = 0;
    for (size_t w : T.vertex(v).walls) {
      sum += wallData[w][mIndex];
      ++n;
    }
    return n ? sum / static_cast<double>(n) : 0.0;
  }

  void rebuild(Tissue &T, Matrix &wallData, Matrix &vertexData) {
    const size_t mIndex = variableIndex(0, 0);
    const double kTether = parameter(0);
    const double betaM = parameter(2);
    const double cosTol = parameter(3);
    const double convexCut = numParameter() == 6 ? parameter(5) : 1.0;
    tethers_.clear();
    if (vertexData.cols() != 2)
      return;

    std::vector<double> nx, ny;
    std::vector<char> attachable;
    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numVertex();
      if (n < 6)
        continue;
      // Signed area fixes which way "into the cell" points.
      double signedArea = T.cellVolume(c, vertexData, true);
      const double orient = signedArea < 0.0 ? -1.0 : 1.0;

      // Mean curvature of a closed outline of this perimeter, 2 pi / P: the
      // reference against which a stretch counts as "too convex to attach".
      double outlinePerimeter = 0.0;
      for (size_t k = 0; k < n; ++k) {
        const size_t v = cell.vertices[k];
        const size_t vNext = cell.vertices[(k + 1) % n];
        outlinePerimeter +=
            std::hypot(vertexData[vNext][0] - vertexData[v][0],
                       vertexData[vNext][1] - vertexData[v][1]);
      }
      if (outlinePerimeter <= 0.0)
        continue;
      const double kappaMean = 2.0 * M_PI / outlinePerimeter;

      nx.assign(n, 0.0);
      ny.assign(n, 0.0);
      attachable.assign(n, 0);
      for (size_t k = 0; k < n; ++k) {
        const size_t vPrev = cell.vertices[(k + n - 1) % n];
        const size_t v = cell.vertices[k];
        const size_t vNext = cell.vertices[(k + 1) % n];
        const double tx = vertexData[vNext][0] - vertexData[vPrev][0];
        const double ty = vertexData[vNext][1] - vertexData[vPrev][1];
        const double tl = std::sqrt(tx * tx + ty * ty);
        if (tl <= 0.0)
          continue;
        // Inward normal: rotate the outline tangent by -90 deg for a
        // counter-clockwise cell, +90 deg for a clockwise one.
        nx[k] = orient * ty / tl;
        ny[k] = -orient * tx / tl;
        // Convexity gate. The signed turning angle per unit arc length is
        // the local signed curvature (positive = bulging out of the cell, a
        // lobe tip). A site is attachable when it is less convex than the
        // cell's own mean curvature, scaled by convex_cut.
        const double ax = vertexData[v][0] - vertexData[vPrev][0];
        const double ay = vertexData[v][1] - vertexData[vPrev][1];
        const double bx = vertexData[vNext][0] - vertexData[v][0];
        const double by = vertexData[vNext][1] - vertexData[v][1];
        const double la = std::hypot(ax, ay);
        const double lb = std::hypot(bx, by);
        if (la <= 0.0 || lb <= 0.0)
          continue;
        double sinTurn = orient * (ax * by - ay * bx) / (la * lb);
        sinTurn = std::max(-1.0, std::min(1.0, sinTurn));
        const double kappaSigned = std::asin(sinTurn) / (0.5 * (la + lb));
        attachable[k] = kappaSigned <= convexCut * kappaMean ? 1 : 0;
      }

      // Pair each attachable vertex with the nearest vertex facing it across
      // the interior. Keep i < j so a pair is created once.
      for (size_t k = 0; k < n; ++k) {
        if (!attachable[k])
          continue;
        const size_t v = cell.vertices[k];
        double best = 0.0;
        size_t bestIdx = n;
        for (size_t l = 0; l < n; ++l) {
          if (l == k || !attachable[l])
            continue;
          // Skip near neighbours along the outline: a chord between them
          // spans wall, not cell interior.
          const size_t sep = k < l ? l - k : k - l;
          if (sep < n / 6 || n - sep < n / 6)
            continue;
          const size_t u = cell.vertices[l];
          const double dx = vertexData[u][0] - vertexData[v][0];
          const double dy = vertexData[u][1] - vertexData[v][1];
          const double d = std::sqrt(dx * dx + dy * dy);
          if (d <= 0.0)
            continue;
          // The chord must run inwards at both ends, i.e. the two vertices
          // face each other across the cell.
          if ((dx * nx[k] + dy * ny[k]) / d < cosTol)
            continue;
          if ((-dx * nx[l] - dy * ny[l]) / d < cosTol)
            continue;
          if (bestIdx == n || d < best) {
            best = d;
            bestIdx = l;
          }
        }
        if (bestIdx == n)
          continue;
        const size_t u = cell.vertices[bestIdx];
        if (v >= u)
          continue; // the partner will add this pair from its own side
        const double m = 0.5 * (vertexReinforcement(T, wallData, v, mIndex) +
                                vertexReinforcement(T, wallData, u, mIndex));
        tethers_.push_back({v, u, kTether * (1.0 + betaM * m)});
      }
    }
  }

  std::vector<Tether> tethers_;
  double elapsed_ = 0.0;
  double nextRefresh_ = 0.0;
};
TISSUE_REGISTER_REACTION(CmtTransverseReinforcement,
                         "CMT::TransverseReinforcement")


// ---------------------------------------------------------------------------
// Pressure2D::AreaElastic
//
// Turgor plus an elastic resistance to areal expansion, with its own slowly
// growing resting area:
//
//   P_eff = P_turgor + K_area * (A0 - A) / A0
//   dA0/dt = k_area * A0
//
// applied to the vertices through the same shoelace area gradient that
// Pressure2D::AreaPotential uses, so P_eff > 0 always inflates regardless of
// the tissue's sorting orientation.
//
// This is the periclinal (outer) face of the cell, which a model of the
// anticlinal outline alone leaves out -- and leaving it out is why such a
// model cannot lobe. With a constant turgor and nothing resisting area, a
// cell whose perimeter grows simply inflates to match, staying smooth. The
// periclinal wall is a real sheet with its own extensibility, and A0 is its
// resting area.
//
// Lobing then follows from a growth-rate mismatch that is directly
// documented: Armour, Barton, Overall & Wasteneys (2015) Plant Physiol 167:
// 1039, "Differential growth in periclinal and anticlinal walls during lobe
// formation in Arabidopsis cotyledon pavement cells". When the anticlinal
// outline extends faster than sqrt(A0), the outline carries more length than
// a smooth curve of that area can hold; the excess goes slack, and with wall
// bending stiffness it buckles at a selected wavelength. Growth then pauses
// while the wall is slack and resumes once areal growth has re-tensioned it,
// which reproduces the punctuated alternation between lobe initiation and
// isotropic expansion that Belteton et al. observe in time lapse.
//
// Where the slack goes is set by CMT::CurvatureRecruitment and
// WallGrowth::StrainWallInhibited: curvature recruits reinforcement, the
// reinforced stretches stop extending, and the undulation locks in.
class Pressure2DAreaElastic : public Reaction {
public:
  Pressure2DAreaElastic(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3)
      throw std::runtime_error(
          "Pressure2D::AreaElastic: uses three parameters "
          "(P_turgor, K_area, k_areaGrowth).");
    if (p[1] < 0.0)
      throw std::runtime_error(
          "Pressure2D::AreaElastic: K_area (parameter 2) must be "
          "non-negative.");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "Pressure2D::AreaElastic: level 0 = cell variable holding the "
          "resting area A0.");
    configure("Pressure2D::AreaElastic", p, i, 3, {1},
              {"P_turgor", "K_area", "k_areaGrowth"});
  }

  // Start at rest: a cell whose stored A0 is not positive takes its current
  // area, so an init file need not know the geometry. A positive stored value
  // is respected, which lets a model prescribe pre-stressed cells.
  void initiate(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const size_t a0Index = variableIndex(0, 0);
    for (size_t c = 0; c < T.numCell(); ++c)
      if (cellData[c][a0Index] <= 0.0)
        cellData[c][a0Index] = T.cellVolume(c, vertexData);
  }

  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error("Pressure2D::AreaElastic requires 2D.");
    const size_t a0Index = variableIndex(0, 0);
    const double pTurgor = parameter(0);
    const double kArea = parameter(1);
    const double kGrowth = parameter(2);

    for (size_t c = 0; c < T.numCell(); ++c)
      cellDerivs[c][a0Index] += kGrowth * cellData[c][a0Index];

    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numVertex();
            const double signedArea = T.cellVolume(c, vertexData, true);
            const double area = std::fabs(signedArea);
            const double a0 = cellData[c][a0Index];
            double pEff = pTurgor;
            if (a0 > 0.0)
              pEff += kArea * (a0 - area) / a0;
            double factor = 0.5 * pEff;
            if (signedArea < 0.0)
              factor = -factor; // P > 0 always inflates
            for (size_t k = 0; k < n; ++k) {
              const size_t v = cell.vertices[k];
              const size_t vPlus = cell.vertices[(k + 1) % n];
              const size_t vMinus = cell.vertices[(k + n - 1) % n];
              out[v][0] +=
                  factor * (vertexData[vPlus][1] - vertexData[vMinus][1]);
              out[v][1] +=
                  factor * (vertexData[vMinus][0] - vertexData[vPlus][0]);
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure2DAreaElastic, "Pressure2D::AreaElastic")


// ---------------------------------------------------------------------------
// GrowthForce::BoundaryDilation
//
// Imposes tissue-level growth by dilating the patch boundary about its
// centroid at a prescribed relative rate:
//
//   dx/dt = g (x - x_centre)      for every vertex on the tissue boundary
//
// overwriting whatever other reactions contributed to those vertices, so the
// reaction must be listed last in the model file (the same convention
// VertexNoUpdateFromIndex relies on). A vertex is on the boundary when any of
// its walls faces the background.
//
// Without this, a free-boundary patch has nowhere for lobing to come from.
// Every interior wall is held taut between two tri-cellular junctions by two
// equally turgid cells, so a straight wall is a stable equilibrium; any
// excess wall length is exported to the free rim, which simply bulges
// outwards. That is exactly what a first attempt here did -- interior cells
// stayed perfectly polygonal for 48 simulated hours while the outer ring
// ballooned. Real leaf epidermis is confined by the surrounding tissue, and
// Sapala et al. (2018) impose the same thing, displacing wall segments
// according to a specified tissue growth.
//
// With the frame imposed, the area available to the interior cells is set by
// the tissue, so any anticlinal wall length synthesised beyond what a smooth
// outline of that area can hold has to go somewhere: it buckles.
//
// g is a LINEAR rate, so areal tissue growth is 2g. For the measured 1.9-fold
// areal growth of cotyledon pavement cells over 48 h (Higaki et al. 2021),
// g = ln(1.9) / (2 x 48) = 0.0067 / h.
//
// The two-parameter form gives separate rates along x and y, so the same
// areal growth can be delivered isotropically or anisotropically. That is the
// control for the central correlation Sapala et al. (2018) report: growth
// anisotropy and lobeyness are negatively correlated across cotyledon cells
// (r = -0.46), so a tissue stretched strongly along one axis should lobe less
// than one growing isotropically at the same areal rate.
class GrowthForceBoundaryDilation : public Reaction {
public:
  GrowthForceBoundaryDilation(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1 && p.size() != 2)
      throw std::runtime_error(
          "GrowthForce::BoundaryDilation: uses one parameter (g_linear, "
          "isotropic) or two (g_x, g_y).");
    std::vector<std::string> ids =
        p.size() == 1 ? std::vector<std::string>{"g_linear"}
                      : std::vector<std::string>{"g_x", "g_y"};
    configure("GrowthForce::BoundaryDilation", p, i, p.size(), {},
              std::move(ids));
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    const size_t dim = vertexData.cols();
    boundary_.clear();
    for (size_t v = 0; v < T.numVertex(); ++v) {
      bool onBoundary = false;
      for (size_t w : T.vertex(v).walls)
        if (Tissue::isBackground(T.wall(w).cell1) ||
            Tissue::isBackground(T.wall(w).cell2))
          onBoundary = true;
      if (onBoundary)
        boundary_.push_back(v);
    }
    centre_.assign(dim, 0.0);
    if (boundary_.empty())
      throw std::runtime_error(
          "GrowthForce::BoundaryDilation: the tissue has no boundary "
          "vertices (every wall is interior).");
    for (size_t v : boundary_)
      for (size_t d = 0; d < dim; ++d)
        centre_[d] += vertexData[v][d];
    for (size_t d = 0; d < dim; ++d)
      centre_[d] /= static_cast<double>(boundary_.size());
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dim = vertexData.cols();
    const bool perAxis = numParameter() == 2;
    for (size_t v : boundary_)
      for (size_t d = 0; d < dim; ++d) {
        const double g = perAxis && d < 2 ? parameter(d) : parameter(0);
        vertexDerivs[v][d] = g * (vertexData[v][d] - centre_[d]); // overwrite
      }
  }

private:
  std::vector<size_t> boundary_;
  std::vector<double> centre_;
};
TISSUE_REGISTER_REACTION(GrowthForceBoundaryDilation,
                         "GrowthForce::BoundaryDilation")

} // namespace
} // namespace tissue
