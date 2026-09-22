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
#include <array>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
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
// CMT::FromFile
//
// The same microtubule occupancy `m`, but read from a file instead of
// computed from a cue. This is the hook that lets an external microtubule
// simulation drive the mechanics: tubulaton runs inside the cell geometry,
// its filament array is reduced to a per-wall occupancy, and that is written
// here. Everything downstream - WallMechanics::SpringModulated,
// WallGrowth::StrainWallInhibited - is unchanged, which is the point. Swap
// CMT::CurvatureRecruitment for CMT::FromFile and the only thing that differs
// between the two models is where `m` came from.
//
// That is what makes an alignment rule a variable rather than a rewrite: the
// rules compared in ../../../../ALIGNMENT.md all end at this same variable,
// so any difference in the result is a difference in the rule and not in the
// hundred other things a second model would change.
//
// The file is called `cmt` and is read from the working directory. A fixed
// name is legacy's convention for this (SisterVertex::InitiateFromFile reads
// `sister` the same way) because a model rule carries numbers, not strings;
// in practice every run already has its own directory.
//
// Format - whitespace, `#` comments, and deliberately dull:
//
//     # anything after a hash is ignored
//     <count> <columns>
//     <value> [...]        x count
//
// Column 0 is the occupancy. `count` must equal the number of walls, and the
// reaction says so if it does not rather than reading whatever is there: a
// silently mismatched file would put one wall's microtubules on another's.
//
//     CMT::FromFile 1 1 1
//       reread_flag   # 0 = read once before the run; 1 = re-read every update
//       m_index       # the wall variable to write
//
// reread_flag 1 is for the alternating loop, where tubulaton rewrites `cmt`
// between growth steps and tissue must pick it up without restarting.
//
class CmtFromFile : public Reaction {
public:
  CmtFromFile(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "CMT::FromFile: level 0 = the wall variable to write the "
          "microtubule occupancy into.");
    configure("CMT::FromFile", p, i, 1, {1}, {"reread_flag"});
  }

  void initiate(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
                Matrix &, Matrix &) override {
    load(T, wallData);
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &, Matrix &wallData, Matrix &,
              double) override {
    if (parameter(0) != 0.0)
      load(T, wallData);
  }

private:
  void load(Tissue &T, Matrix &wallData) {
    std::ifstream in("cmt");
    if (!in)
      throw std::runtime_error(
          "CMT::FromFile: cannot open 'cmt' in the working directory. The "
          "filename is fixed, as for SisterVertex::InitiateFromFile.");
    // Strip comments first, then read numbers: simpler than interleaving the
    // two, and the header is just the first two numbers in the file.
    std::ostringstream body;
    std::string line;
    while (std::getline(in, line)) {
      const size_t hash = line.find('#');
      body << (hash == std::string::npos ? line : line.substr(0, hash)) << ' ';
    }
    std::istringstream num(body.str());
    double countValue = 0.0, columnValue = 0.0;
    if (!(num >> countValue >> columnValue))
      throw std::runtime_error(
          "CMT::FromFile: 'cmt' has no '<count> <columns>' header.");
    const size_t count = static_cast<size_t>(countValue);
    const size_t columns = static_cast<size_t>(columnValue);
    if (count != T.numWall())
      throw std::runtime_error(
          "CMT::FromFile: 'cmt' describes " + std::to_string(count) +
          " walls but the tissue has " + std::to_string(T.numWall()) +
          ". Refusing to read it: a mismatched file would put one wall's "
          "microtubules on another's.");
    if (columns == 0)
      throw std::runtime_error("CMT::FromFile: 'cmt' declares zero columns.");
    const size_t mIndex = variableIndex(0, 0);
    for (size_t w = 0; w < count; ++w) {
      double first = 0.0;
      for (size_t c = 0; c < columns; ++c) {
        double x;
        if (!(num >> x))
          throw std::runtime_error(
              "CMT::FromFile: 'cmt' ends early - expected " +
              std::to_string(count * columns) + " values.");
        if (c == 0)
          first = x;
      }
      wallData[w][mIndex] = first;
    }
  }

};
TISSUE_REGISTER_REACTION(CmtFromFile, "CMT::FromFile")

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
// MEASURED CAVEAT: in the reference pavement model this term does almost
// nothing -- beta = 0 costs 0.9% of the neck count (5.58 -> 5.53) and 0.07%
// of lobeyness over 48 h, and with the growth-inhibition channel also off,
// nothing at all (2.95 necks either way).
//
// The reason is arithmetic, not physical, and it is specific to how that run
// is configured. It uses WallGrowth::StrainWallInhibited in ACTIVE mode
// (growth_mode = 0, drive = 1.0) with the gate disabled (strain_th = -1), so
// the growth rate is k L / (1 + gamma m) and the elastic strain -- the only
// quantity beta affects -- appears nowhere in it. The direct path from
// stiffness to growth is closed by construction.
//
// What remains is an indirect geometric path: stiffness changes the
// instantaneous elastic shape, which changes local curvature, which changes m
// through CMT::CurvatureRecruitment, which changes growth. The two controls
// separate the contributions. With gamma = 0, m cannot reach growth at all
// and beta moves the neck count by 0.0% and lobeyness by 0.03% -- the pure
// elastic-shape effect. With gamma = 9 the curvature path opens and the neck
// count moves by 0.9%. So the geometric path is real but weak.
//
// Do NOT generalise this to "wall stiffness does not matter". In a
// strain-gated model the direct path is open, and the effect is then roughly
// the term's share of the force balance: on the apical-hook model in this
// repository Y_fiber is ~6% of wall stiffness and contributes ~4.5% of the
// opening, under both RK5Adaptive and QuasiStatic. Here beta m raises wall
// stiffness by ~160% at typical m and still buys 0.9% -- because the path is
// shut, not because the term is small.
//
// An earlier version of this comment blamed drag lag. That was wrong: the
// apical-hook model's Y_fiber effect survives QuasiStatic, where there is no
// lag by construction. Lag IS present in this model -- a rate-scaling sweep
// (all rates x s, t_end / s) has RK5Adaptive converging on QuasiStatic as
// growth slows, with the lobeyness gap going 14.7% -> 9.6% -> -0.3% for
// s = 1, 1/2, 1/4, against a QuasiStatic result that moves by under 0.3%
// across the same range. So lag is real and measurable here; it simply is
// not the channel through which stiffness acts.
//
// Keep the term: it is physically real, and it is the right place for
// cellulose stiffening in any variant with strain-gated growth. But a fitted
// beta in the reference configuration is not evidence about cellulose
// stiffness -- those observables cannot constrain it.
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
// Armour et al. (2015) Plant Cell 27:2484. With an independent synthesis
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
// documented: Armour, Barton, Law & Overall (2015) Plant Cell 27:2484-2500,
// "Differential growth in periclinal and anticlinal walls during lobe
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
    // A per-face mesh gives every wall cell2 = background, so the test below
    // would call every vertex a boundary vertex and this reaction would
    // overwrite every derivative -- the tissue then only dilates rigidly and
    // looks frozen. In such a mesh an interior vertex is one that has a
    // sister; the true tissue boundary is the vertices with none.
    std::vector<char> hasSister(T.numVertex(), 0);
    for (size_t i = 0; i < T.numSisterVertex(); ++i) {
      hasSister[T.sisterVertex(i, 0)] = 1;
      hasSister[T.sisterVertex(i, 1)] = 1;
    }
    for (size_t v = 0; v < T.numVertex(); ++v) {
      if (hasSister[v])
        continue;
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

  // This is a prescribed velocity, not a force: g (x - centre) never vanishes,
  // so a force-balance solver can never satisfy it and instead accelerates
  // the boundary outward without bound. (Measured before the interface
  // existed: cell area 296 -> 2.8e6 um^2 over 48 h under QuasiStatic, with
  // every growth step hitting the relaxation cap.) Declaring it here lets
  // QuasiStatic integrate this part over the growth step and hold it out of
  // the relaxation; explicit solvers are unaffected.
  bool prescribesVelocity() const override { return true; }
  void velocityDerivs(Tissue &T, Matrix &cellData, Matrix &wallData,
                      Matrix &vertexData, Matrix &vertexVel) override {
    Matrix ignoredCell, ignoredWall;
    ignoredCell.reshapeLike(cellData);
    ignoredWall.reshapeLike(wallData);
    derivs(T, cellData, wallData, vertexData, ignoredCell, ignoredWall,
           vertexVel);
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


// ---------------------------------------------------------------------------
// WallMechanics::Bending
//
// Wall bending stiffness with the correct continuum limit, so that the
// selected wavelength is set by the physics and not by the mesh.
//
// WallMechanics::BendingChain applies F_v = k (midpoint - x_v), a discrete
// Laplacian stencil, conservative for E = sum_v (k/2)|s_v|^2. For a sinusoid
// of amplitude A and wavelength lambda at vertex spacing h that energy is
//
//   E / length = k pi^4 A^2 h^3 / lambda^4
//
// (verified against the exact discrete sum: the ratio per halving of h
// converges on 1/8, and the closed form agrees to four significant figures).
//
// So the WAVELENGTH exponent is already correct -- it penalises short
// wavelengths as lambda^-4, exactly as a bending term should. The defect is
// the h^3 prefactor: the effective stiffness is proportional to the cube of
// the vertex spacing, so halving the mesh weakens bending EIGHTFOLD.
//
// (An earlier version of this comment said the energy went as h/lambda^2 and
// blamed a wrong wavelength exponent. Both were wrong: that probe summed
// |F_v| A, a force density, which does scale as h^1. The remedy below was
// derived from the energy and is unaffected -- B = k h^3 / 4 carries the
// correct h^3 -- but the size of the effect is 8x per halving, not 2x.)
//
// That is not academic. Refining this model's wall discretisation from 0.6 to
// 0.3 um, on an identical initial outline, moves lobeyness from 1.295 to
// 1.127 (-12.9%) and the neck count from 5.29 to 10.14 per cell (+91.9%):
// weaker effective bending, shorter selected wavelength, nearly twice as many
// lobes. The lobe spacing was largely a property of the mesh. Since
// lambda ~ B^(1/4), an eightfold bending change predicts 8^(1/4) = 1.68x more
// necks, against 1.92x measured -- the same order, so bending plausibly is
// what sets the wavelength here.
//
// "On an identical initial outline" is load-bearing. The first version of
// this study seeded the initial perturbation with a per-vertex random draw of
// fixed amplitude, whose SLOPE therefore scales as 1/h: refining the mesh
// silently roughened the starting condition (total initial wall length
// 1.000 -> 1.019 -> 1.091 over a 0.6/0.3/0.15 um series), making it a
// two-variable comparison. It reported +30.7% where the single-variable
// answer is +91.9%, i.e. the confound was hiding two thirds of the effect.
//
// A real bending energy for a chain is
//
//   E = sum_v (B/2) kappa_v^2 l_v,   kappa_v = 2 |s_v| / h^2,   l_v = h
//     = sum_v 2 B |s_v|^2 / h^3
//
// with s_v = (x_a + x_b)/2 - x_v the sagitta vector. Substituting a sinusoid
// gives an energy per unit length proportional to B A^2 / lambda^4 with no h
// in it -- the same fourth-power law BendingChain already has, but now with a
// mesh-independent prefactor.
//
// The force is the same three-point stencil BendingChain uses, which is why
// this is a small change rather than a new operator:
//
//   F_v = 4 B s_v / h^3,   F_a = F_b = -F_v / 2
//
// so it is exactly BendingChain with k replaced by 4B/h^3, evaluated per
// vertex from the current spacing. B is a bending modulus (force x length^2)
// and is a material property, unlike k_bend which had to be retuned for every
// discretisation.
//
// Measured, same 0.6 -> 0.3 um refinement on an identical initial outline:
// lobeyness drift -2.5% against BendingChain's -12.9%, neck-count drift
// +12.2% against +91.9%. Roughly 5x better on lobeyness and 7.5x on the neck
// count, which is the wavelength observable and the one that matters.
//
// NOT converged, though -- 12.2% over a single halving is much better than
// 92% and is still not convergence. The residual could be the neck-tip
// self-intersection (this model has no self-avoidance), the curvature
// estimator, or simply needing h = 0.15.
//
// BendingChain is left alone: other models may be calibrated against it, and
// B = k_bend * h^3 / 4 converts a calibrated k to the equivalent modulus.
class WallMechanicsBending : public Reaction {
public:
  WallMechanicsBending(const ParameterList &p, const IndexLevels &i) {
    configure("WallMechanics::Bending", p, i, 1, {1}, {"B_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t flagIndex = variableIndex(0, 0);
    const double bBend = parameter(0);
    const size_t dim = vertexData.cols();
    parallelScatter1(
        T.numVertex(), vertexDerivs, [&](size_t begin, size_t end,
                                         Matrix &out) {
          for (size_t v = begin; v < end; ++v) {
            size_t chain[2];
            size_t found = 0;
            for (size_t w : T.vertex(v).walls) {
              if (wallData[w][flagIndex] != 0.0) {
                if (found < 2)
                  chain[found] = w;
                ++found;
              }
            }
            if (found != 2)
              continue; // chain ends and junctions carry no bending force
            const size_t a = T.wall(chain[0]).otherVertex(v);
            const size_t b = T.wall(chain[1]).otherVertex(v);
            double da = 0.0, db = 0.0;
            for (size_t d = 0; d < dim; ++d) {
              da += (vertexData[a][d] - vertexData[v][d]) *
                    (vertexData[a][d] - vertexData[v][d]);
              db += (vertexData[b][d] - vertexData[v][d]) *
                    (vertexData[b][d] - vertexData[v][d]);
            }
            const double h = 0.5 * (std::sqrt(da) + std::sqrt(db));
            if (h <= 0.0)
              continue;
            const double k = 4.0 * bBend / (h * h * h);
            for (size_t d = 0; d < dim; ++d) {
              const double f =
                  k * (0.5 * (vertexData[a][d] + vertexData[b][d]) -
                       vertexData[v][d]);
              out[v][d] += f;
              out[a][d] -= 0.5 * f;
              out[b][d] -= 0.5 * f;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(WallMechanicsBending, "WallMechanics::Bending")


// ---------------------------------------------------------------------------
// Differential wall synthesis across a shared wall
//
// Why this exists. With one growth rule per wall, a vertex model can only
// deepen a lobe by taking material from the neck: the cell's area is set by
// turgor against an area-elastic term whose resting area grows at the same
// rate for every cell (dA0/dt = k_area A0), so no cell can advance into its
// neighbour and the perimeter must buckle inward within a fixed budget.
// Measured against 370 tracked Arabidopsis pavement cells, that model
// reproduces circularity (0.399 against 0.397), solidity (0.734 against
// 0.747) and lobe number (4 against 5) at the end of morphogenesis, and
// misses minimum neck width by a factor of five (0.054 against 0.279) -- the
// one shape statistic of the four that cannot be written in terms of
// perimeter and area, and so the one a perimeter-and-area mechanism cannot
// be expected to get right.
//
// What a real cell does instead is deposit wall material on one face of a
// shared wall and not the other. The two cells do not end up with different
// lengths of the same wall -- there is only one wall -- they end up with a
// wall that is *bent*, because material added to the outer face lengthens the
// outer arc. So differential synthesis appears in a vertex model not as a
// second resting length but as a non-zero rest curvature.
//
// The consequence is the point. An elastic wall bent away from a flat rest
// state pulls back, and that tension is what a lobe must be held out against
// -- paid for by the neck. A wall whose rest state *is* the bend exerts no
// such pull, so the lobe is held out by the material in it rather than by
// tension stolen from elsewhere.
//
// Rest curvature is stored per wall and is signed relative to the wall's
// cell1: positive means bowed away from cell1, into cell2. Both reactions
// below take the reference direction from cell1's centroid, so the sign is
// unambiguous and the two agree.
// ---------------------------------------------------------------------------

namespace {

// Geometry at an interior chain vertex: the sagitta, the spacing, and a unit
// reference direction pointing from cell1's centroid towards the vertex.
// Returns false at chain ends, junctions, and walls with no cell on either
// side, none of which carry bending.
struct ChainGeom {
  double sDotU = 0.0;   // sagitta projected on the reference direction
  double h = 0.0;       // mean spacing to the two neighbours
  double u[3] = {0, 0, 0};
  double s[3] = {0, 0, 0};
  size_t a = 0, b = 0;
};

inline bool chainGeometry(const Tissue &T, const Matrix &wallData,
                          const Matrix &vertexData, size_t v,
                          size_t flagIndex,
                          const std::vector<std::array<double, 3>> &centre,
                          ChainGeom &g) {
  size_t chain[2];
  size_t found = 0;
  for (size_t w : T.vertex(v).walls) {
    if (wallData[w][flagIndex] != 0.0) {
      if (found < 2)
        chain[found] = w;
      ++found;
    }
  }
  if (found != 2)
    return false;
  size_t ref = T.wall(chain[0]).cell1;
  if (ref == kBackground)
    ref = T.wall(chain[0]).cell2;
  if (ref == kBackground)
    return false;

  g.a = T.wall(chain[0]).otherVertex(v);
  g.b = T.wall(chain[1]).otherVertex(v);
  const size_t dim = vertexData.cols();
  double da = 0.0, db = 0.0, un = 0.0;
  for (size_t d = 0; d < dim; ++d) {
    const double ea = vertexData[g.a][d] - vertexData[v][d];
    const double eb = vertexData[g.b][d] - vertexData[v][d];
    da += ea * ea;
    db += eb * eb;
    g.s[d] = 0.5 * (vertexData[g.a][d] + vertexData[g.b][d]) - vertexData[v][d];
    g.u[d] = vertexData[v][d] - centre[ref][d];
    un += g.u[d] * g.u[d];
  }
  g.h = 0.5 * (std::sqrt(da) + std::sqrt(db));
  if (g.h <= 0.0 || un <= 0.0)
    return false;
  un = std::sqrt(un);
  g.sDotU = 0.0;
  for (size_t d = 0; d < dim; ++d) {
    g.u[d] /= un;
    g.sDotU += g.s[d] * g.u[d];
  }
  return true;
}

inline std::vector<std::array<double, 3>>
cellCentres(const Tissue &T, const Matrix &vertexData) {
  std::vector<std::array<double, 3>> c(T.numCell(), {0.0, 0.0, 0.0});
  for (size_t i = 0; i < T.numCell(); ++i) {
    const Vec3 p = T.cellPosition(i, vertexData);
    c[i] = {p[0], p[1], p[2]};
  }
  return c;
}

} // namespace

// ---------------------------------------------------------------------------
// WallMechanics::BendingSpontaneous
//
// WallMechanics::Bending with a rest curvature read from a wall variable
// instead of assumed zero. The energy is
//
//   E = sum_v (B/2) (kappa_v - kappa0_v)^2 l_v
//
// which, with kappa = 2 |s| / h^2 as in Bending, gives the same three-point
// stencil displaced by a rest sagitta s0 = (kappa0 h^2 / 2) u:
//
//   F_v = 4 B (s_v - s0) / h^3,   F_a = F_b = -F_v / 2
//
// At kappa0 = 0 it is Bending exactly, so a model can be switched between
// them to isolate what the rest curvature does.
class WallMechanicsBendingSpontaneous : public Reaction {
public:
  WallMechanicsBendingSpontaneous(const ParameterList &p,
                                  const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "WallMechanics::BendingSpontaneous: level 0 = wall bend flag, "
          "level 1 = wall rest-curvature index.");
    configure("WallMechanics::BendingSpontaneous", p, i, 1, {1, 1},
              {"B_bend"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t flagIndex = variableIndex(0, 0);
    const size_t k0Index = variableIndex(1, 0);
    const double bBend = parameter(0);
    const size_t dim = vertexData.cols();
    const auto centre = cellCentres(T, vertexData);

    parallelScatter1(
        T.numVertex(), vertexDerivs,
        [&](size_t begin, size_t end, Matrix &out) {
          ChainGeom g;
          for (size_t v = begin; v < end; ++v) {
            if (!chainGeometry(T, wallData, vertexData, v, flagIndex, centre,
                               g))
              continue;
            // rest curvature at the vertex: mean over its two chain walls
            double k0 = 0.0;
            size_t n = 0;
            for (size_t w : T.vertex(v).walls)
              if (wallData[w][flagIndex] != 0.0) {
                k0 += wallData[w][k0Index];
                ++n;
              }
            if (n)
              k0 /= double(n);
            const double k = 4.0 * bBend / (g.h * g.h * g.h);
            const double s0 = 0.5 * k0 * g.h * g.h;
            for (size_t d = 0; d < dim; ++d) {
              const double f = k * (g.s[d] - s0 * g.u[d]);
              out[v][d] += f;
              out[g.a][d] -= 0.5 * f;
              out[g.b][d] -= 0.5 * f;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(WallMechanicsBendingSpontaneous,
                         "WallMechanics::BendingSpontaneous")

// ---------------------------------------------------------------------------
// WallGrowth::CurvaturePlasticity
//
// The rest curvature follows the actual curvature, with a rate:
//
//   dkappa0/dt = k_plastic (amplify * kappa - kappa0) / (1 + gamma_m m)
//
// amplify = 1 is pure plasticity: the rest curvature relaxes to whatever
// curvature the wall is held at, the elastic restoring moment decays to
// zero, and the wall neither pulls back nor pushes on. That turns out not to
// be enough. Removing the restoring moment lets a lobe stay where it is
// without the neck paying to hold it, but it does not make the lobe advance,
// and the arc-length budget is untouched -- so the neck still pays for
// getting there. Measured, it makes necks very slightly thinner, because a
// wall with no restoring moment is also a wall that buckles more easily.
//
// amplify > 1 is the part that matters. Material laid down on the convex
// face of a wall that is already bent does not merely record the bend, it
// drives it: the rest curvature overshoots the actual curvature and the wall
// carries an active bending moment pushing it further. That is what makes a
// lobe advance into the neighbour rather than be carved out of the cell, and
// it is the half of differential synthesis a relaxation term cannot express.
//
// This is what differential synthesis does, stated at the level a vertex
// model can represent. A wall that is held bent has its outer face longer
// than its inner one; material deposited there while it is bent takes on that
// shape, and the bend stops being elastic strain and becomes the wall's
// resting form. It is the bending analogue of the irreversible extension
// WallGrowth::StrainWallInhibited already applies to length, and it is
// inhibited by the same reinforcement variable for the same reason: where
// microtubules and the cellulose they guide are dense, the wall yields less.
//
// Setting k_plastic = 0 recovers a purely elastic wall, so the mechanism can
// be switched off without changing anything else in the model.
class WallGrowthCurvaturePlasticity : public Reaction {
public:
  WallGrowthCurvaturePlasticity(const ParameterList &p, const IndexLevels &i) {
    if (p.size() < 1 || p.size() > 3)
      throw std::runtime_error(
          "WallGrowth::CurvaturePlasticity: uses one to three parameters "
          "(k_plastic, [gamma_m, default 0], [amplify, default 1]).");
    if (p.size() == 3 && p[2] < 0.0)
      throw std::runtime_error(
          "WallGrowth::CurvaturePlasticity: amplify (parameter 3) must be "
          "non-negative.");
    if (p[0] < 0.0)
      throw std::runtime_error(
          "WallGrowth::CurvaturePlasticity: k_plastic must be non-negative.");
    if (i.size() < 2 || i[0].size() != 1 || i[1].size() != 1 ||
        (i.size() == 3 && i[2].size() != 1) || i.size() > 3)
      throw std::runtime_error(
          "WallGrowth::CurvaturePlasticity: level 0 = wall bend flag, "
          "level 1 = wall rest-curvature index, "
          "optional level 2 = wall reinforcement index.");
    std::vector<std::string> ids{"k_plastic"};
    if (p.size() >= 2)
      ids.push_back("gamma_m");
    if (p.size() >= 3)
      ids.push_back("amplify");
    const std::vector<size_t> counts =
        i.size() == 3 ? std::vector<size_t>{1, 1, 1}
                      : std::vector<size_t>{1, 1};
    configure("WallGrowth::CurvaturePlasticity", p, i, p.size(), counts,
              std::move(ids));
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t flagIndex = variableIndex(0, 0);
    const size_t k0Index = variableIndex(1, 0);
    const bool useM = numVariableIndexLevel() == 3;
    const size_t mIndex = useM ? variableIndex(2, 0) : 0;
    const double kPlastic = parameter(0);
    const double gammaM = numParameter() >= 2 ? parameter(1) : 0.0;
    const double amplify = numParameter() >= 3 ? parameter(2) : 1.0;
    if (kPlastic == 0.0)
      return;
    const auto centre = cellCentres(T, vertexData);

    // Walked per wall, not per vertex, so each wall is written exactly once
    // and the loop needs no synchronisation.
    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      ChainGeom g;
      for (size_t w = begin; w < end; ++w) {
        if (wallData[w][flagIndex] == 0.0)
          continue;
        double kappa = 0.0;
        size_t n = 0;
        for (size_t v : {T.wall(w).vertex1, T.wall(w).vertex2}) {
          if (!chainGeometry(T, wallData, vertexData, v, flagIndex, centre, g))
            continue;
          // kappa = 2 (s . u) / h^2, signed: positive bows away from cell1
          kappa += 2.0 * g.sDotU / (g.h * g.h);
          ++n;
        }
        if (!n)
          continue;
        kappa /= double(n);
        const double m = useM ? wallData[w][mIndex] : 0.0;
        wallDerivs[w][k0Index] += kPlastic *
            (amplify * kappa - wallData[w][k0Index]) / (1.0 + gammaM * m);
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthCurvaturePlasticity,
                         "WallGrowth::CurvaturePlasticity")

// ---------------------------------------------------------------------------
// WallMechanics::ActiveWallExpansion
//
// A turgor-independent outward force on the wall, set by a signed per-wall
// variable and directed along the wall normal.
//
// Why this is not the rest-curvature reaction above. Haas et al. (2020,
// Science 367:1003) show that de-esterified homogalacturonan nanofilaments
// expand, and that the wall reshapes an epidermal cell with the turgor
// removed -- the wall is not a passive material being inflated, it pushes.
// Rest curvature expresses the geometry an asymmetric swelling would leave
// behind: it is a couple, and a couple bends a wall in place. Swelling is a
// force, and a force translates the wall, so the lobe advances into the
// neighbour instead of being folded out of the cell's own boundary. Measured
// here, the couple makes necks 3-12% thinner and the whole point is that it
// never moved material outward.
//
// The force on a segment is
//
//   F = k_push * p * l * n,    n the unit normal, l the segment length
//
// split between its two vertices, with n pointing away from cell1 so a
// positive p drives the wall into cell2. Like turgor, it is unbalanced on a
// single wall and balanced over a cell; unlike turgor it is local, signed,
// and does not scale with the cell's area.
//
// p is read rather than computed, so the patterning is a separate question
// from the mechanics. Pointing it at the variable
// WallGrowth::CurvaturePlasticity maintains gives a swelling that follows
// the wall's own curvature, which is the geometric feedback; pointing it at
// a variable driven by anything else -- a polarity field, a prescribed
// pattern, a pectin chemistry -- tests that instead without touching this
// reaction.
class WallMechanicsActiveWallExpansion : public Reaction {
public:
  WallMechanicsActiveWallExpansion(const ParameterList &p,
                                   const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "WallMechanics::ActiveWallExpansion: level 0 = the signed per-wall "
          "variable driving the expansion (positive pushes into cell2).");
    if (p.size() < 1 || p.size() > 3)
      throw std::runtime_error(
          "WallMechanics::ActiveWallExpansion: uses one to three parameters "
          "(k_push, [p_offset, default 0], [subtract_cell_mean, default 0]). "
          "The offset is subtracted from the driver before use, so a variable "
          "that is positive everywhere -- a morphogen about its steady state, "
          "say -- can drive a signed push. subtract_cell_mean = 1 uses the "
          "mean of the driver over the cell's own walls instead, which "
          "matters more than it sounds: a normal force applied around a "
          "closed loop with a non-zero mean mostly INFLATES the loop, exactly "
          "like extra turgor, and only the variation along it changes shape. "
          "Measured with a fixed offset, the push grew to dominate every "
          "other force (1.36x their sum) while correlating with neither wall "
          "strain (-0.05) nor wall growth (-0.06): almost all of it went into "
          "isotropic inflation that the area-elastic term then fought.");
    std::vector<std::string> ids{"k_push"};
    if (p.size() >= 2)
      ids.push_back("p_offset");
    if (p.size() >= 3)
      ids.push_back("subtract_cell_mean");
    configure("WallMechanics::ActiveWallExpansion", p, i, p.size(), {1},
              std::move(ids));
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "WallMechanics::ActiveWallExpansion requires 2D.");
    const size_t pIndex = variableIndex(0, 0);
    const double kPush = parameter(0);
    const double pOff = numParameter() >= 2 ? parameter(1) : 0.0;
    const bool perCellMean = numParameter() >= 3 && parameter(2) != 0.0;
    if (kPush == 0.0)
      return;
    const auto centre = cellCentres(T, vertexData);

    // Per-cell mean of the driver, so the push has no net inflating part.
    std::vector<double> cellMean;
    if (perCellMean) {
      cellMean.assign(T.numCell(), 0.0);
      std::vector<double> len(T.numCell(), 0.0);
      for (size_t w = 0; w < T.numWall(); ++w) {
        const double l = T.wallLengthFromVertices(w, vertexData);
        for (size_t c : {T.wall(w).cell1, T.wall(w).cell2}) {
          if (c == kBackground)
            continue;
          cellMean[c] += wallData[w][pIndex] * l;
          len[c] += l;
        }
      }
      for (size_t c = 0; c < T.numCell(); ++c)
        if (len[c] > 0.0)
          cellMean[c] /= len[c];
    }

    parallelScatter1(
        T.numWall(), vertexDerivs, [&](size_t begin, size_t end, Matrix &out) {
          for (size_t w = begin; w < end; ++w) {
            size_t ref = T.wall(w).cell1;
            if (ref == kBackground)
              continue;   // a free edge has no neighbour to advance into
            const double p = wallData[w][pIndex] -
                             (perCellMean ? cellMean[ref] : pOff);
            if (p == 0.0)
              continue;
            const size_t v1 = T.wall(w).vertex1;
            const size_t v2 = T.wall(w).vertex2;
            const double tx = vertexData[v2][0] - vertexData[v1][0];
            const double ty = vertexData[v2][1] - vertexData[v1][1];
            const double len = std::sqrt(tx * tx + ty * ty);
            if (len <= 0.0)
              continue;
            // Normal, oriented away from cell1's centroid.
            double nx = -ty / len, ny = tx / len;
            const double mx = 0.5 * (vertexData[v1][0] + vertexData[v2][0]);
            const double my = 0.5 * (vertexData[v1][1] + vertexData[v2][1]);
            if (nx * (mx - centre[ref][0]) + ny * (my - centre[ref][1]) < 0.0) {
              nx = -nx;
              ny = -ny;
            }
            const double f = 0.5 * kPush * p * len;
            out[v1][0] += f * nx;
            out[v1][1] += f * ny;
            out[v2][0] += f * nx;
            out[v2][1] += f * ny;
          }
        });
  }
};
TISSUE_REGISTER_REACTION(WallMechanicsActiveWallExpansion,
                         "WallMechanics::ActiveWallExpansion")

// ---------------------------------------------------------------------------
// WallPattern::Turing
//
// An activator-inhibitor system living on the wall chain, so that a cell can
// pattern its own boundary without reference to the boundary's shape.
//
// Why. Every spatial pattern in this model currently comes from curvature:
// the reinforcement is recruited by it, the rest curvature follows it, the
// active expansion reads it, and nothing is upstream of it. Three mechanisms
// built to widen necks -- a passive couple, an active couple, and an active
// force -- all made them narrower, and what they share is that driver. A
// curvature-driven feedback amplifies the curvature already present, which
// sharpens a feature instead of broadening it, and sharpening is the failure.
//
// In a real cell the pattern does not come from the shape. De-esterified
// homogalacturonan is enriched at necks *before* the cortical microtubules
// are, and ROP2/ROP6 mutual antagonism is the documented candidate for
// placing the domains. Both are reaction-diffusion in character: local
// self-activation with longer-range inhibition, which selects a wavelength
// from kinetics and diffusion rather than from the geometry it then imposes.
//
// Schnakenberg kinetics, which is the smallest system that Turing-patterns:
//
//   du/dt = Du L(u) + a - u + u^2 v
//   dv/dt = Dv L(v) + b     - u^2 v
//
// with L the discrete Laplacian along the chain of walls. Patterns need
// Dv >> Du; the selected wavelength then scales as sqrt(Dv), which is the
// knob that sets lobe spacing here, in place of the bending modulus.
//
// Junction handling, which turned out to matter more than the kinetics.
//
// With no flux at junctions (junction_flux = 0) the pattern is not free: a
// no-flux boundary pins the phase, and the measured result is an activator
// depleted at junctions (u = 0.75 within 1 um) and peaking at mid-wall
// (u = 1.23 at 6-10 um, half the selected wavelength away). The pattern
// then places its peaks exactly where the cell's existing junction geometry
// dictates, which defeats the purpose of using a driver independent of the
// geometry.
//
// junction_flux > 0 couples every wall meeting at a vertex to the mean of
// the others, so the pattern can choose its own phase across a junction.
// This blurs cell identity, because a wall is shared and one variable cannot
// carry two cells' morphogens -- the same shared-wall limit that defeated
// expressing differential synthesis as two resting lengths. It is a
// simplification, and the honest version needs a per-cell-side variable.
class WallPatternTuring : public Reaction {
public:
  WallPatternTuring(const ParameterList &p, const IndexLevels &i) {
    if (p.size() < 4 || p.size() > 5)
      throw std::runtime_error(
          "WallPattern::Turing: uses four or five parameters "
          "(D_u, D_v, a, b, [junction_flux, default 0]).");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "WallPattern::Turing: level 0 = activator wall index, "
          "level 1 = inhibitor wall index.");
    std::vector<std::string> ids{"D_u", "D_v", "a", "b"};
    if (p.size() == 5)
      ids.push_back("junction_flux");
    configure("WallPattern::Turing", p, i, p.size(), {1, 1},
              std::move(ids));
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t uI = variableIndex(0, 0);
    const size_t vI = variableIndex(1, 0);
    const double Du = parameter(0), Dv = parameter(1);
    const double a = parameter(2), b = parameter(3);
    const double jflux = numParameter() == 5 ? parameter(4) : 0.0;

    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        const double u = wallData[w][uI], v = wallData[w][vI];
        double lapU = 0.0, lapV = 0.0;
        // neighbours: the other wall at each of this wall's vertices, when
        // that vertex joins exactly two walls of the same chain
        for (size_t vtx : {T.wall(w).vertex1, T.wall(w).vertex2}) {
          const auto &ws = T.vertex(vtx).walls;
          if (ws.size() == 2) {
            size_t other = kBackground;
            for (size_t x : ws)
              if (x != w)
                other = x;
            if (other == kBackground)
              continue;
            const double h =
                0.5 * (T.wallLengthFromVertices(w, vertexData) +
                       T.wallLengthFromVertices(other, vertexData));
            if (h <= 0.0)
              continue;
            lapU += (wallData[other][uI] - u) / (h * h);
            lapV += (wallData[other][vI] - v) / (h * h);
          } else if (jflux > 0.0 && ws.size() > 2) {
            // junction: exchange with the mean of the other walls meeting
            // here, weighted so that jflux = 1 matches the strength of an
            // ordinary neighbour exchange
            double mu = 0.0, mv = 0.0, hsum = 0.0;
            size_t n = 0;
            for (size_t x : ws) {
              if (x == w)
                continue;
              mu += wallData[x][uI];
              mv += wallData[x][vI];
              hsum += T.wallLengthFromVertices(x, vertexData);
              ++n;
            }
            if (!n)
              continue;
            const double h =
                0.5 * (T.wallLengthFromVertices(w, vertexData) + hsum / n);
            if (h <= 0.0)
              continue;
            lapU += jflux * (mu / n - u) / (h * h);
            lapV += jflux * (mv / n - v) / (h * h);
          }
        }
        const double uuv = u * u * v;
        wallDerivs[w][uI] += Du * lapU + a - u + uuv;
        wallDerivs[w][vI] += Dv * lapV + b - uuv;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallPatternTuring, "WallPattern::Turing")

// ---------------------------------------------------------------------------
// CMT::Spread
//
// Lets the reinforcement diffuse along the wall, so that it carries a length
// scale.
//
// Every recruitment law here is pointwise: m on a wall segment is set by that
// segment's own curvature, or its own strain. A pointwise law can reinforce
// an arbitrarily sharp feature, and that is the signature of the model's
// failure -- three mechanisms built to broaden necks all sharpened them
// instead, and the shapes end up as a round body carrying thin appendages.
//
// A microtubule cannot do that. It is a stiff filament several microns long
// -- 4.4 um mean in the calibrated arrays here -- and it reinforces the wall
// along its whole length, so it cannot reinforce a neck half a micron across
// while leaving the wall on either side unreinforced. It bridges the neck.
// The filament length is a minimum feature size for reinforcement, and
// nothing in a pointwise law expresses that.
//
//   dm/dt += D_m * L(m)
//
// with L the discrete Laplacian along the chain, as in WallPattern::Turing.
// Against the turnover k_off of the recruitment reaction this sets a length
// scale sqrt(D_m / k_off), which should be set to the mean filament length
// rather than fitted: at k_off = 0.2 /h and 4.4 um, D_m = 3.9 um^2/h.
//
// This is a reduced stand-in for running the filament model in the loop. It
// captures the one property of an array that a pointwise law cannot -- that
// reinforcement has a minimum spatial extent -- and none of the others:
// orientation, bundling, severing at crossovers, or the array's own response
// to stress. If it broadens the necks, the filament length is the missing
// ingredient and the full coupling is worth its cost; if it does not, the
// coupling would be unlikely to help for this reason and the search should
// go elsewhere.
class CMTSpread : public Reaction {
public:
  CMTSpread(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1)
      throw std::runtime_error("CMT::Spread: uses one parameter (D_m).");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "CMT::Spread: level 0 = wall reinforcement index.");
    configure("CMT::Spread", p, i, 1, {1}, {"D_m"});
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t mI = variableIndex(0, 0);
    const double Dm = parameter(0);
    if (Dm == 0.0)
      return;
    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        const double m = wallData[w][mI];
        double lap = 0.0;
        for (size_t vtx : {T.wall(w).vertex1, T.wall(w).vertex2}) {
          const auto &ws = T.vertex(vtx).walls;
          if (ws.size() != 2)
            continue;   // a junction is where one filament's reach ends
          size_t other = kBackground;
          for (size_t x : ws)
            if (x != w)
              other = x;
          if (other == kBackground)
            continue;
          const double h =
              0.5 * (T.wallLengthFromVertices(w, vertexData) +
                     T.wallLengthFromVertices(other, vertexData));
          if (h <= 0.0)
            continue;
          lap += (wallData[other][mI] - m) / (h * h);
        }
        wallDerivs[w][mI] += Dm * lap;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CMTSpread, "CMT::Spread")

// ---------------------------------------------------------------------------
// CMT::SignedCurvatureRecruitment
//
// Recruitment driven by curvature *as the owning cell sees it*, so that the
// two faces of a wall can reinforce differently.
//
// On a shared-wall mesh this reaction is meaningless: one wall has one m, and
// the two cells cannot disagree about it. On a per-face mesh each wall
// belongs to exactly one cell, and the same physical wall is concave for one
// cell and convex for the other. That sign is the only input the two faces
// do not share -- geometry, length, strain and unsigned curvature are all
// identical by construction -- so it is the whole of what a per-face mesh
// buys, and without it the two faces evolve identically and the mesh is
// wasted. Measured on the first per-face run: m differed across the faces by
// exactly zero.
//
// The asymmetry is measured, not assumed. Running tubulaton in these cell
// geometries gives cortical density relative to the cell's mean of
//
//   convex (lobe tip)  0.670 +/- 0.109
//   concave (neck)     1.191 +/- 0.220
//
// so a stiff filament is depleted about 1.8x on convex cortex relative to
// concave. The shared-wall model could only represent the sum of the two
// faces, g(kappa) + g(-kappa), which is even in kappa and throws the
// asymmetry away; that is why its recruitment had to be written as a
// function of |kappa|.
//
//   dm/dt = k_on * S * (1 - m) - k_off * m
//   S     = 1 / (1 + exp(-kappa_signed / kappa_half))
//
// A logistic rather than a Hill, because the cue is signed and a Hill in
// |kappa| cannot express "more on one side than the other". S -> 1 on
// strongly concave wall, -> 0 on strongly convex, and 0.5 on flat, so
// kappa_half sets how sharply the two faces separate. Positive curvature is
// measured away from the owning cell's centroid, i.e. convex for that cell,
// so the sign is negated here to make concave the recruiting side.
class CMTSignedCurvatureRecruitment : public Reaction {
public:
  CMTSignedCurvatureRecruitment(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3)
      throw std::runtime_error(
          "CMT::SignedCurvatureRecruitment: uses three parameters "
          "(k_on, k_off, kappa_half).");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "CMT::SignedCurvatureRecruitment: level 0 = wall bend flag, "
          "level 1 = wall reinforcement index.");
    configure("CMT::SignedCurvatureRecruitment", p, i, 3, {1, 1},
              {"k_on", "k_off", "kappa_half"});
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t flagIndex = variableIndex(0, 0);
    const size_t mIndex = variableIndex(1, 0);
    const double kOn = parameter(0), kOff = parameter(1);
    const double kHalf = parameter(2);
    if (kHalf <= 0.0)
      throw std::runtime_error(
          "CMT::SignedCurvatureRecruitment: kappa_half must be positive.");
    const auto centre = cellCentres(T, vertexData);

    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      ChainGeom g;
      for (size_t w = begin; w < end; ++w) {
        if (wallData[w][flagIndex] == 0.0)
          continue;
        double kappa = 0.0;
        size_t n = 0;
        for (size_t v : {T.wall(w).vertex1, T.wall(w).vertex2}) {
          if (!chainGeometry(T, wallData, vertexData, v, flagIndex, centre, g))
            continue;
          kappa += 2.0 * g.sDotU / (g.h * g.h);
          ++n;
        }
        const double m = wallData[w][mIndex];
        double S = 0.5;
        if (n) {
          kappa /= double(n);
          S = 1.0 / (1.0 + std::exp(kappa / kHalf));   // concave recruits
        }
        wallDerivs[w][mIndex] += kOn * S * (1.0 - m) - kOff * m;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CMTSignedCurvatureRecruitment,
                         "CMT::SignedCurvatureRecruitment")

// ---------------------------------------------------------------------------
// CMT::PatternRecruitment
//
// Recruits reinforcement from a patterning variable instead of from local
// curvature, with the same first-order kinetics as
// CMT::CurvatureRecruitment so that the two are a controlled comparison:
// identical downstream mechanics, identical dm/dt, different cue.
//
// The reason for wanting that comparison is in stage3's figure 52. A
// curvature cue is pure positive feedback with no length scale of its own,
// so the number of lobes is whatever the initial condition seeded and
// growing the cell only stretches them -- measured, the model's perimeter
// per lobe rises from 1.56 to 2.85 while real cells refine from 2.50 to
// 1.39. A reaction-diffusion pattern has an intrinsic wavelength and
// inserts new peaks as its domain grows, which is the behaviour the data
// show.
//
// Why the cue needs an offset, and why that is not a detail. Feeding a
// Turing activator straight into the reinforcement slot fails for a
// reason this model has hit before: the Schnakenberg steady state is
// u* = a + b, so u carries a large constant term. With u* = 1 the wall is
// uniformly stiffened (1 + beta*u = 5x) and uniformly growth-arrested
// (1/(1 + gamma*u) = 0.1) before the pattern contributes anything, and the
// tissue simply stops -- measured, circularity moved 0.765 to 0.737 in
// 48 h. The same mistake as driving ActiveWallExpansion with a force of
// non-zero mean, which inflated instead of shaping. The pattern's
// *variation* is the signal, so the cue subtracts the field's own mean:
//
//     S = clamp((u - mean(u)) / width, 0, 1)
//
// Subtracting the running mean rather than the analytic u* also keeps the
// cue honest if the pattern's mean drifts, which it does once the domain
// grows.
class CMTPatternRecruitment : public Reaction {
public:
  CMTPatternRecruitment(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3)
      throw std::runtime_error(
          "CMT::PatternRecruitment: uses three parameters "
          "(k_on, k_off, width).");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "CMT::PatternRecruitment: level 0 = wall pattern index, "
          "level 1 = wall reinforcement index.");
    configure("CMT::PatternRecruitment", p, i, 3, {1, 1},
              {"k_on", "k_off", "width"});
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t uIndex = variableIndex(0, 0);
    const size_t mIndex = variableIndex(1, 0);
    const double kOn = parameter(0), kOff = parameter(1);
    const double width = parameter(2);
    if (width <= 0.0)
      throw std::runtime_error(
          "CMT::PatternRecruitment: width must be positive.");
    const size_t n = T.numWall();
    if (!n)
      return;
    double sum = 0.0;
    for (size_t w = 0; w < n; ++w)
      sum += wallData[w][uIndex];
    const double mean = sum / double(n);

    parallelFor(n, [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        double S = (wallData[w][uIndex] - mean) / width;
        S = S < 0.0 ? 0.0 : (S > 1.0 ? 1.0 : S);
        const double m = wallData[w][mIndex];
        wallDerivs[w][mIndex] += kOn * S * (1.0 - m) - kOff * m;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(CMTPatternRecruitment, "CMT::PatternRecruitment")

// ---------------------------------------------------------------------------
// CMT::WidthRecruitment
//
// Recruits reinforcement where the cell is LOCALLY WIDE, with the same
// first-order kinetics as the other cues so the three are interchangeable.
//
// Why width. Figure 52 of stage3 showed the model nucleating too few lobes
// and then deepening them: real cells climb 2,3,4,5,7 lobes as circularity
// falls while every variant here saturates at four, and the model's
// perimeter per lobe *rises* from 1.56 to 2.85 while real cells refine from
// 2.50 to 1.39. A cue with no length scale of its own cannot do otherwise --
// curvature feedback amplifies whatever the initial condition seeded, and
// growing the cell just stretches it. Replacing curvature with a Turing
// pattern supplies a length scale and still does not help, because that
// length is fixed by diffusion constants and does not grow with the cell.
//
// The cell's own width does grow with the cell. Sapala et al. 2018 (eLife
// 7:e32794) measure that turgor stress on the outer periclinal wall scales
// with the largest inscribed circle, and argue lobing is how a cell keeps
// that circle small as it enlarges. That is a mechanism whose length scale
// is the cell itself: as a region widens, its periclinal stress rises, and
// splitting it in two is what relieves it. It is also the one piece of
// physics a strictly 2D model of the periclinal plane leaves out, because
// the stress in question is carried by the face this model does not have.
//
// The local width at a boundary vertex is the radius of the largest circle
// lying inside the cell and tangent to the boundary there. For a circle of
// radius r centred at v + r*n (n the inward normal), staying inside requires
// |q - v - r n|^2 >= r^2 for every other boundary point q, i.e.
//
//     r <= |q-v|^2 / (2 (q-v).n)     for every q with (q-v).n > 0
//
// so r is the minimum of that over the cell's own boundary. Recruitment is
// then a Hill function of r, which makes "wide" mean wide relative to
// w_half rather than in absolute microns.
//
// Cost: quadratic in the vertices of one cell, which is nothing at ~50 per
// cell, but it is quadratic per derivative evaluation if computed there, so
// it is refreshed on an interval like the self-avoidance neighbour list.
class CMTWidthRecruitment : public Reaction {
public:
  CMTWidthRecruitment(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 5)
      throw std::runtime_error(
          "CMT::WidthRecruitment: uses five parameters "
          "(k_on, k_off, w_half, n_hill, refresh_interval).");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "CMT::WidthRecruitment: level 0 = wall reinforcement index.");
    if (p[2] <= 0.0)
      throw std::runtime_error(
          "CMT::WidthRecruitment: w_half (parameter 3) must be positive.");
    configure("CMT::WidthRecruitment", p, i, 5, {1},
              {"k_on", "k_off", "w_half", "n_hill", "refresh_interval"});
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    rebuild(T, vertexData);
    elapsed_ = 0.0;
  }

  void update(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              double h) override {
    elapsed_ += h;
    if (elapsed_ < parameter(4))
      return;
    elapsed_ = 0.0;
    rebuild(T, vertexData);
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &, Matrix &,
              Matrix &wallDerivs, Matrix &) override {
    const size_t mI = variableIndex(0, 0);
    const double kOn = parameter(0), kOff = parameter(1);
    const double wHalf = parameter(2), nHill = parameter(3);
    if (width_.size() != T.numWall())
      return;                       // between a topology change and a refresh
    parallelFor(T.numWall(), [&](size_t begin, size_t end) {
      for (size_t w = begin; w < end; ++w) {
        const double r = width_[w];
        double S = 0.0;
        if (r > 0.0) {
          const double x = std::pow(r / wHalf, nHill);
          S = x / (1.0 + x);
        }
        const double m = wallData[w][mI];
        wallDerivs[w][mI] += kOn * S * (1.0 - m) - kOff * m;
      }
    });
  }

private:
  std::vector<double> width_;   // per wall, the local inscribed radius
  double elapsed_ = 0.0;

  void rebuild(Tissue &T, Matrix &vertexData) {
    width_.assign(T.numWall(), 0.0);
    if (vertexData.cols() != 2)
      return;
    std::vector<double> acc(T.numWall(), 0.0);
    std::vector<int> cnt(T.numWall(), 0);
    for (size_t c = 0; c < T.numCell(); ++c) {
      const auto &cell = T.cell(c);
      const size_t n = cell.vertices.size();
      if (n < 4)
        continue;
      const Vec3 mid = T.cellPosition(c, vertexData);
      for (size_t k = 0; k < n; ++k) {
        const size_t v = cell.vertices[k];
        const double vx = vertexData[v][0], vy = vertexData[v][1];
        // Inward normal: perpendicular to the local chord, sense fixed by
        // the cell centre rather than by winding, which a growing mesh does
        // not reliably preserve.
        const size_t vp = cell.vertices[(k + n - 1) % n];
        const size_t vn = cell.vertices[(k + 1) % n];
        double tx = vertexData[vn][0] - vertexData[vp][0];
        double ty = vertexData[vn][1] - vertexData[vp][1];
        const double tl = std::sqrt(tx * tx + ty * ty);
        if (tl <= 0.0)
          continue;
        tx /= tl; ty /= tl;
        double nx = -ty, ny = tx;
        if (nx * (mid[0] - vx) + ny * (mid[1] - vy) < 0.0) {
          nx = -nx; ny = -ny;
        }
        double r = std::numeric_limits<double>::max();
        for (size_t j = 0; j < n; ++j) {
          if (j == k)
            continue;
          const size_t q = cell.vertices[j];
          const double dx = vertexData[q][0] - vx;
          const double dy = vertexData[q][1] - vy;
          const double proj = dx * nx + dy * ny;
          if (proj <= 1e-12)
            continue;               // behind the tangent plane, no constraint
          const double cand = (dx * dx + dy * dy) / (2.0 * proj);
          if (cand < r)
            r = cand;
        }
        if (r == std::numeric_limits<double>::max())
          continue;
        // Carried on the two walls meeting at this vertex, since the
        // reinforcement variable lives on walls.
        for (size_t wI : T.vertex(v).walls) {
          if (wI < acc.size()) {
            acc[wI] += r;
            cnt[wI] += 1;
          }
        }
      }
    }
    for (size_t w = 0; w < width_.size(); ++w)
      width_[w] = cnt[w] ? acc[w] / cnt[w] : 0.0;
  }
};
TISSUE_REGISTER_REACTION(CMTWidthRecruitment, "CMT::WidthRecruitment")

// ---------------------------------------------------------------------------
// WallMechanics::PericlinalFoundation
//
// The restoring force the periclinal wall exerts on the anticlinal wall, and
// the reason a strictly 2D pavement model cannot select a lobe wavelength.
//
// Take a sinusoidal wrinkle of amplitude a and wavenumber q on a cell's
// outline. Its area change is zero to first order -- what goes in at one
// phase comes out at the next -- so neither turgor nor the area-elastic term
// resists it at all, at any wavelength. Bending resists it as a^2 q^4 and
// wall tension as a^2 q^2. Both of those penalise SHORT wavelengths and
// nothing in the model penalises long ones, so the longest wavelength the
// cell can hold always wins and the lobe count sits at the smallest number
// the geometry allows. Measured, it sits at four however the parameters,
// the cue, the mesh topology or the pattern boundary conditions are changed,
// against six or seven in tracked real cells.
//
// Selecting a wavelength needs a term in u^2 -- a force proportional to
// displacement rather than to its derivatives. Physically that is the
// periclinal wall: the anticlinal wall is clamped along the top and bottom
// faces, so moving it inward at one point has to shear the sheet attached
// to it, which costs energy locally and in proportion to how far it moved.
// With such a foundation of stiffness k the energy is
//
//     E ~ integral [ B kappa^2 + sigma u'^2 + k u^2 ]
//
// and minimising over q gives lambda = 2 pi (B/k)^(1/4), a wavelength set by
// the wall's own mechanics. Lobe number then grows with perimeter, which is
// what real cells do: their perimeter per lobe holds near 1.4 sqrt(A0) while
// perimeter doubles, where this model's rises from 1.56 to 2.85.
//
// The reference the displacement is measured from is the attachment the
// periclinal wall currently rests at. That cannot be the initial outline --
// the cell grows, and growth is not a wrinkle -- so it relaxes toward the
// current position with time constant tau. Deformations faster than tau are
// resisted; the slow expansion of the whole cell is not. tau is therefore a
// remodelling time for the periclinal attachment, and should sit between the
// mechanical relaxation time and the growth time.
class WallMechanicsPericlinalFoundation : public Reaction {
public:
  WallMechanicsPericlinalFoundation(const ParameterList &p,
                                    const IndexLevels &i) {
    if (p.size() != 2)
      throw std::runtime_error(
          "WallMechanics::PericlinalFoundation: uses two parameters "
          "(k_foundation, tau_remodel).");
    if (p[1] <= 0.0)
      throw std::runtime_error(
          "WallMechanics::PericlinalFoundation: tau_remodel must be "
          "positive.");
    configure("WallMechanics::PericlinalFoundation", p, i, 2, {},
              {"k_foundation", "tau_remodel"});
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    ref_.assign(vertexData.rows(),
                std::array<double, 2>{0.0, 0.0});
    for (size_t v = 0; v < vertexData.rows(); ++v)
      ref_[v] = {vertexData[v][0], vertexData[v][1]};
  }

  // The reference remodels between steps rather than inside derivs, so it is
  // a property of the trajectory and not of whichever Runge-Kutta stage
  // happens to evaluate it.
  void update(Tissue &, Matrix &, Matrix &, Matrix &vertexData,
              double h) override {
    if (ref_.size() != vertexData.rows()) {
      initiateFrom(vertexData);
      return;
    }
    const double f = h / parameter(1);
    const double w = f > 1.0 ? 1.0 : f;
    for (size_t v = 0; v < vertexData.rows(); ++v) {
      ref_[v][0] += w * (vertexData[v][0] - ref_[v][0]);
      ref_[v][1] += w * (vertexData[v][1] - ref_[v][1]);
    }
  }

  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "WallMechanics::PericlinalFoundation requires a 2D tissue.");
    if (ref_.size() != vertexData.rows())
      return;
    const double k = parameter(0);
    if (k == 0.0)
      return;
    parallelFor(vertexData.rows(), [&](size_t begin, size_t end) {
      for (size_t v = begin; v < end; ++v) {
        vertexDerivs[v][0] -= k * (vertexData[v][0] - ref_[v][0]);
        vertexDerivs[v][1] -= k * (vertexData[v][1] - ref_[v][1]);
      }
    });
  }

private:
  std::vector<std::array<double, 2>> ref_;
  void initiateFrom(Matrix &vertexData) {
    ref_.resize(vertexData.rows());
    for (size_t v = 0; v < vertexData.rows(); ++v)
      ref_[v] = {vertexData[v][0], vertexData[v][1]};
  }
};
TISSUE_REGISTER_REACTION(WallMechanicsPericlinalFoundation,
                         "WallMechanics::PericlinalFoundation")


// ---------------------------------------------------------------------------
// WallMechanics::SelfAvoidance
//
// Steric exclusion between wall segments. A vertex model has no notion that
// matter cannot occupy the same place twice: nothing stops one stretch of
// wall passing through another, and in this model the deepest necks do
// exactly that. Measured before this reaction existed: 1.7-2.0 crossing pairs
// per cell out of ~107 outline segments, appearing from frame 3 of 24 onward
// in every lobing run. The crossings are local and small (~0.1 um overlaps),
// but they sit precisely at the necks, which is the feature the whole model
// is about, and a self-intersecting outline is not a cell.
//
// Force: any two non-adjacent segments closer than d_min repel along the line
// between their closest points,
//
//   F = k_repel (d_min - d) n,     n = (P - Q) / |P - Q|
//
// distributed to each segment's two vertices by the barycentric position of
// its closest point, so the force is momentum conserving and applies no net
// torque about the contact. It is one-sided: segments further apart than
// d_min feel nothing, so the reaction is inert until walls actually approach
// and cannot perturb the shape the rest of the model produces.
//
// d_min should be read as the closest approach two anticlinal walls can make,
// i.e. roughly the sum of their half-thicknesses plus whatever cytoplasm must
// remain between them.
//
// Cost. Brute force is O(numWall^2) -- 7 million pairs per derivative
// evaluation on a 3746-wall tissue, evaluated several times per solver step,
// which is not affordable. Candidate pairs are therefore found once per
// update() with a uniform grid over segment bounding boxes and cached; the
// derivative pass only walks that list. The grid is rebuilt every
// refresh_interval of simulated time, and the candidate margin is widened by
// the distance walls could travel in that interval, so a pair cannot approach
// contact between rebuilds without having been listed.
class WallMechanicsSelfAvoidance : public Reaction {
public:
  WallMechanicsSelfAvoidance(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3 && p.size() != 4)
      throw std::runtime_error(
          "WallMechanics::SelfAvoidance: uses three or four parameters "
          "(k_repel, d_min, refresh_interval, [margin, default 2*d_min]).");
    if (p[1] <= 0.0)
      throw std::runtime_error(
          "WallMechanics::SelfAvoidance: d_min (parameter 2) must be "
          "positive.");
    std::vector<std::string> ids{"k_repel", "d_min", "refresh_interval"};
    if (p.size() == 4)
      ids.push_back("margin");
    configure("WallMechanics::SelfAvoidance", p, i, p.size(), {},
              std::move(ids));
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    rebuild(T, vertexData);
    elapsed_ = 0.0;
  }

  void update(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              double h) override {
    elapsed_ += h;
    if (elapsed_ < parameter(2))
      return;
    elapsed_ = 0.0;
    rebuild(T, vertexData);
  }

  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "WallMechanics::SelfAvoidance requires a 2D tissue.");
    const double k = parameter(0);
    const double dMin = parameter(1);
    // Walked in parallel with a per-thread scatter buffer. This loop is the
    // dominant cost of the whole model -- 81% of derivative time on an
    // 11980-wall mesh -- and it was the only reaction here still running
    // serially while every other one used the thread pool. The pairs write
    // into shared vertex rows, so it needs the scattering variant rather
    // than a plain parallel for.
    parallelScatter1(
        pairs_.size(), vertexDerivs,
        [&](size_t begin, size_t end, Matrix &out) {
          for (size_t i = begin; i < end; ++i) {
            const auto &pr = pairs_[i];
            const Wall &wa = T.wall(pr.first);
            const Wall &wb = T.wall(pr.second);
            const size_t a1 = wa.vertex1, a2 = wa.vertex2;
            const size_t b1 = wb.vertex1, b2 = wb.vertex2;
            double sc, tc, nx, ny;
            const double d = segmentDistance(vertexData, a1, a2, b1, b2, sc,
                                             tc, nx, ny);
            if (d >= dMin || d <= 0.0)
              continue;
            const double f = k * (dMin - d);
            const double fx = f * nx, fy = f * ny;
            out[a1][0] += (1.0 - sc) * fx;
            out[a1][1] += (1.0 - sc) * fy;
            out[a2][0] += sc * fx;
            out[a2][1] += sc * fy;
            out[b1][0] -= (1.0 - tc) * fx;
            out[b1][1] -= (1.0 - tc) * fy;
            out[b2][0] -= tc * fx;
            out[b2][1] -= tc * fy;
          }
        });
  }

  // Number of candidate pairs currently cached (diagnostic).
  size_t numPairs() const { return pairs_.size(); }

private:
  // vertex -> its sister partners, built once per refresh
  std::vector<std::vector<size_t>> sisterOf_;

  void buildSisterMap(const Tissue &T) {
    sisterOf_.assign(T.numVertex(), {});
    for (size_t i = 0; i < T.numSisterVertex(); ++i) {
      const size_t a = T.sisterVertex(i, 0), b = T.sisterVertex(i, 1);
      if (a < sisterOf_.size() && b < sisterOf_.size()) {
        sisterOf_[a].push_back(b);
        sisterOf_[b].push_back(a);
      }
    }
  }

  bool isSister(size_t v, size_t w) const {
    if (v >= sisterOf_.size())
      return false;
    for (size_t x : sisterOf_[v])
      if (x == w)
        return true;
    return false;
  }

  // Any endpoint paired. Requiring both endpoints catches a wall's exact
  // sister and nothing else, which is not enough: the sister's neighbours
  // along the same face are coincident too, about one segment away, and so
  // sit well inside the search margin. They share exactly one sister
  // endpoint, so an AND test lets them through and they dominate the pair
  // list. Two walls sharing any sister endpoint are parts of the same
  // physical wall seen from the two sides, and are never a collision.
  bool isSisterWall(const Wall &wa, const Wall &wb) const {
    if (sisterOf_.empty())
      return false;
    return isSister(wa.vertex1, wb.vertex1) ||
           isSister(wa.vertex1, wb.vertex2) ||
           isSister(wa.vertex2, wb.vertex1) ||
           isSister(wa.vertex2, wb.vertex2);
  }

public:

private:
  // Closest approach between segments a1-a2 and b1-b2, with the barycentric
  // positions of the closest points and the unit vector from the b-point to
  // the a-point. Standard clamped-parameter construction.
  static double segmentDistance(const Matrix &x, size_t a1, size_t a2,
                                size_t b1, size_t b2, double &sc, double &tc,
                                double &nx, double &ny) {
    const double ux = x[a2][0] - x[a1][0], uy = x[a2][1] - x[a1][1];
    const double vx = x[b2][0] - x[b1][0], vy = x[b2][1] - x[b1][1];
    const double wx = x[a1][0] - x[b1][0], wy = x[a1][1] - x[b1][1];
    const double a = ux * ux + uy * uy;
    const double b = ux * vx + uy * vy;
    const double c = vx * vx + vy * vy;
    const double d = ux * wx + uy * wy;
    const double e = vx * wx + vy * wy;
    const double den = a * c - b * b;
    if (den > 1e-12) {
      sc = (b * e - c * d) / den;
      tc = (a * e - b * d) / den;
    } else { // near-parallel: pin one parameter and solve the other
      sc = 0.0;
      tc = c > 0.0 ? e / c : 0.0;
    }
    sc = sc < 0.0 ? 0.0 : (sc > 1.0 ? 1.0 : sc);
    tc = c > 0.0 ? (e + b * sc) / c : 0.0;
    tc = tc < 0.0 ? 0.0 : (tc > 1.0 ? 1.0 : tc);
    sc = a > 0.0 ? (b * tc - d) / a : 0.0;
    sc = sc < 0.0 ? 0.0 : (sc > 1.0 ? 1.0 : sc);
    const double px = x[a1][0] + sc * ux, py = x[a1][1] + sc * uy;
    const double qx = x[b1][0] + tc * vx, qy = x[b1][1] + tc * vy;
    double dx = px - qx, dy = py - qy;
    const double dist = std::sqrt(dx * dx + dy * dy);
    if (dist > 0.0) {
      nx = dx / dist;
      ny = dy / dist;
    } else {
      nx = 0.0;
      ny = 0.0;
    }
    return dist;
  }

  void rebuild(Tissue &T, Matrix &vertexData) {
    pairs_.clear();
    buildSisterMap(T);
    if (vertexData.cols() != 2)
      return;
    const size_t n = T.numWall();
    const double margin =
        numParameter() == 4 ? parameter(3) : 2.0 * parameter(1);
    const double reach = parameter(1) + margin;
    static const bool oneSided = [] {
      const char *e = std::getenv("TISSUE_AVOID_ONESIDED");
      return e && *e && *e != '0';
    }();

    // Uniform grid sized to the search reach, over segment bounding boxes.
    double lo[2] = {1e30, 1e30}, hi[2] = {-1e30, -1e30};
    for (size_t v = 0; v < T.numVertex(); ++v)
      for (size_t d = 0; d < 2; ++d) {
        lo[d] = std::min(lo[d], vertexData[v][d]);
        hi[d] = std::max(hi[d], vertexData[v][d]);
      }
    const double cell = reach > 0.0 ? reach : 1.0;
    const size_t nx =
        std::max<size_t>(1, static_cast<size_t>((hi[0] - lo[0]) / cell) + 1);
    const size_t ny =
        std::max<size_t>(1, static_cast<size_t>((hi[1] - lo[1]) / cell) + 1);
    std::vector<std::vector<size_t>> grid(nx * ny);

    // Both the insertion box and the query box are inflated by reach, so the
    // list actually holds every pair within about 2*reach -- six times d_min
    // at the default margin, and 313k pairs against 119k for a single
    // inflation on an 11980-wall mesh. Every one of those extras is
    // distance-tested on every derivative evaluation, and dropping them runs
    // 1.25x faster.
    //
    // TISSUE_AVOID_ONESIDED=1 inflates the query side only. Whether that is
    // safe is an empirical question and not the one it first appears to be:
    // the two agree bit-for-bit until the first neighbour-list rebuild and
    // diverge right after it, which looks like dropped contacts but is also
    // exactly what a changed pair ordering does to a model whose vertex
    // trajectories are chaotic. The test that distinguishes them is whether
    // the scored population statistics move, not whether the coordinates do.
    // Left off by default until that has been run on a mesh where
    // self-avoidance is actually load-bearing.
    auto cellsOf = [&](size_t w, double grow, size_t &i0, size_t &i1,
                       size_t &j0, size_t &j1) {
      const Wall &wall = T.wall(w);
      const double x1 = vertexData[wall.vertex1][0];
      const double y1 = vertexData[wall.vertex1][1];
      const double x2 = vertexData[wall.vertex2][0];
      const double y2 = vertexData[wall.vertex2][1];
      const double xa = std::min(x1, x2) - grow, xb = std::max(x1, x2) + grow;
      const double ya = std::min(y1, y2) - grow, yb = std::max(y1, y2) + grow;
      auto clampIdx = [](double v, size_t m) {
        long idx = static_cast<long>(v);
        if (idx < 0)
          idx = 0;
        if (idx >= static_cast<long>(m))
          idx = static_cast<long>(m) - 1;
        return static_cast<size_t>(idx);
      };
      i0 = clampIdx((xa - lo[0]) / cell, nx);
      i1 = clampIdx((xb - lo[0]) / cell, nx);
      j0 = clampIdx((ya - lo[1]) / cell, ny);
      j1 = clampIdx((yb - lo[1]) / cell, ny);
    };

    for (size_t w = 0; w < n; ++w) {
      size_t i0, i1, j0, j1;
      cellsOf(w, oneSided ? 0.0 : reach, i0, i1, j0, j1);
      for (size_t i = i0; i <= i1; ++i)
        for (size_t j = j0; j <= j1; ++j)
          grid[i * ny + j].push_back(w);
    }

    std::vector<char> seen(n, 0);
    std::vector<size_t> touched;
    for (size_t w = 0; w < n; ++w) {
      size_t i0, i1, j0, j1;
      cellsOf(w, reach, i0, i1, j0, j1);
      const Wall &wa = T.wall(w);
      touched.clear();
      for (size_t i = i0; i <= i1; ++i)
        for (size_t j = j0; j <= j1; ++j)
          for (size_t o : grid[i * ny + j]) {
            if (o <= w || seen[o])
              continue;
            seen[o] = 1;
            touched.push_back(o);
            const Wall &wb = T.wall(o);
            // Sister walls are the two faces of one physical wall in a
            // per-face mesh. They do not share vertices -- their vertices
            // are distinct and coincident -- so the shared-vertex test below
            // does not catch them, and they sit at separation zero, which is
            // maximum repulsion. Left in, self-avoidance fights the sister
            // constraint on every wall in the tissue: measured at 73% of all
            // derivative time, and physically wrong, since the two faces of
            // one wall are not two walls colliding.
            if (isSisterWall(wa, wb))
              continue;
            // Segments sharing a vertex are neighbours along a wall chain and
            // are always in contact by construction.
            if (wa.vertex1 == wb.vertex1 || wa.vertex1 == wb.vertex2 ||
                wa.vertex2 == wb.vertex1 || wa.vertex2 == wb.vertex2)
              continue;
            pairs_.emplace_back(w, o);
          }
      for (size_t o : touched)
        seen[o] = 0;
    }
  }

  std::vector<std::pair<size_t, size_t>> pairs_;
  double elapsed_ = 0.0;
};
TISSUE_REGISTER_REACTION(WallMechanicsSelfAvoidance,
                         "WallMechanics::SelfAvoidance")

} // namespace
} // namespace tissue
