//
// Strain-driven irreversible wall growth with length saturation and
// concentration (hormone) gating, after the apical-hook opening model of
// Walia, Carter et al. (2024) Dev Cell 59:3245 (Eq. 6 + growth saturation)
// and the strain-based growth of Bozorg et al. (2016) Phys Biol 13:065002:
//
//   dL/dt = k_growth * (eps - eps_th)+ * L * (1 - L/L_max)+ * H(c)
//
// with eps = (d - L)/L the elastic strain of the wall, L_max an absolute
// resting-length saturation (0 disables), and H a Hill function of the mean
// concentration c of the wall's adjacent cells:
//   hill_flag=0 (repressing, e.g. auxin): H = K^n / (K^n + c^n)
//   hill_flag=1 (activating):             H = c^n / (K^n + c^n)
// An optional wall flag variable (level 2) restricts growth to a subset of
// walls (e.g. longitudinal walls only).
//
// An optional 7th parameter strain_flag (default 1) selects the drive:
//   1: strain-gated yielding as above (passive wall creep)
//   0: active expansive growth, dL/dt = k * L * (1-L/L_max)+ * H(c) —
//      hormone-gated turgor-driven elongation whose RATE is strain-
//      independent (the subepidermal longitudinal force of Walia et al.
//      2024). The threshold still acts as an on/off gate (set it slightly
//      negative to let cells push, but stop them growing into compression).
//   2: stress-driven yielding (Lockhart in stress form, as in the paper's
//      Eq. 6 where strain derives from stress via wall compliance):
//      dL/dt = k * (sigma - threshold)+ * L * (1-L/L_max)+ * H(c), with
//      sigma read from a wall variable (level 3), e.g. the tension saved by
//      WallMechanics::Spring.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class WallGrowthStrainSaturationHill : public Reaction {
public:
  WallGrowthStrainSaturationHill(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 6 && p.size() != 7)
      throw std::runtime_error(
          "WallGrowth::StrainSaturationHill: uses six or seven parameters "
          "(k_growth, strain_threshold, L_max (0=off), K_Hill, n_Hill, "
          "hill_flag (0=repressing, 1=activating), [strain_flag "
          "(1=strain-gated, 0=active constant-rate)]).");
    if (p.size() == 7 && p[6] != 0.0 && p[6] != 1.0 && p[6] != 2.0)
      throw std::runtime_error(
          "WallGrowth::StrainSaturationHill: strain_flag must be 0, 1 or 2.");
    if (p[5] != 0.0 && p[5] != 1.0)
      throw std::runtime_error(
          "WallGrowth::StrainSaturationHill: hill_flag must be 0 or 1.");
    bool okLevels = i.size() >= 2 && i.size() <= 4;
    for (auto &level : i)
      okLevels = okLevels && level.size() == 1;
    if (!okLevels)
      throw std::runtime_error(
          "WallGrowth::StrainSaturationHill: level 0 = wall length index, "
          "level 1 = cell concentration index, optional level 2 = wall "
          "grow-flag index, optional level 3 = wall stress index (mode 2).");
    if (p.size() == 7 && p[6] == 2.0 && i.size() != 4)
      throw std::runtime_error(
          "WallGrowth::StrainSaturationHill: stress mode (strain_flag=2) "
          "needs level 3 = wall stress variable index.");
    std::vector<std::string> ids{"k_growth", "strain_threshold", "L_max",
                                 "K_Hill", "n_Hill", "hill_flag"};
    if (p.size() == 7)
      ids.push_back("strain_flag");
    configure("WallGrowth::StrainSaturationHill", p, i, p.size(),
              std::vector<size_t>(i.size(), 1), std::move(ids));
  }

  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &, Matrix &wallDerivs,
              Matrix &) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const size_t concIndex = variableIndex(1, 0);
    const bool typed = numVariableIndexLevel() >= 3;
    const size_t typeIndex = typed ? variableIndex(2, 0) : 0;
    const double k = parameter(0);
    const double threshold = parameter(1);
    const double lMax = parameter(2);
    const double kHillN = std::pow(parameter(3), parameter(4));
    const double n = parameter(4);
    const bool activating = parameter(5) == 1.0;
    const double mode = numParameter() < 7 ? 1.0 : parameter(6);
    const bool strainGated = mode == 1.0;
    const bool stressGated = mode == 2.0;
    const size_t stressIndex =
        stressGated ? variableIndex(3, 0) : 0;
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w) {
        if (typed && wallData[w][typeIndex] == 0.0)
          continue;
        const double L = wallData[w][lengthIndex];
        const double d = T.wallLengthFromVertices(w, vertexData);
        const double eps = (d - L) / L;
        const double gateValue =
            stressGated ? wallData[w][stressIndex] : eps;
        if (gateValue <= threshold)
          continue; // no growth below yield / into compression
        // Mean concentration of adjacent (non-background) cells.
        const Wall &wall = T.wall(w);
        double c = 0.0;
        size_t nc = 0;
        if (!Tissue::isBackground(wall.cell1)) {
          c += cellData[wall.cell1][concIndex];
          ++nc;
        }
        if (!Tissue::isBackground(wall.cell2)) {
          c += cellData[wall.cell2][concIndex];
          ++nc;
        }
        if (nc)
          c /= static_cast<double>(nc);
        const double cN = std::pow(c, n);
        const double hill = activating ? cN / (kHillN + cN) : kHillN / (kHillN + cN);
        double sat = 1.0;
        if (lMax > 0.0) {
          sat = 1.0 - L / lMax;
          if (sat < 0.0)
            sat = 0.0;
        }
        double drive = 1.0;
        if (strainGated)
          drive = eps - threshold;
        else if (stressGated)
          drive = wallData[w][stressIndex] - threshold;
        wallDerivs[w][lengthIndex] += k * drive * L * sat * hill;
      }
    });
  }
};
TISSUE_REGISTER_REACTION(WallGrowthStrainSaturationHill,
                         "WallGrowth::StrainSaturationHill")

} // namespace
} // namespace tissue
