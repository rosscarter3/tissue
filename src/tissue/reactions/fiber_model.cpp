//
// Mechanical-property feedback ported from legacy fiberModel.cc: reactions
// that read an anisotropy measure computed elsewhere (typically the stress or
// strain anisotropy written by one of the VertexFromTRBS* reactions) and move
// a cell's Young's moduli towards a target set by that anisotropy.
//
// Only General (also registered under the bare name "FiberModel") and
// Deposition are ported, because they are the only two legacy registers: the
// Linear, LinearEquilibrium, Hill and HillEquilibrium classes exist in
// fiberModel.cc but their entries in baseReaction.cc are commented out, so no
// model file can ever have reached them.
//
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Y_L target under the Hill feedback of Bozorg et al. (2014), eq. 8-9:
//   Y_M + 0.5 (1 + a^n / ((1-a)^n K^n + a^n)) Y_F
// The caller supplies the matrix offset because General's direct branch
// (flag 2) leaves it out; see the note there.
double hillTarget(double a, double Kh, double Nh, double youngFiber) {
  const double aN = std::pow(a, Nh);
  return 0.5 * (1.0 + aN / (std::pow(1.0 - a, Nh) * std::pow(Kh, Nh) + aN)) *
         youngFiber;
}

///
/// Moves the longitudinal Young's modulus towards a target set by an
/// anisotropy input, either gradually (an Euler step at rate k_rate) or in one
/// step. Used for the mechanical feedback in Bozorg et al. (2014).
///
///   FiberModel::General 8 3 1 1[2] 1
///     k_rate, velocity_threshold, linear-hill_flag (0=linear/1=Hill/
///     2=Hill_direct), k_hill, n_hill, Y_matrix, Y_fiber, initiate_flag
///
///     anisotropy_index
///     Young_L_index [FiberL_index]
///     velocity_index
///
/// All the work happens in update(); derivs() contributes nothing, so this
/// reaction changes the material between solver steps rather than within one.
class FiberModelGeneral : public Reaction {
public:
  FiberModelGeneral(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 3 || i[0].size() != 1 ||
        (i[1].size() != 1 && i[1].size() != 2) || i[2].size() != 1)
      throw std::runtime_error(
          "FiberModel::General: level 0 = the anisotropy to read; level 1 = "
          "the longitudinal Young's modulus to update (a second index is "
          "accepted for compatibility but, as in legacy, never written); "
          "level 2 = the velocity to compare against the threshold.");
    configure("FiberModel::General", p, i, 8,
              {1, i[1].size(), 1},
              {"k_rate", "velocity_threshold", "linear-hill_flag", "k_hill",
               "n_hill", "Y_matrix", "Y_fiber", "init_flag"});
  }

  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &, Matrix &,
                Matrix &, Matrix &) override {
    const size_t anisoIndex = variableIndex(0, 0);
    const size_t youngLIndex = variableIndex(1, 0);
    const double youngMatrix = parameter(5);
    const double youngFiber = parameter(6);

    if (parameter(7) == 1) { // isotropic: zero anisotropy, Y_L = Y_M + Y_F/2
      for (size_t n = 0; n < cellData.rows(); ++n) {
        cellData[n][anisoIndex] = 0.0;
        cellData[n][youngLIndex] = youngMatrix + 0.5 * youngFiber;
      }
    } else if (parameter(7) == 2) { // start from the anisotropy already stored
      for (size_t n = 0; n < cellData.rows(); ++n) {
        const double a = cellData[n][anisoIndex];
        if (parameter(2) == 0)
          cellData[n][youngLIndex] =
              youngMatrix + 0.5 * (1.0 + a) * youngFiber;
        else if (parameter(2) == 1 || parameter(2) == 2)
          cellData[n][youngLIndex] =
              youngMatrix + hillTarget(a, parameter(3), parameter(4),
                                       youngFiber);
      }
    }
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}

  void update(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              double h) override {
    if (parameter(0) == 0.0)
      return;

    const size_t anisoIndex = variableIndex(0, 0);
    const size_t youngLIndex = variableIndex(1, 0);
    const size_t velocityIndex = variableIndex(2, 0);
    const double kRate = parameter(0);
    const double threshold = parameter(1);
    const double Kh = parameter(3);
    const double Nh = parameter(4);
    const double youngMatrix = parameter(5);
    const double youngFiber = parameter(6);
    const int flag = static_cast<int>(parameter(2));

    for (size_t n = 0; n < cellData.rows(); ++n) {
      const double a = cellData[n][anisoIndex];
      const double youngL = cellData[n][youngLIndex];

      // The two gradual branches only fire near mechanical equilibrium (the
      // velocity gate) and stop at Y_M + Y_F, the value an anisotropy of 1
      // would give: with an *absolute* anisotropy measure the target can
      // otherwise run past the material limits the mechanics assume.
      if (flag == 0 && cellData[n][velocityIndex] < threshold &&
          youngL < youngMatrix + youngFiber) {
        const double target = youngMatrix + 0.5 * (1.0 + a) * youngFiber;
        cellData[n][youngLIndex] += kRate * h * (target - youngL);
      } else if (flag == 1 && cellData[n][velocityIndex] < threshold &&
                 youngL < youngMatrix + youngFiber) {
        const double target =
            youngMatrix + hillTarget(a, Kh, Nh, youngFiber);
        cellData[n][youngLIndex] += kRate * h * (target - youngL);
      } else if (flag == 2) {
        // Direct update, ungated: no velocity threshold and no ceiling. Two
        // legacy quirks are kept deliberately. The matrix term Y_M is *not*
        // added here, unlike every other branch and unlike the documented
        // equation, so this branch writes a fiber contribution rather than a
        // modulus; and variableIndex(1,1) ("FiberL_index") is read by legacy
        // but never written, the result going to Young_L_index like the rest.
        cellData[n][youngLIndex] = hillTarget(a, Kh, Nh, youngFiber);
      }
    }
  }
};
TISSUE_REGISTER_REACTION(FiberModelGeneral, "FiberModel::General", "FiberModel")

///
/// Redistributes a fixed total amount of fiber between cells in proportion to
/// how far each cell's maximal stress exceeds a floor, capped at a ceiling,
/// weighted by cell area.
///
///   FiberModel::Deposition 7 1 6
///     k_rate, randomness, velocity_threshold, init_flag, k_hill, n_hill,
///     Y_fiber
///
///     anisotropy, max_stress, area, Y_fiber_cell, Y_fiberL_cell, velocity
///
/// The stress window is hard-coded in legacy as [8, 14] and kept that way
/// here; it is not exposed as a parameter.
class FiberModelDeposition : public Reaction {
public:
  FiberModelDeposition(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 6)
      throw std::runtime_error(
          "FiberModel::Deposition: level 0 lists six cell variables: "
          "anisotropy, maximal stress, area, fiber, longitudinal fiber, "
          "velocity.");
    configure("FiberModel::Deposition", p, i, 7, {6},
              {"k_rate", "randomness", "velocity_threshold", "init_flag",
               "k_Hill", "n_Hill", "Y_fiber"});
  }

  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &, Matrix &,
                Matrix &, Matrix &) override {
    const size_t fiberIndex = variableIndex(0, 3);
    const size_t fiberLIndex = variableIndex(0, 4);
    const double youngFiber = parameter(6);
    const int initFlag = static_cast<int>(parameter(3));

    if (initFlag == 1) {
      for (size_t n = 0; n < cellData.rows(); ++n) {
        cellData[n][fiberIndex] = youngFiber;
        cellData[n][fiberLIndex] = 0.5 * youngFiber;
      }
    } else if (initFlag == 2 || initFlag == 3) {
      // Legacy seeds the global rand() from the clock here, so these two
      // flags give a different tissue on every run and cannot be compared
      // against legacy step for step. Kept as-is rather than given a seeded
      // generator: models using them were never reproducible either.
      std::srand(static_cast<unsigned>(std::time(nullptr)));
      const double randCoef = parameter(1);
      for (size_t n = 0; n < cellData.rows(); ++n) {
        cellData[n][fiberIndex] =
            youngFiber *
            (1.0 + randCoef * 2.0 *
                       (0.5 - static_cast<double>(std::rand()) / RAND_MAX));
        cellData[n][fiberLIndex] = 0.5 * cellData[n][fiberIndex];
      }
    }
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}

  void update(Tissue &, Matrix &cellData, Matrix &, Matrix &,
              double /*h*/) override {
    if (parameter(0) == 0.0)
      return;

    const size_t anisoIndex = variableIndex(0, 0);
    const size_t stressIndex = variableIndex(0, 1);
    const size_t areaIndex = variableIndex(0, 2);
    const size_t fiberIndex = variableIndex(0, 3);
    const size_t fiberLIndex = variableIndex(0, 4);
    const size_t velocityIndex = variableIndex(0, 5);
    const double kRate = parameter(0);
    const double vThresh = parameter(2);
    const double Kh = parameter(4);
    const double Nh = parameter(5);
    const double youngFiber = parameter(6);
    const size_t numCell = cellData.rows();

    constexpr double kStressMax = 14.0;
    constexpr double kStressMin = 8.0;

    // Deposition is a tissue-wide redistribution, so it waits for *every*
    // cell to be near equilibrium rather than gating cell by cell.
    for (size_t n = 0; n < numCell; ++n)
      if (cellData[n][velocityIndex] > vThresh)
        return;

    // Stress above the floor, area-weighted, is what each cell competes with.
    double totalStressArea = 0.0;
    double totalArea = 0.0;
    for (size_t n = 0; n < numCell; ++n) {
      const double s = cellData[n][stressIndex];
      if (s > kStressMin)
        totalStressArea += (std::min(s, kStressMax) - kStressMin) *
                           cellData[n][areaIndex];
      totalArea += cellData[n][areaIndex];
    }

    for (size_t n = 0; n < numCell; ++n) {
      const double s = cellData[n][stressIndex];
      // Legacy divides unguarded, so a tissue entirely below the floor gives
      // a non-finite target and poisons the fiber field; kept, but only
      // reachable when nothing is above the floor, in which case every share
      // is zero anyway.
      double targetFiber = 0.0;
      if (s > kStressMin)
        targetFiber = youngFiber * totalArea *
                      (std::min(s, kStressMax) - kStressMin) / totalStressArea;

      // k_rate here is a fraction of the gap per update, not a rate: legacy
      // does not multiply by the step h, so the approach speed depends on how
      // often update() is called rather than on simulated time.
      cellData[n][fiberIndex] += kRate * (targetFiber - cellData[n][fiberIndex]);

      if (parameter(3) == 3) {
        const double a = cellData[n][anisoIndex];
        cellData[n][fiberLIndex] +=
            kRate * (hillTarget(a, Kh, Nh, cellData[n][fiberIndex]) -
                     cellData[n][fiberLIndex]);
      } else {
        cellData[n][fiberLIndex] = 0.5 * cellData[n][fiberIndex];
      }
    }
  }
};
TISSUE_REGISTER_REACTION(FiberModelDeposition, "FiberModel::Deposition",
                         "FiberDeposition")

} // namespace
} // namespace tissue
