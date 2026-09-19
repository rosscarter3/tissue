//
// Threshold and flag bookkeeping, ported from legacy adhocReaction.cc.
//
// Small reactions that turn a continuous cell variable into a discrete flag,
// or move values between cell variables. Like the Boolean gates they work
// entirely in update(), between solver steps, so a flag is not smeared across
// the stages of an adaptive step.
//
#include <stdexcept>
#include <string>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

///
/// Raises a flag in every cell whose input reaches a threshold.
///
///   ThresholdSwitch 2 2 1 1
///     threshold, resetFlag (0 = also lower the flag below the threshold)
///     input_index
///     flag_index
///
class ThresholdSwitch : public Reaction {
public:
  ThresholdSwitch(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "ThresholdSwitch: level 0 = the input; level 1 = the flag to set.");
    configure("ThresholdSwitch", p, i, 2, {1, 1}, {"threshold", "resetFlag"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t in = variableIndex(0, 0), out = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      if (cellData[n][in] >= parameter(0))
        cellData[n][out] = 1.0;
      else if (parameter(1) == 0.0)
        cellData[n][out] = 0.0;
    }
  }
};
TISSUE_REGISTER_REACTION(ThresholdSwitch, "ThresholdSwitch")

///
/// Clears a cell variable where the input reaches a threshold.
///
///   ThresholdReset 2 2 1 1
///     threshold, resetFlag
///     input_index
///     target_index
///
/// Both of legacy's branches write zero, so `resetFlag` 0 clears the target
/// in every cell and any other value clears it only above the threshold. The
/// second branch looks like it was copied from ThresholdSwitch and not
/// finished - the shape suggests it was meant to write something else below
/// the threshold - but it is what the published runs did, so it is kept as
/// written rather than guessed at.
class ThresholdReset : public Reaction {
public:
  ThresholdReset(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "ThresholdReset: level 0 = the input; level 1 = the variable to "
          "clear.");
    configure("ThresholdReset", p, i, 2, {1, 1}, {"threshold", "resetFlag"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t in = variableIndex(0, 0), out = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      if (cellData[n][in] >= parameter(0))
        cellData[n][out] = 0.0;
      else if (parameter(1) == 0.0)
        cellData[n][out] = 0.0;
    }
  }
};
TISSUE_REGISTER_REACTION(ThresholdReset, "ThresholdReset")

///
/// Adds a fixed amount to a cell variable wherever a flag is raised - an
/// accumulator gated on discrete state.
///
///   FlagAddValue 1 2 1 1
///     value
///     flag_index
///     target_index
///
/// The flag is compared with `== 1` on a double, as in legacy and the Boolean
/// gates: it is meant to be read from a flag those reactions wrote.
class FlagAddValue : public Reaction {
public:
  FlagAddValue(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "FlagAddValue: level 0 = the flag; level 1 = the variable to add "
          "to.");
    configure("FlagAddValue", p, i, 1, {1, 1}, {"value"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t flag = variableIndex(0, 0), out = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n)
      if (cellData[n][flag] == 1.0)
        cellData[n][out] += parameter(0);
  }
};
TISSUE_REGISTER_REACTION(FlagAddValue, "FlagAddValue")

///
/// Copies one cell variable to another, once per solver step.
///
///   CopyVariable 0 2 1 1
///     source_index
///     destination_index
///
/// Useful for snapshotting a value another reaction is about to overwrite,
/// or for feeding a diagnostic into a rule that reads a different index.
class CopyVariable : public Reaction {
public:
  CopyVariable(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "CopyVariable: level 0 = the source; level 1 = the destination.");
    configure("CopyVariable", p, i, 0, {1, 1}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t in = variableIndex(0, 0), out = variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n)
      cellData[n][out] = cellData[n][in];
  }
};
TISSUE_REGISTER_REACTION(CopyVariable, "CopyVariable")

} // namespace
} // namespace tissue
