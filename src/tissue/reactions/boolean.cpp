//
// Boolean gates and counters, ported from legacy boolean.cc.
//
// These read and write cell variables holding 0/1 flags, so a model can carry
// a small piece of discrete state alongside the continuous chemistry - "this
// cell has been through stage A and not yet stage B". They contribute nothing
// to the derivatives and do all their work in update(), between solver steps,
// which is what keeps a flag from being smeared by an adaptive step.
//
// Legacy compares against 1 and 0 with `==` on doubles throughout. That is
// exact for flags written by these same reactions, and it is reproduced
// rather than replaced by a tolerance: a model where the inputs are not
// exactly 0 or 1 behaves differently, and changing that silently would change
// results.
//
#include <stdexcept>
#include <string>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Most of these differ only in the condition they test and in whether a false
// condition clears the output. `Latch` = true means the output is only ever
// set, never cleared, so the flag records "this has happened at some point".
template <class Condition>
void gateLoop(Tissue &T, Matrix &cellData, size_t outIndex, bool clearOnFalse,
              Condition condition) {
  for (size_t n = 0; n < T.numCell(); ++n) {
    if (condition(n))
      cellData[n][outIndex] = 1.0;
    else if (clearOnFalse)
      cellData[n][outIndex] = 0.0;
  }
}

// Shared shape for the two-input gates: one parameter selecting whether a
// false condition clears the output, two inputs, one output.
//
//   Boolean::AndGate 1 2 2 1
//     gatetype (0 = clear the output when false, otherwise latch)
//     input1 input2
//     output
template <bool kNegateSecond> class TwoInputGate : public Reaction {
public:
  TwoInputGate(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
      throw std::runtime_error(
          std::string(kName) +
          ": level 0 = the two input flags; level 1 = the output flag.");
    configure(kName, p, i, 1, {2, 1}, {"gatetype"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t a = variableIndex(0, 0), b = variableIndex(0, 1);
    gateLoop(T, cellData, variableIndex(1, 0), parameter(0) == 0.0,
             [&](size_t n) {
               return cellData[n][a] == 1.0 &&
                      cellData[n][b] == (kNegateSecond ? 0.0 : 1.0);
             });
  }

private:
  static constexpr const char *kName =
      kNegateSecond ? "Boolean::AndNotGate" : "Boolean::AndGate";
};
using BooleanAndGate = TwoInputGate<false>;
using BooleanAndNotGate = TwoInputGate<true>;
TISSUE_REGISTER_REACTION(BooleanAndGate, "Boolean::AndGate", "AndGate")
TISSUE_REGISTER_REACTION(BooleanAndNotGate, "Boolean::AndNotGate", "AndNotGate")

// input1 && !input2 && !input3, clearing the output whenever false. The only
// gate here that always clears - the other three-input forms latch.
//
//   Boolean::AndSpecialGate 0 2 3 1
class BooleanAndSpecialGate : public Reaction {
public:
  BooleanAndSpecialGate(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 3 || i[1].size() != 1)
      throw std::runtime_error(
          "Boolean::AndSpecialGate: level 0 = three input flags; level 1 = "
          "the output flag.");
    configure("Boolean::AndSpecialGate", p, i, 0, {3, 1}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t a = variableIndex(0, 0), b = variableIndex(0, 1),
                 c = variableIndex(0, 2);
    gateLoop(T, cellData, variableIndex(1, 0), /*clearOnFalse=*/true,
             [&](size_t n) {
               return cellData[n][a] == 1.0 && cellData[n][b] == 0.0 &&
                      cellData[n][c] == 0.0;
             });
  }
};
TISSUE_REGISTER_REACTION(BooleanAndSpecialGate, "Boolean::AndSpecialGate",
                         "AndSpecialGate")

// input1 && input2 && !input3, latching.
//
//   Boolean::AndSpecialGate2 0 2 3 1
class BooleanAndSpecialGate2 : public Reaction {
public:
  BooleanAndSpecialGate2(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 3 || i[1].size() != 1)
      throw std::runtime_error(
          "Boolean::AndSpecialGate2: level 0 = three input flags; level 1 = "
          "the output flag.");
    configure("Boolean::AndSpecialGate2", p, i, 0, {3, 1}, {});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t a = variableIndex(0, 0), b = variableIndex(0, 1),
                 c = variableIndex(0, 2);
    gateLoop(T, cellData, variableIndex(1, 0), /*clearOnFalse=*/false,
             [&](size_t n) {
               return cellData[n][a] == 1.0 && cellData[n][b] == 1.0 &&
                      cellData[n][c] == 0.0;
             });
  }
};
TISSUE_REGISTER_REACTION(BooleanAndSpecialGate2, "Boolean::AndSpecialGate2",
                         "AndSpecialGate2")

// As AndSpecialGate2, but the first input is a concentration tested against a
// threshold rather than a flag. Latching.
//
//   Boolean::AndSpecialGate3 1 2 3 1
//     threshold
class BooleanAndSpecialGate3 : public Reaction {
public:
  BooleanAndSpecialGate3(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 3 || i[1].size() != 1)
      throw std::runtime_error(
          "Boolean::AndSpecialGate3: level 0 = a concentration and two input "
          "flags; level 1 = the output flag.");
    configure("Boolean::AndSpecialGate3", p, i, 1, {3, 1}, {"threshold"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t a = variableIndex(0, 0), b = variableIndex(0, 1),
                 c = variableIndex(0, 2);
    gateLoop(T, cellData, variableIndex(1, 0), /*clearOnFalse=*/false,
             [&](size_t n) {
               return cellData[n][a] > parameter(0) &&
                      cellData[n][b] == 1.0 && cellData[n][c] == 0.0;
             });
  }
};
TISSUE_REGISTER_REACTION(BooleanAndSpecialGate3, "Boolean::AndSpecialGate3",
                         "AndSpecialGate3")

// Two concentrations, each over its own threshold. Latching.
//
//   Boolean::AndThresholdGate 2 2 2 1
//     threshold1, threshold2
class BooleanAndThresholdsGate : public Reaction {
public:
  BooleanAndThresholdsGate(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
      throw std::runtime_error(
          "Boolean::AndThresholdGate: level 0 = the two concentrations; "
          "level 1 = the output flag.");
    configure("Boolean::AndThresholdGate", p, i, 2, {2, 1},
              {"threshold1", "threshold2"});
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t a = variableIndex(0, 0), b = variableIndex(0, 1);
    gateLoop(T, cellData, variableIndex(1, 0), /*clearOnFalse=*/false,
             [&](size_t n) {
               return cellData[n][a] > parameter(0) &&
                      cellData[n][b] > parameter(1);
             });
  }
};
// Legacy registers the class AndThresholdsGate under the name
// "Boolean::AndThresholdGate" (no s) with the alias "AndThresholdsGate"
// (with one). Both spellings are kept so either model file loads.
TISSUE_REGISTER_REACTION(BooleanAndThresholdsGate, "Boolean::AndThresholdGate",
                         "AndThresholdsGate")

// The counters. Each adds one per update() to a cell variable when its
// condition holds, so the variable ends up counting *solver steps* rather
// than time - an adaptive solver takes a different number of them for the
// same simulated duration. That is legacy's behaviour and models using these
// were tuned against it, so it is kept.
enum class CountKind { Always, Flag, And, Or, OrSpecial };

template <CountKind kKind> class BooleanCounter : public Reaction {
public:
  BooleanCounter(const ParameterList &p, const IndexLevels &i) {
    if constexpr (kKind == CountKind::Always) {
      if (i.size() != 1 || i[0].size() != 1)
        throw std::runtime_error(std::string(kName) +
                                 ": level 0 = the counter to increment.");
      configure(kName, p, i, 0, {1}, {});
    } else if constexpr (kKind == CountKind::Flag) {
      if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
        throw std::runtime_error(
            std::string(kName) +
            ": level 0 = the input flag; level 1 = the counter.");
      configure(kName, p, i, 0, {1, 1}, {});
    } else {
      if (i.size() != 2 || i[0].size() != 2 || i[1].size() != 1)
        throw std::runtime_error(
            std::string(kName) +
            ": level 0 = the two input flags; level 1 = the counter.");
      configure(kName, p, i, 0, {2, 1}, {});
    }
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
  void update(Tissue &T, Matrix &cellData, Matrix &, Matrix &,
              double) override {
    const size_t out =
        kKind == CountKind::Always ? variableIndex(0, 0) : variableIndex(1, 0);
    for (size_t n = 0; n < T.numCell(); ++n) {
      bool hit;
      if constexpr (kKind == CountKind::Always)
        hit = true;
      else if constexpr (kKind == CountKind::Flag)
        hit = cellData[n][variableIndex(0, 0)] == 1.0;
      else if constexpr (kKind == CountKind::And)
        hit = cellData[n][variableIndex(0, 0)] == 1.0 &&
              cellData[n][variableIndex(0, 1)] == 1.0;
      else if constexpr (kKind == CountKind::Or)
        hit = cellData[n][variableIndex(0, 0)] == 1.0 ||
              cellData[n][variableIndex(0, 1)] == 1.0;
      else // OrSpecial: the second input is a concentration, not a flag
        hit = cellData[n][variableIndex(0, 0)] == 1.0 ||
              cellData[n][variableIndex(0, 1)] > 0.0;
      if (hit)
        cellData[n][out] += 1.0;
    }
  }

private:
  static constexpr const char *kName =
      kKind == CountKind::Always  ? "Boolean::Count"
      : kKind == CountKind::Flag  ? "Boolean::FlagCount"
      : kKind == CountKind::And   ? "Boolean::AndGateCount"
      : kKind == CountKind::Or    ? "Boolean::OrGateCount"
                                  : "Boolean::OrSpecialGateCount";
};
using BooleanCount = BooleanCounter<CountKind::Always>;
using BooleanFlagCount = BooleanCounter<CountKind::Flag>;
using BooleanAndGateCount = BooleanCounter<CountKind::And>;
using BooleanOrGateCount = BooleanCounter<CountKind::Or>;
using BooleanOrSpecialGateCount = BooleanCounter<CountKind::OrSpecial>;
TISSUE_REGISTER_REACTION(BooleanCount, "Boolean::Count", "Count")
TISSUE_REGISTER_REACTION(BooleanFlagCount, "Boolean::FlagCount", "FlagCount")
TISSUE_REGISTER_REACTION(BooleanAndGateCount, "Boolean::AndGateCount",
                         "AndGateCount")
TISSUE_REGISTER_REACTION(BooleanOrGateCount, "Boolean::OrGateCount",
                         "OrGateCount")
TISSUE_REGISTER_REACTION(BooleanOrSpecialGateCount,
                         "Boolean::OrSpecialGateCount", "OrSpecialGateCount")

} // namespace
} // namespace tissue
