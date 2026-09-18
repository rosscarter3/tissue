//
// Direction rules: the third block of a legacy model file, after reactions and
// compartment changes. A model declares at most one direction, described by a
// pair of rules - how it changes during a simulation (DirectionUpdate) and
// what happens to it when a cell divides (DirectionDivision).
//
// The direction itself lives in cellData like any other cell variable, so a
// daughter cell inherits its mother's copy for free; these rules exist for the
// cases where that is not what should happen.
//
#ifndef TISSUE2_CORE_DIRECTION_H
#define TISSUE2_CORE_DIRECTION_H

#include <memory>
#include <string>
#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/core/registry.h"

namespace tissue {

class Tissue;

// Common parameter/index storage for both halves of a direction rule.
class DirectionRuleBase {
public:
  using ParameterList = std::vector<double>;
  using IndexLevels = std::vector<std::vector<size_t>>;

  virtual ~DirectionRuleBase() = default;

  const std::string &id() const { return id_; }
  size_t numParameter() const { return parameter_.size(); }
  double parameter(size_t i) const { return parameter_[i]; }
  size_t numVariableIndexLevel() const { return variableIndex_.size(); }
  size_t numVariableIndex(size_t level) const {
    return variableIndex_[level].size();
  }
  size_t variableIndex(size_t i, size_t j) const { return variableIndex_[i][j]; }

protected:
  // Validates the parameter count and index layout and stores them, the same
  // way Reaction::configure does. kAnyCount accepts any count at a level.
  static constexpr size_t kAnyCount = static_cast<size_t>(-1);
  void configure(const std::string &name, const ParameterList &parameters,
                 const IndexLevels &indices, size_t parameterCount,
                 const std::vector<size_t> &indexCounts);

private:
  std::string id_;
  ParameterList parameter_;
  IndexLevels variableIndex_;
};

// How the direction evolves: initiate() once before the run, update() after
// each accepted solver step.
class DirectionUpdate : public DirectionRuleBase {
public:
  virtual void initiate(Tissue &, Matrix & /*cellData*/, Matrix & /*wallData*/,
                        Matrix & /*vertexData*/, Matrix & /*cellDerivs*/,
                        Matrix & /*wallDerivs*/, Matrix & /*vertexDerivs*/) {}
  virtual void update(Tissue &, double /*h*/, Matrix & /*cellData*/,
                      Matrix & /*wallData*/, Matrix & /*vertexData*/,
                      Matrix & /*cellDerivs*/, Matrix & /*wallDerivs*/,
                      Matrix & /*vertexDerivs*/) {}

  using Factory =
      Registry<DirectionUpdate, const ParameterList &, const IndexLevels &>;
  static std::unique_ptr<DirectionUpdate> create(const std::string &name,
                                                 const ParameterList &parameters,
                                                 const IndexLevels &indices);
};

// What happens to the direction when cell `cellI` divides. Called with the
// mother's index, immediately after the split; the daughter is the last cell.
class DirectionDivision : public DirectionRuleBase {
public:
  virtual void update(Tissue &, size_t /*cellI*/, Matrix & /*cellData*/,
                      Matrix & /*wallData*/, Matrix & /*vertexData*/,
                      Matrix & /*cellDerivs*/, Matrix & /*wallDerivs*/,
                      Matrix & /*vertexDerivs*/) {}

  using Factory =
      Registry<DirectionDivision, const ParameterList &, const IndexLevels &>;
  static std::unique_ptr<DirectionDivision> create(const std::string &name,
                                                   const ParameterList &parameters,
                                                   const IndexLevels &indices);
};

} // namespace tissue

#define TISSUE_REGISTER_DIRECTION_UPDATE(Type, ...)                            \
  namespace {                                                                  \
  const bool tissueDirU_##Type = [] {                                          \
    ::tissue::DirectionUpdate::Factory::instance().add(                        \
        {__VA_ARGS__},                                                         \
        [](const ::tissue::DirectionUpdate::ParameterList &p,                  \
           const ::tissue::DirectionUpdate::IndexLevels &i)                    \
            -> std::unique_ptr<::tissue::DirectionUpdate> {                    \
          return std::make_unique<Type>(p, i);                                 \
        });                                                                    \
    return true;                                                               \
  }();                                                                         \
  }

#define TISSUE_REGISTER_DIRECTION_DIVISION(Type, ...)                          \
  namespace {                                                                  \
  const bool tissueDirD_##Type = [] {                                          \
    ::tissue::DirectionDivision::Factory::instance().add(                      \
        {__VA_ARGS__},                                                         \
        [](const ::tissue::DirectionDivision::ParameterList &p,                \
           const ::tissue::DirectionDivision::IndexLevels &i)                  \
            -> std::unique_ptr<::tissue::DirectionDivision> {                  \
          return std::make_unique<Type>(p, i);                                 \
        });                                                                    \
    return true;                                                               \
  }();                                                                         \
  }

#endif
