//
// CompartmentChange: discrete topology-changing rules (cell division and
// removal), checked once per accepted solver step. Mirrors the legacy
// BaseCompartmentChange contract: flag() decides per cell, update() performs
// the change, and numChange tells the sweep how the cell count moved
// (+1 division, -1 one removal, -2 many removals).
//
#ifndef TISSUE2_COMPARTMENT_COMPARTMENT_CHANGE_H
#define TISSUE2_COMPARTMENT_COMPARTMENT_CHANGE_H

#include <memory>
#include <string>
#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/core/registry.h"

namespace tissue {

class Tissue;

class CompartmentChange {
public:
  using ParameterList = std::vector<double>;
  using IndexLevels = std::vector<std::vector<size_t>>;

  virtual ~CompartmentChange() = default;

  const std::string &id() const { return id_; }
  int numChange() const { return numChange_; }
  size_t numParameter() const { return parameter_.size(); }
  double parameter(size_t i) const { return parameter_[i]; }
  size_t numVariableIndexLevel() const { return variableIndex_.size(); }
  size_t numVariableIndex(size_t level) const { return variableIndex_[level].size(); }
  size_t variableIndex(size_t i, size_t j) const { return variableIndex_[i][j]; }
  const std::vector<size_t> &variableIndexLevel(size_t i) const {
    return variableIndex_[i];
  }

  // Non-zero: update() should run for cell i.
  virtual int flag(Tissue &T, size_t i, Matrix &cellData, Matrix &wallData,
                   Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                   Matrix &vertexDerivs) = 0;
  virtual void update(Tissue &T, size_t i, Matrix &cellData, Matrix &wallData,
                      Matrix &vertexData, Matrix &cellDerivs,
                      Matrix &wallDerivs, Matrix &vertexDerivs) = 0;

  using Factory =
      Registry<CompartmentChange, const ParameterList &, const IndexLevels &>;

  static std::unique_ptr<CompartmentChange>
  create(const std::string &name, const ParameterList &parameters,
         const IndexLevels &indices);

protected:
  static constexpr size_t kAnyCount = static_cast<size_t>(-1);
  void configure(const std::string &name, const ParameterList &parameters,
                 const IndexLevels &indices, size_t parameterCount,
                 const std::vector<size_t> &indexCounts, int numChange);

private:
  std::string id_;
  int numChange_ = 0;
  ParameterList parameter_;
  IndexLevels variableIndex_;
};

#define TISSUE_REGISTER_COMPARTMENT_CHANGE(Type, ...)                          \
  namespace {                                                                  \
  const bool tissueRegCC_##Type = [] {                                         \
    ::tissue::CompartmentChange::Factory::instance().add(                      \
        {__VA_ARGS__},                                                         \
        [](const ::tissue::CompartmentChange::ParameterList &p,                \
           const ::tissue::CompartmentChange::IndexLevels &i)                  \
            -> std::unique_ptr<::tissue::CompartmentChange> {                  \
          return std::make_unique<Type>(p, i);                                 \
        });                                                                    \
    return true;                                                               \
  }();                                                                         \
  }

} // namespace tissue

#endif
