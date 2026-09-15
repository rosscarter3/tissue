//
// Reaction base class: one term of the ODE right-hand side. Mirrors the
// legacy BaseReaction contract (parameters + layered variable indices,
// derivs() accumulating into cell/wall/vertex derivative matrices) so model
// files remain fully compatible, but uses a self-registering factory and the
// flat Matrix state.
//
#ifndef TISSUE2_REACTIONS_REACTION_H
#define TISSUE2_REACTIONS_REACTION_H

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/core/registry.h"

namespace tissue {

class Tissue;

class Reaction {
public:
  using ParameterList = std::vector<double>;
  using IndexLevels = std::vector<std::vector<size_t>>;

  virtual ~Reaction() = default;

  const std::string &id() const { return id_; }
  size_t numParameter() const { return parameter_.size(); }
  double parameter(size_t i) const { return parameter_[i]; }
  double &parameterRef(size_t i) { return parameter_[i]; }
  const std::string &parameterId(size_t i) const { return parameterId_[i]; }

  size_t numVariableIndexLevel() const { return variableIndex_.size(); }
  size_t numVariableIndex(size_t level) const {
    return variableIndex_[level].size();
  }
  size_t variableIndex(size_t i, size_t j) const {
    return variableIndex_[i][j];
  }

  void setId(std::string value) { id_ = std::move(value); }
  void setParameter(ParameterList value) { parameter_ = std::move(value); }
  void setParameterIds(std::vector<std::string> value) {
    parameterId_ = std::move(value);
  }
  void setVariableIndex(IndexLevels value) { variableIndex_ = std::move(value); }

  // Optional one-time setup before the simulation starts (may write state).
  virtual void initiate(Tissue &, Matrix & /*cellData*/, Matrix & /*wallData*/,
                        Matrix & /*vertexData*/, Matrix & /*cellDerivs*/,
                        Matrix & /*wallDerivs*/, Matrix & /*vertexDerivs*/) {}

  // Adds this reaction's contribution to the derivative matrices.
  virtual void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
                      Matrix &vertexData, Matrix &cellDerivs,
                      Matrix &wallDerivs, Matrix &vertexDerivs) = 0;

  // Variant also accumulating |contribution| (used by noise-aware solvers).
  // Default: fall back to derivs() without absolute-value tracking.
  virtual void derivsWithAbs(Tissue &T, Matrix &cellData, Matrix &wallData,
                             Matrix &vertexData, Matrix &cellDerivs,
                             Matrix &wallDerivs, Matrix &vertexDerivs,
                             Matrix & /*sdydtCell*/, Matrix & /*sdydtWall*/,
                             Matrix & /*sdydtVertex*/) {
    derivs(T, cellData, wallData, vertexData, cellDerivs, wallDerivs,
           vertexDerivs);
  }

  // Optional state update between solver steps (h = step just taken).
  // Cell variables this reaction drives as positional degrees of freedom
  // rather than as concentrations - a center-triangulation vertex stored in
  // the cell row, for instance, whose "derivative" is a force. Solvers that
  // separate mechanics from growth must not integrate these on the growth
  // step; the default (none) is right for every ordinary reaction.
  virtual void positionalCellVariables(std::vector<size_t> &) const {}

  // True if this reaction's vertex contribution is a prescribed *velocity*
  // rather than a force - a term that does not vanish at mechanical
  // equilibrium (imposed tissue-level dilation, a moving boundary). Explicit
  // solvers cannot tell the difference and do not need to, but a solver that
  // relaxes to force balance must not try to zero such a term: it never can,
  // and will instead drive the prescribed motion without bound. Reactions
  // that say true here must also implement velocityDerivs().
  //
  // Setting a derivative to *zero* (a clamp) is not prescribed velocity:
  // zero is consistent with force balance, so clamps need nothing here.
  virtual bool prescribesVelocity() const { return false; }

  // The prescribed velocity contribution, written the same way derivs()
  // writes vertexDerivs. Only called by force-balance solvers, which apply it
  // as motion over the growth step and hold it out of the relaxation.
  virtual void velocityDerivs(Tissue &T, Matrix &cellData, Matrix &wallData,
                              Matrix &vertexData, Matrix &vertexVel) {
    (void)T; (void)cellData; (void)wallData; (void)vertexData; (void)vertexVel;
  }

  virtual void update(Tissue &, Matrix & /*cellData*/, Matrix & /*wallData*/,
                      Matrix & /*vertexData*/, double /*h*/) {}

  // --- factory -------------------------------------------------------------

  using Factory = Registry<Reaction, const ParameterList &, const IndexLevels &>;

  // Creates a reaction by (possibly legacy alias) name; throws
  // std::runtime_error on unknown name or bad parameter/index layout.
  static std::unique_ptr<Reaction> create(const std::string &name,
                                          const ParameterList &parameters,
                                          const IndexLevels &indices);

protected:
  // One-call setup for concrete constructors: validates the parameter count
  // and index-level layout (indexCounts[i] = required count at level i; the
  // value kAnyCount accepts any number at that level), then stores id,
  // parameters, indices and parameter names. Throws std::runtime_error with a
  // model-author-friendly message on mismatch.
  static constexpr size_t kAnyCount = static_cast<size_t>(-1);
  void configure(const std::string &name, const ParameterList &parameters,
                 const IndexLevels &indices, size_t parameterCount,
                 const std::vector<size_t> &indexCounts,
                 std::vector<std::string> parameterIds);

private:
  std::string id_;
  ParameterList parameter_;
  std::vector<std::string> parameterId_;
  IndexLevels variableIndex_;
};

// Registers class Type under one or more model-file names (legacy aliases
// included). Usage (at namespace scope in the .cpp):
//   TISSUE_REGISTER_REACTION(CreationZero, "Creation::Zero", "CreationZero");
#define TISSUE_REGISTER_REACTION(Type, ...)                                    \
  namespace {                                                                  \
  const bool tissueReg_##Type = [] {                                           \
    ::tissue::Reaction::Factory::instance().add(                               \
        {__VA_ARGS__},                                                         \
        [](const ::tissue::Reaction::ParameterList &p,                         \
           const ::tissue::Reaction::IndexLevels &i)                           \
            -> std::unique_ptr<::tissue::Reaction> {                           \
          return std::make_unique<Type>(p, i);                                 \
        });                                                                    \
    return true;                                                               \
  }();                                                                         \
  }

} // namespace tissue

#endif
