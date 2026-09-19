#include "tissue/reactions/reaction.h"

#include <sstream>
#include <stdexcept>

namespace tissue {

namespace {

// Names legacy itself withdrew: its factory refuses them and points at the
// replacement. Saying "the legacy simulator may still support it" for these
// is simply wrong - a model using one is stale for both simulators, not
// waiting on this port - and that mislabelling turned up in the benchmark
// corpus, where 8 models were counted as coverage gaps on this basis.
// Text follows legacy's own guidance in baseReaction.cc.
struct Withdrawn {
  const char *name;
  const char *advice;
};
constexpr Withdrawn kWithdrawn[] = {
    {"Calculate::MaxVelocity", "use Calculate::VertexVelocity"},
    {"maxVelocity", "use Calculate::VertexVelocity"},
    {"WallGrowthExponentialTruncated", "use WallGrowth::Constant"},
    {"WallGrowthExponentialStressTruncated",
     "use WallGrowth::Stress with stretch_flag 1 and an L_th"},
    {"WallLengthGrowExperimental",
     "use WallGrowth::Stress with stretch_flag 0, or WallGrowth::Force"},
    {"WallGrowthConstantStress", "use WallGrowth::Stress"},
    {"WallGrowthConstantStressConcentrationHill",
     "use WallGrowth::StressConcentrationHill"},
    {"VertexFromHypocotylGrowth", "use Force::Axial"},
    {"EpidermalVertexForce", "use Force::EpidermalCoordinate"},
    {"VertexFromWallSpringExperimental",
     "use WallMechanics::Spring with parameter 2 set to 1.0"},
    {"VertexFromWallSpringAsymmetric",
     "the *SpringAsymmetric reactions were renamed *Spring"},
    {"VertexFromEpidermalWallSpringAsymmetric",
     "the *SpringAsymmetric reactions were renamed *Spring"},
    {"VertexFromEpidermalCellWallSpringAsymmetric",
     "the *SpringAsymmetric reactions were renamed *Spring"},
    {"VertexFromCellPressure", "use Pressure2D::EdgeForce"},
    {"VertexFromCellPressureVolumeNormalized",
     "use Pressure2D::AreaPotential with parameter 1 set to 1"},
    {"VertexFromCellPressureThresholdFromMaxPos",
     "use Pressure2D::AreaPotentialSpatialThreshold"},
    {"VertexFromCellInternalPressure",
     "use Pressure2D::AreaPotential with parameter 2 set to 1"},
    {"VertexFromPressureExperimental",
     "renamed to Pressure2D::AreaPotentialTargetArea"},
    {"CellVolumeExperimental", "renamed to TargetAreaFromPressure"},
    {"MembraneDiffusionSimple2",
     "removed upstream in 2021; it only ever worked on a single cell"},
};

} // namespace

std::unique_ptr<Reaction> Reaction::create(const std::string &name,
                                           const ParameterList &parameters,
                                           const IndexLevels &indices) {
  auto reaction = Factory::instance().create(name, parameters, indices);
  if (!reaction) {
    std::ostringstream msg;
    for (const Withdrawn &w : kWithdrawn) {
      if (name == w.name) {
        msg << "Reaction::create: '" << name
            << "' was withdrawn upstream and is refused by the legacy "
               "simulator too - "
            << w.advice << ".";
        throw std::runtime_error(msg.str());
      }
    }
    msg << "Reaction::create: unknown reaction '" << name << "'.\n"
        << "This reaction either does not exist or has not been ported to "
           "tissue v2 yet; the legacy simulator may still support it.";
    throw std::runtime_error(msg.str());
  }
  return reaction;
}

void Reaction::configure(const std::string &name,
                         const ParameterList &parameters,
                         const IndexLevels &indices, size_t parameterCount,
                         const std::vector<size_t> &indexCounts,
                         std::vector<std::string> parameterIds) {
  if (parameters.size() != parameterCount) {
    std::ostringstream msg;
    msg << name << ": expects " << parameterCount << " parameter(s), got "
        << parameters.size() << ".";
    throw std::runtime_error(msg.str());
  }
  if (indices.size() != indexCounts.size()) {
    std::ostringstream msg;
    msg << name << ": expects " << indexCounts.size()
        << " variable index level(s), got " << indices.size() << ".";
    throw std::runtime_error(msg.str());
  }
  for (size_t level = 0; level < indexCounts.size(); ++level) {
    if (indexCounts[level] != kAnyCount &&
        indices[level].size() != indexCounts[level]) {
      std::ostringstream msg;
      msg << name << ": expects " << indexCounts[level]
          << " variable index/indices at level " << level << ", got "
          << indices[level].size() << ".";
      throw std::runtime_error(msg.str());
    }
  }
  setId(name);
  setParameter(parameters);
  setVariableIndex(indices);
  setParameterIds(std::move(parameterIds));
}

} // namespace tissue
