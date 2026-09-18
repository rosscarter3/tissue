//
// Direction rules. Only the two rules the published models actually use are
// ported; every other legacy direction rule name is still refused at read
// time, so a model needing real direction machinery fails loudly rather than
// running with a direction that silently never moves.
//
#include "tissue/core/direction.h"

#include <sstream>
#include <stdexcept>

#include "tissue/core/tissue.h"

namespace tissue {

void DirectionRuleBase::configure(const std::string &name,
                                  const ParameterList &parameters,
                                  const IndexLevels &indices,
                                  size_t parameterCount,
                                  const std::vector<size_t> &indexCounts) {
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
  id_ = name;
  parameter_ = parameters;
  variableIndex_ = indices;
}

namespace {

std::string unknownMessage(const std::string &kind, const std::string &name) {
  return "Direction" + kind + "::create: unknown rule '" + name +
         "'. This rule either does not exist or has not been ported to tissue "
         "v2 yet; the legacy simulator may still support it.";
}

///
/// A direction that does not change during the simulation. The cell variables
/// holding it are written once, by whatever set them up, and left alone.
///
///   StaticDirection 0 0
///
class StaticDirection : public DirectionUpdate {
public:
  StaticDirection(const ParameterList &p, const IndexLevels &i) {
    configure("StaticDirection", p, i, 0, {});
  }
};
TISSUE_REGISTER_DIRECTION_UPDATE(StaticDirection, "StaticDirection")

///
/// Both daughters keep the mother's direction.
///
///   ParallellDirection 0 1 1
///     direction_index
///
/// Nothing to do: the direction lives in cellData, which division already
/// copies to the daughter. Legacy's body is a no-op for the same reason, and
/// its one branch that does work re-picks each cell's "directional wall",
/// which only exists when the update rule is WallDirection. That rule is not
/// ported, so the branch is unreachable here; porting WallDirection is what
/// would require writing it. (The spelling with two l's is legacy's, and is
/// what model files say.)
class ParallellDirection : public DirectionDivision {
public:
  ParallellDirection(const ParameterList &p, const IndexLevels &i) {
    configure("ParallellDirection", p, i, 0, {1});
  }
};
TISSUE_REGISTER_DIRECTION_DIVISION(ParallellDirection, "ParallellDirection")

} // namespace

std::unique_ptr<DirectionUpdate>
DirectionUpdate::create(const std::string &name,
                        const ParameterList &parameters,
                        const IndexLevels &indices) {
  auto rule = Factory::instance().create(name, parameters, indices);
  if (!rule)
    throw std::runtime_error(unknownMessage("Update", name));
  return rule;
}

std::unique_ptr<DirectionDivision>
DirectionDivision::create(const std::string &name,
                          const ParameterList &parameters,
                          const IndexLevels &indices) {
  auto rule = Factory::instance().create(name, parameters, indices);
  if (!rule)
    throw std::runtime_error(unknownMessage("Division", name));
  return rule;
}

} // namespace tissue
