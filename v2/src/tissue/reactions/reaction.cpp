#include "tissue/reactions/reaction.h"

#include <sstream>
#include <stdexcept>

namespace tissue {

std::unique_ptr<Reaction> Reaction::create(const std::string &name,
                                           const ParameterList &parameters,
                                           const IndexLevels &indices) {
  auto reaction = Factory::instance().create(name, parameters, indices);
  if (!reaction) {
    std::ostringstream msg;
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
