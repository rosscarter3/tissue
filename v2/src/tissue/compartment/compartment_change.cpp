#include "tissue/compartment/compartment_change.h"

#include <sstream>
#include <stdexcept>

namespace tissue {

std::unique_ptr<CompartmentChange>
CompartmentChange::create(const std::string &name,
                          const ParameterList &parameters,
                          const IndexLevels &indices) {
  auto rule = Factory::instance().create(name, parameters, indices);
  if (!rule) {
    std::ostringstream msg;
    msg << "CompartmentChange::create: unknown rule '" << name << "'.\n"
        << "This division/removal rule either does not exist or has not been "
           "ported to tissue v2 yet; the legacy simulator may still support "
           "it.";
    throw std::runtime_error(msg.str());
  }
  return rule;
}

void CompartmentChange::configure(const std::string &name,
                                  const ParameterList &parameters,
                                  const IndexLevels &indices,
                                  size_t parameterCount,
                                  const std::vector<size_t> &indexCounts,
                                  int numChange) {
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
  numChange_ = numChange;
}

} // namespace tissue
