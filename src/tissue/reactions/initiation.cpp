//
// Initiation reactions: set variables once before the simulation.
// Ported from legacy initiation.cc.
//
#include <stdexcept>

#include "tissue/core/random.h"
#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Initializes a cell variable to uniform random values in [0, maxVal), after
// reseeding the GLOBAL ran3 generator with the given seed (legacy side
// effect, preserved).
class InitiationRandom : public Reaction {
public:
  InitiationRandom(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2 || p[0] < 0.0)
      throw std::runtime_error(
          "Initiation::Random: uses two parameters (maxVal >= 0, seed).");
    configure("Initiation::Random", p, i, 2, {1}, {"maxVal", "seed"});
  }
  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &, Matrix &,
                Matrix &, Matrix &) override {
    random::sran3(static_cast<long>(parameter(1)));
    const size_t cIndex = variableIndex(0, 0);
    for (size_t cell = 0; cell < cellData.rows(); ++cell)
      cellData[cell][cIndex] = parameter(0) * random::ran3();
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(InitiationRandom, "Initiation::Random")

} // namespace
} // namespace tissue
