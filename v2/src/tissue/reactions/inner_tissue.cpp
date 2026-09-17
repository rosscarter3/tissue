//
// Longitudinal loading of the epidermis by sub-epidermal growth.
//
// The epidermis is the growth-limiting layer: inner tissues are under
// compression, the epidermis under tension, and when the inner tissue
// elongates it stretches the epidermis along the organ axis. In Walia, Carter
// et al. (2024) this is the term that drives the microtubule switch - their
// hook_part.calc_long_stress() is
//
//     sigma_long = p r / (2 t)  +  long_internal_stress_rate * time
//
// against a hoop stress that is fixed by the toroid geometry,
//
//     sigma_hoop = p r / (2 t) * (2 lambda + sin phi) / (lambda + sin phi)
//
// i.e. 2.53x the reference on the inner flank and 1.74x on the outer at
// lambda = 2.88. Longitudinal stress therefore overtakes hoop on the *outer*
// flank first, which is what makes the outer microtubules switch to
// longitudinal while the inner ones stay circumferential. Without this term
// the only longitudinal load is turgor, the hoop direction always wins at
// force balance, and no switch occurs.
//
// Applying this to a cell-resolved shell is not simply a matter of adding an
// axial force. Tension along the surface walls of a *curved* shell exerts a
// straightening moment, so the load needed to flip the anisotropy opens the
// hook far faster than observed - measured at 23 degrees in 3 minutes. The
// paper's analytic parts have no such coupling: each carries a prescribed
// stress independent of the organ's shape.
//
// What the inner tissue actually does to the epidermis is stretch it. So the
// term is applied as a progressive extension of the axial walls' reference
// length rather than as an external force:
//
//     dL_rest/dt = -rate * L_rest     (axial walls only)
//
// The epidermis is thereby held progressively short of the length the growing
// core imposes, which raises longitudinal strain - and hence longitudinal
// stress - without any applied moment. Hoop walls are untouched, so the
// anisotropy shifts exactly the way the paper's does. The shape then responds
// only through the mechanics and the growth law, which is the classic
// epidermal-growth-control picture: inner tissue in compression, epidermis in
// tension, epidermal yielding setting the rate.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class InnerTissueAxialGrowth : public Reaction {
public:
  InnerTissueAxialGrowth(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1)
      throw std::runtime_error(
          "InnerTissue::AxialGrowth: uses one parameter, the rate at which "
          "the inner tissue outgrows the epidermis (per unit time).");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "InnerTissue::AxialGrowth: level 0 = wall variable flagging the "
          "axial walls (nonzero = stretched); level 1 = wall resting-length "
          "index.");
    configure("InnerTissue::AxialGrowth", p, i, 1, {1, 1}, {"rate"});
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &,
              Matrix &, Matrix &wallDerivs, Matrix &) override {
    const size_t flagIndex = variableIndex(0, 0);
    const size_t lengthIndex = variableIndex(1, 0);
    const double rate = parameter(0);
    if (rate <= 0.0)
      return;
    parallelFor(T.numWall(), [&](size_t b, size_t e) {
      for (size_t w = b; w < e; ++w)
        if (wallData[w][flagIndex] != 0.0)
          wallDerivs[w][lengthIndex] -= rate * wallData[w][lengthIndex];
    });
  }

};
TISSUE_REGISTER_REACTION(InnerTissueAxialGrowth, "InnerTissue::AxialGrowth")

} // namespace
} // namespace tissue
