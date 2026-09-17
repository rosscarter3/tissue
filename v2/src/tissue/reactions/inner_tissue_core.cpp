//
// The sub-epidermal tissue as a growing load-bearing core.
//
// Two earlier attempts to load the epidermis longitudinally failed in ways
// that between them specify what is actually needed (see
// InnerTissue::AxialGrowth and the hypocotyl README):
//
//   - a prescribed tension along the axial walls raises longitudinal stress,
//     but tension along the surface of a curved shell exerts a straightening
//     moment, so the load that flips the anisotropy also opens the hook in
//     minutes;
//   - shortening the axial walls' own reference length releases stress rather
//     than building it, because a shell free to contract simply contracts.
//
// What the inner tissue really provides is a body the epidermis cannot
// contract past. So it is modelled as a core with its own reference length
// per axial wall, which
//
//   - grows *proportionally*, dLc/dt = g Lc, so every segment lengthens by the
//     same factor. Uniform proportional elongation is a similarity transform:
//     inner and outer arcs grow alike, curvature is preserved, and no
//     straightening moment is applied. This is what the constant-tension
//     version got wrong.
//   - pushes *one-sidedly*, F = K (Lc - d) only while d < Lc. The core resists
//     being compressed and cannot pull, which is what a tissue does. The force
//     is self-limiting: the shell stretches until d ~ Lc and the force falls
//     away, so it drives the shell to a definite length rather than loading it
//     without bound.
//
// Applied to the axial walls only, this loads the epidermis longitudinally
// while leaving hoop alone - the anisotropy shift the paper's
// long_internal_stress_rate produces. Because hoop stress is fixed by the
// toroid geometry at 2.53x the reference on the inner flank and only 1.74x on
// the outer, longitudinal overtakes hoop on the *outer* flank first, which is
// the measured inner/outer switch asymmetry.
//
// The reference length starts at the wall's current length, so the core is
// stress-free at t = 0 whatever geometry it is handed.
//
#include <cmath>
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class InnerTissueGrowingCore : public Reaction {
public:
  InnerTissueGrowingCore(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 2)
      throw std::runtime_error(
          "InnerTissue::GrowingCore: uses two parameters (relative elongation "
          "rate of the core, penalty stiffness).");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "InnerTissue::GrowingCore: level 0 = wall variable flagging the "
          "walls the core bears against (the axial ones).");
    configure("InnerTissue::GrowingCore", p, i, 2, {1}, {"rate", "stiffness"});
  }

  // Core reference length starts at the current geometry: stress-free at t=0.
  void initiate(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    const size_t flagIndex = variableIndex(0, 0);
    core_.assign(T.numWall(), 0.0);
    for (size_t w = 0; w < T.numWall(); ++w)
      if (wallData[w][flagIndex] != 0.0)
        core_[w] = T.wallLengthFromVertices(w, vertexData);
  }

  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    if (core_.size() != T.numWall())
      return; // topology changed under us; initiate() will reset
    const size_t flagIndex = variableIndex(0, 0);
    const double K = parameter(1);
    const size_t dim = vertexData.cols();
    parallelScatter1(
        T.numWall(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            if (wallData[w][flagIndex] == 0.0)
              continue;
            const Wall &wall = T.wall(w);
            const size_t v1 = wall.vertex1, v2 = wall.vertex2;
            double u[3] = {0, 0, 0}, d = 0.0;
            for (size_t k = 0; k < dim; ++k) {
              u[k] = vertexData[v1][k] - vertexData[v2][k];
              d += u[k] * u[k];
            }
            d = std::sqrt(d);
            if (d <= 0.0 || d >= core_[w])
              continue; // core only resists compression
            const double f = K * (core_[w] - d) / d;
            for (size_t k = 0; k < dim; ++k) {
              out[v1][k] += f * u[k];
              out[v2][k] -= f * u[k];
            }
          }
        });
  }

  void update(Tissue &T, Matrix &, Matrix &, Matrix &, double h) override {
    if (core_.size() != T.numWall())
      return;
    const double g = parameter(0);
    for (double &c : core_)
      c *= (1.0 + g * h);
  }

private:
  std::vector<double> core_; // core reference length per wall
};
TISSUE_REGISTER_REACTION(InnerTissueGrowingCore, "InnerTissue::GrowingCore")

} // namespace
} // namespace tissue
