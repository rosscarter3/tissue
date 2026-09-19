//
// VertexFromConstStressBoundary, ported from legacy adhocReaction.cc.
//
// Applies a constant tensile stress to a rectangular template by holding each
// of its four sides flat and moving it rigidly: every vertex on a side gets
// the same velocity, the average of what the mechanics gave that side plus
// the applied stress times the template's current width or height. Used to
// stretch a patch of tissue at fixed stress rather than fixed displacement.
//
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

///
///   VertexFromConstStressBoundary 7 0
///     stress_x, stress_y, vertex_sensitivity,
///     right_x, left_x, top_y, bottom_y
///
/// The last four give where the sides start; a vertex joins a side when it
/// lies within vertex_sensitivity of that plane. 3D only - the reaction zeroes
/// the z velocity of every boundary vertex, so the template stays planar.
class VertexFromConstStressBoundary : public Reaction {
public:
  VertexFromConstStressBoundary(const ParameterList &p, const IndexLevels &i) {
    if (!i.empty())
      throw std::runtime_error(
          "VertexFromConstStressBoundary: uses no variable indices.");
    configure("VertexFromConstStressBoundary", p, i, 7, {},
              {"stress_x", "stress_y", "vertex_sensitivity", "right_x",
               "left_x", "top_y", "bottom_y"});
  }

  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "VertexFromConstStressBoundary: only implemented for three "
          "dimensions.");
    const double eps = parameter(2);
    // Each side collects the vertices near its plane, ordered by the
    // coordinate running along it. The order only decides which vertex is
    // read for the side's position, but it is legacy's, so it is kept:
    // descending, and stable, which is what its bubble sort produces.
    collect(vertexData, right_, 0, parameter(3), eps, 1);
    collect(vertexData, left_, 0, parameter(4), eps, 1);
    collect(vertexData, top_, 1, parameter(5), eps, 0);
    collect(vertexData, bottom_, 1, parameter(6), eps, 0);
    for (auto *side : {&right_, &left_, &top_, &bottom_})
      if (side->empty())
        throw std::runtime_error(
            "VertexFromConstStressBoundary: no vertex lies within "
            "vertex_sensitivity of one of the four boundary planes; check "
            "the four position parameters against the template.");
    numOldVertices_ = T.numVertex();
    // Legacy also measures the gaps between neighbouring vertices on each
    // side here. Nothing reads them - the only code that did is commented
    // out - so they are not computed.
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    // The template's current size, measured between the first vertex of
    // opposite sides. Stress times this is the force on that side.
    const double deltaY = vertexData[top_[0]][1] - vertexData[bottom_[0]][1];
    const double deltaX = vertexData[right_[0]][0] - vertexData[left_[0]][0];

    // Each side moves rigidly at the mean of its vertices' velocities plus
    // the applied load, so the side stays flat however the interior pulls.
    //
    // All four means are taken before any of them is written back. That only
    // matters when a vertex belongs to two sides along the same axis - which
    // needs a vertex_sensitivity wide enough for the sides to overlap - but
    // there it decides whether the second side averages the original
    // velocities or the ones the first side just overwrote. Legacy reads the
    // originals.
    const double fRight = sideMean(vertexDerivs, right_, 0,
                                   parameter(0) * deltaY);
    const double fLeft = sideMean(vertexDerivs, left_, 0,
                                  -parameter(0) * deltaY);
    const double fTop = sideMean(vertexDerivs, top_, 1,
                                 parameter(1) * deltaX);
    const double fBottom = sideMean(vertexDerivs, bottom_, 1,
                                    -parameter(1) * deltaX);
    writeSide(vertexDerivs, right_, 0, fRight);
    writeSide(vertexDerivs, left_, 0, fLeft);
    writeSide(vertexDerivs, top_, 1, fTop);
    writeSide(vertexDerivs, bottom_, 1, fBottom);

    // Hold the template planar.
    for (const auto *side : {&right_, &left_, &top_, &bottom_})
      for (size_t v : *side)
        vertexDerivs[v][2] = 0.0;
  }

  void update(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              double) override {
    // Only vertices created since the last call are considered, and only
    // against where the sides are now: a vertex already on a side stays on
    // it, and the sides are never re-scanned.
    const size_t numVertices = T.numVertex();
    if (numVertices <= numOldVertices_) {
      numOldVertices_ = numVertices;
      return;
    }
    const double eps = parameter(2);
    const double rx = vertexData[right_[0]][0];
    const double lx = vertexData[left_[0]][0];
    const double ty = vertexData[top_[0]][1];
    const double by = vertexData[bottom_[0]][1];
    for (size_t v = numOldVertices_; v < numVertices; ++v) {
      if (near(vertexData[v][0], rx, eps))
        right_.push_back(v);
      if (near(vertexData[v][0], lx, eps))
        left_.push_back(v);
      if (near(vertexData[v][1], ty, eps))
        top_.push_back(v);
      if (near(vertexData[v][1], by, eps))
        bottom_.push_back(v);
    }
    numOldVertices_ = numVertices;
  }

private:
  static bool near(double x, double plane, double eps) {
    return x > plane - eps && x < plane + eps;
  }

  // Vertices within eps of `plane` along `axis`, ordered by `alongAxis`
  // descending.
  static void collect(const Matrix &vertexData, std::vector<size_t> &out,
                      size_t axis, double plane, double eps,
                      size_t alongAxis) {
    out.clear();
    for (size_t v = 0; v < vertexData.rows(); ++v)
      if (near(vertexData[v][axis], plane, eps))
        out.push_back(v);
    std::stable_sort(out.begin(), out.end(), [&](size_t a, size_t b) {
      return vertexData[a][alongAxis] > vertexData[b][alongAxis];
    });
  }

  // The velocity a side moves at: the mean of its vertices' current
  // velocities along `axis`, plus the applied load spread over the side.
  static double sideMean(const Matrix &vertexDerivs,
                         const std::vector<size_t> &side, size_t axis,
                         double load) {
    double total = load;
    for (size_t v : side)
      total += vertexDerivs[v][axis];
    return total / static_cast<double>(side.size());
  }

  // Assignment, not accumulation: whatever the mechanics put there has been
  // absorbed into the mean already.
  static void writeSide(Matrix &vertexDerivs, const std::vector<size_t> &side,
                        size_t axis, double value) {
    for (size_t v : side)
      vertexDerivs[v][axis] = value;
  }

  std::vector<size_t> right_, left_, top_, bottom_;
  size_t numOldVertices_ = 0;
};
TISSUE_REGISTER_REACTION(VertexFromConstStressBoundary,
                         "VertexFromConstStressBoundary")

} // namespace
} // namespace tissue
