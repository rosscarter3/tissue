//
// Hypocotyl3D::limitZdis, ported from legacy hypocotyl3D.cc.
//
// Written for one paper's template (Bou Daher et al. 2018, eLife) and, as
// legacy's own error message says, hard-coded to it: the cell variables that
// label a cell, the label values, and the z plane separating the two ends of
// the cylinder are all constants in the source rather than parameters. Kept
// that way, because the models that use it are that paper's.
//
#include <stdexcept>
#include <vector>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// The hypocotyl template's own conventions, from legacy.
constexpr size_t kLayerIndex = 37;     // cell layer label
constexpr size_t kBoundaryIndex = 38;  // -1 marks a boundary cell
constexpr double kBoundaryFlag = -1.0;
constexpr double kZSplit = -50.0;      // between the cylinder's two ends

///
/// Holds each end of a growing cylinder flat in z: every vertex on the end
/// takes the mean z velocity of that end, so the end translates rather than
/// deforming. Without it the boundary cells at the two ends distort and the
/// axial growth measured from them is wrong.
///
///   Hypocotyl3D::limitZdis 0 1 2
///     index index
///
/// Both indices are required by legacy and read by neither it nor this: the
/// class documentation describes a copy-from/copy-to pair belonging to some
/// other reaction. They are accepted so existing model files load.
class Hypocotyl3DLimitZdis : public Reaction {
public:
  Hypocotyl3DLimitZdis(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 2)
      throw std::runtime_error(
          "Hypocotyl3D::limitZdis: one index level with two indices (legacy "
          "requires them and neither it nor this build reads them).");
    configure("Hypocotyl3D::limitZdis", p, i, 0, {2}, {});
  }

  void initiate(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "Hypocotyl3D::limitZdis: only meaningful for a 3D template.");
    if (cellData.cols() <= kBoundaryIndex)
      throw std::runtime_error(
          "Hypocotyl3D::limitZdis: needs at least " +
          std::to_string(kBoundaryIndex + 1) +
          " cell variables; it reads the hypocotyl template's labels at "
          "fixed indices " + std::to_string(kLayerIndex) + " and " +
          std::to_string(kBoundaryIndex) + ".");

    top_.clear();
    bottom_.clear();
    for (size_t n = 0; n < T.numCell(); ++n) {
      if (cellData[n][kBoundaryIndex] != kBoundaryFlag)
        continue;
      const double layer = cellData[n][kLayerIndex];
      if (layer != -4.0 && layer != 0.0)
        continue;
      // Vertices shared between two qualifying cells are added once per
      // cell, so they count twice in the mean below. That is legacy's
      // weighting and changes the result, so it is kept rather than
      // de-duplicated.
      for (size_t v : T.cell(n).vertices) {
        if (vertexData[v][2] > kZSplit)
          top_.push_back(v);
        if (vertexData[v][2] < kZSplit)
          bottom_.push_back(v);
      }
    }
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &vertexDerivs) override {
    flatten(vertexDerivs, top_);
    flatten(vertexDerivs, bottom_);
  }

private:
  // Replaces the z velocity of every listed vertex with their mean. Legacy
  // divides by the list size unguarded, so an empty end gives NaN; this says
  // what went wrong instead, since an empty end means the labels or the z
  // plane do not match the template.
  static void flatten(Matrix &vertexDerivs, const std::vector<size_t> &side) {
    if (side.empty())
      throw std::runtime_error(
          "Hypocotyl3D::limitZdis: one end of the cylinder has no vertices; "
          "the template's cell labels or its position relative to z = " +
          std::to_string(static_cast<int>(kZSplit)) + " do not match what "
          "this reaction expects.");
    double total = 0.0;
    for (size_t v : side)
      total += vertexDerivs[v][2];
    total /= static_cast<double>(side.size());
    for (size_t v : side)
      vertexDerivs[v][2] = total;
  }

  std::vector<size_t> top_, bottom_;
};
TISSUE_REGISTER_REACTION(Hypocotyl3DLimitZdis, "Hypocotyl3D::limitZdis")

} // namespace
} // namespace tissue
