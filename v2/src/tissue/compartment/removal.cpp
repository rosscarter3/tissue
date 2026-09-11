//
// Cell removal rules. Ported from legacy compartmentRemoval.cc.
//
#include <cmath>
#include <iostream>

#include "tissue/compartment/compartment_change.h"
#include "tissue/core/tissue.h"

namespace tissue {
namespace {

// Removes cells whose center of mass lies outside R_threshold from the origin.
class RemovalOutsideRadius : public CompartmentChange {
public:
  RemovalOutsideRadius(const ParameterList &p, const IndexLevels &i) {
    configure("RemovalOutsideRadius", p, i, 1, {}, -1);
  }
  int flag(Tissue &T, size_t i, Matrix &, Matrix &, Matrix &vertexData,
           Matrix &, Matrix &, Matrix &) override {
    Vec3 center = T.cellPosition(i, vertexData);
    double r2 = 0.0;
    for (size_t d = 0; d < vertexData.cols(); ++d)
      r2 += center[d] * center[d];
    if (std::sqrt(r2) > parameter(0)) {
      std::cerr << "Cell " << i << " marked for removal (R=" << std::sqrt(r2)
                << ")" << std::endl;
      return 1;
    }
    return 0;
  }
  void update(Tissue &T, size_t i, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
              Matrix &vertexDerivs) override {
    T.removeCell(i, cellData, wallData, vertexData, cellDerivs, wallDerivs,
                 vertexDerivs);
  }
};
TISSUE_REGISTER_COMPARTMENT_CHANGE(RemovalOutsideRadius, "RemovalOutsideRadius")

} // namespace
} // namespace tissue
