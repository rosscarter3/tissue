//
// Axial cap force for 2D median-section models of pressurized tubular organs.
//
// A 2D longitudinal section misses the 3D pressure force on the organ's end
// caps (P * pi * r^2 pulling the tube axially): in 3D this force makes the
// longitudinal shell stress uniform around a curved tube, so a pressurized
// curved tube experiences a net straightening (Bourdon tube) moment. Without
// it, a 2D section shows an artificial closing moment from the outer-face
// tension. This reaction adds the cap force back: an outward axial force on
// the tissue's end-cap walls, distributed per wall in proportion to wall
// length: F = P_cap * wallLength along the outward cap normal.
// NOTE: this compensates the 2D artifact only for OPEN curved sections
// (C-shapes); for a closed hook (arms antiparallel) the two cap forces are
// nearly parallel and load the outer arc instead. The hoop-stress asymmetry
// that straightens a closed pressurized hook in 3D has no in-plane 2D
// equivalent (cf. Walia et al. 2024, who use the analytic toroid formula).
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class Pressure2DCapForce : public Reaction {
public:
  Pressure2DCapForce(const ParameterList &p, const IndexLevels &i) {
    // level 0: wall variable that is ZERO for transverse (cap-candidate)
    // walls and non-zero for longitudinal walls (e.g. the bendFlag chain
    // marker); caps are transverse walls with a background neighbor.
    configure("Pressure2D::CapForce", p, i, 1, {1}, {"P_cap"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &, Matrix &,
              Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error("Pressure2D::CapForce requires 2D.");
    (void)cellData;
    const size_t longFlagIndex = variableIndex(0, 0);
    const double pCap = parameter(0);
    for (size_t w = 0; w < T.numWall(); ++w) {
      const Wall &wall = T.wall(w);
      if (wallData[w][longFlagIndex] != 0.0)
        continue; // longitudinal wall, not a cap
      bool bg1 = Tissue::isBackground(wall.cell1);
      bool bg2 = Tissue::isBackground(wall.cell2);
      if (bg1 == bg2)
        continue; // interior wall (or degenerate)
      size_t inside = bg1 ? wall.cell2 : wall.cell1;
      auto a = vertexData[wall.vertex1];
      auto b = vertexData[wall.vertex2];
      // Outward normal: perpendicular to the wall, pointing away from the
      // adjacent cell's centroid.
      Vec3 center = T.cellPosition(inside, vertexData);
      double ex = b[0] - a[0], ey = b[1] - a[1];
      double nx = -ey, ny = ex;
      double mx = 0.5 * (a[0] + b[0]) - center[0];
      double my = 0.5 * (a[1] + b[1]) - center[1];
      if (nx * mx + ny * my < 0.0) {
        nx = -nx;
        ny = -ny;
      }
      double len = std::sqrt(nx * nx + ny * ny);
      if (len <= 0.0)
        continue;
      // Total force P_cap * wallLength split between the two vertices.
      double f = 0.5 * pCap * std::sqrt(ex * ex + ey * ey) / len;
      vertexDerivs[wall.vertex1][0] += f * nx;
      vertexDerivs[wall.vertex1][1] += f * ny;
      vertexDerivs[wall.vertex2][0] += f * nx;
      vertexDerivs[wall.vertex2][1] += f * ny;
    }
  }
};
TISSUE_REGISTER_REACTION(Pressure2DCapForce, "Pressure2D::CapForce")

} // namespace
} // namespace tissue
