#ifndef TISSUE_CORE_MESHING_H
#define TISSUE_CORE_MESHING_H

#include <array>
#include <cstddef>
#include <vector>

namespace tissue {

// A triangulation of one cell's face.
//
// The first `numBoundary` points are the cell outline, in order and unmoved:
// they are shared wall vertices and nothing here may touch them. Anything
// after that is an interior point this code added.
struct PolygonMesh {
  std::vector<std::array<double, 2>> points;
  std::vector<std::array<size_t, 3>> tris; // counter-clockwise
  std::size_t numBoundary = 0;

  bool interiorOnlyAdded() const { return points.size() >= numBoundary; }
};

// Triangulate a simple polygon, adding interior points until no triangle is
// worse than `qualityBound` in radius ratio (circumradius over twice the
// inradius; 1 is equilateral).
//
// Why this exists. CenterTriangulation meshes a cell as a fan from one
// centre, which is a valid triangulation only where the cell is star-shaped
// about that centre. Once pavement cells lobe that is 4% of them, so the fan
// inverts, a TRBS membrane is handed a rest state no triangle can adopt, and
// the largest stiffness eigenvalue reaches 1e19. Measured on real lobed
// cells, against the fan: inverted triangles 377 in 30 cells against none,
// and worst radius ratio 5e7 against 4.2.
//
// Ear clipping alone -- a triangulation using only the outline's own
// vertices, which is always possible for a simple polygon and adds no
// degrees of freedom -- removes every inversion and is still not enough. It
// leaves slivers of radius ratio 1e6 and upward, and a relaxation's stable
// step tracks the worst element rather than the average one. Interior points
// are the requirement, not a refinement.
//
// `spacing` <= 0 picks 1.6x the median outline edge, which is what was
// measured; the result holds the worst element near 4 and the cotangent
// stiffness bound at 10 to 12 whatever the cell's shape, for 1.5x to 2.3x
// the vertices.
PolygonMesh triangulatePolygon(const std::vector<std::array<double, 2>> &outline,
                               double qualityBound = 5.0, double spacing = 0.0);

// Radius ratio of each triangle; infinity for an inverted or degenerate one,
// so a caller can test quality and validity with the same number.
std::vector<double> triangleQuality(const PolygonMesh &m);

} // namespace tissue

#endif
