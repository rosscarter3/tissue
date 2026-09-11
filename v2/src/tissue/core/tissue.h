//
// Tissue: the vertex-based cell mesh, its state, and the model rules
// (reactions, compartment changes, directions) read from a model file.
//
// Modern re-architecture of the legacy Tissue/Cell/Wall/Vertex classes:
//  - topology is index-based (no pointer webs); the background is the
//    sentinel index kBackground, matching the legacy int(-1) -> size_t cast
//  - all numeric state lives in flat Matrix tables (cache friendly, SIMD
//    friendly, cheap to hand to solvers), with a mirror copy held here for
//    init/VTK output like the legacy Cell/Wall/Vertex "object" state
//  - file formats (init, model) and observable semantics (sorting, division,
//    removal, connectivity checking) match the legacy implementation
//
#ifndef TISSUE2_CORE_TISSUE_H
#define TISSUE2_CORE_TISSUE_H

#include <array>
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/core/types.h"

namespace tissue {

class Reaction;
class CompartmentChange;

inline constexpr size_t kBackground = static_cast<size_t>(-1);

struct Wall {
  size_t vertex1 = kBackground;
  size_t vertex2 = kBackground;
  size_t cell1 = kBackground; // kBackground = tissue exterior
  size_t cell2 = kBackground;
  // Orientation of (vertex1 -> vertex2) relative to each adjacent cell's
  // sorted cycle: +1 along, -1 against, 0 not yet sorted.
  int cellSort1 = 0;
  int cellSort2 = 0;

  bool hasVertex(size_t v) const { return vertex1 == v || vertex2 == v; }
  size_t otherVertex(size_t v) const { return vertex1 == v ? vertex2 : vertex1; }
  size_t otherCell(size_t c) const { return cell1 == c ? cell2 : cell1; }
};

struct CellTopo {
  std::vector<size_t> walls;    // cyclic after sorting: wall k joins
  std::vector<size_t> vertices; // vertex k and vertex (k+1)%n
  // Center-triangulation data read by readInitCenterTri (otherwise empty).
  std::vector<double> centerPosition;
  std::vector<double> edgeLength;

  size_t numWall() const { return walls.size(); }
  size_t numVertex() const { return vertices.size(); }
  bool hasVertex(size_t v) const {
    for (size_t x : vertices)
      if (x == v)
        return true;
    return false;
  }
};

struct VertexTopo {
  std::vector<size_t> cells; // never contains kBackground
  std::vector<size_t> walls;
};

class Tissue {
public:
  Tissue();
  ~Tissue();

  // --- topology ------------------------------------------------------------
  size_t numCell() const { return cells_.size(); }
  size_t numWall() const { return walls_.size(); }
  size_t numVertex() const { return vertices_.size(); }
  CellTopo &cell(size_t i) { return cells_[i]; }
  const CellTopo &cell(size_t i) const { return cells_[i]; }
  Wall &wall(size_t i) { return walls_[i]; }
  const Wall &wall(size_t i) const { return walls_[i]; }
  VertexTopo &vertex(size_t i) { return vertices_[i]; }
  const VertexTopo &vertex(size_t i) const { return vertices_[i]; }
  static bool isBackground(size_t cellIndex) { return cellIndex == kBackground; }

  size_t dimension() const { return vertexState_.cols(); }
  size_t numCellVariable() const { return numCellVariable_; }

  // --- state mirror (legacy Cell/Wall/Vertex variable storage) --------------
  // Row layouts: cellState rows = cell variables; wallState rows =
  // [restLength, var...]; vertexState rows = positions.
  Matrix &cellState() { return cellState_; }
  Matrix &wallState() { return wallState_; }
  Matrix &vertexState() { return vertexState_; }
  const Matrix &cellState() const { return cellState_; }
  const Matrix &wallState() const { return wallState_; }
  const Matrix &vertexState() const { return vertexState_; }

  // --- geometry (positions always read from the given vertex table) ---------
  // 2D: |shoelace| (signedFlag returns the signed value); 3D: triangle-fan
  // surface area about the cell center.
  double cellVolume(size_t c, const Matrix &vertexData, bool signedArea = false) const;
  // 2D: area-weighted polygon centroid; 3D: wall-length-weighted mean of wall
  // midpoints. Unused components are zero.
  Vec3 cellPosition(size_t c, const Matrix &vertexData) const;
  double wallLengthFromVertices(size_t w, const Matrix &vertexData) const;
  // Rejection sampling of a uniform point inside the (2D, convexly tested)
  // cell; nullopt after numberOfTries failures. Consumes 2 Rnd() per try.
  std::optional<std::array<double, 2>>
  randomPositionInCell(size_t c, const Matrix &vertexData, int numberOfTries = 10000) const;

  // --- sister vertices -------------------------------------------------------
  size_t numSisterVertex() const { return sisterVertexIndex_.size(); }
  size_t sisterVertex(size_t i, size_t k) const { return sisterVertexIndex_[i][k]; }
  void addSisterVertex(size_t v1, size_t v2) { sisterVertexIndex_.push_back({v1, v2}); }
  std::vector<std::array<size_t, 2>> &sisterVertices() { return sisterVertexIndex_; }

  // --- model rules -----------------------------------------------------------
  size_t numReaction() const { return reactions_.size(); }
  Reaction &reaction(size_t i) { return *reactions_[i]; }
  size_t numCompartmentChange() const { return compartmentChanges_.size(); }
  CompartmentChange &compartmentChange(size_t i) { return *compartmentChanges_[i]; }

  // --- file I/O ---------------------------------------------------------------
  // Both use '#'-to-end-of-line comments and whitespace-delimited tokens.
  void readModel(const std::string &file, int verbose);
  void readModel(std::istream &in, int verbose);
  void readInit(const std::string &file, int verbose);
  void readInit(std::istream &in, int verbose);
  // Init variant with center position + edge lengths appended per cell row
  // (legacy -centerTri_init; 3D only, no comment support in legacy — comments
  // are accepted here).
  void readInitCenterTri(const std::string &file, int verbose);
  void readInitCenterTri(std::istream &in, int verbose);

  // Prints state in init format (round-trippable with readInit).
  void printInit(const Matrix &cellData, const Matrix &wallData,
                 const Matrix &vertexData, std::ostream &os) const;

  // --- simulation hooks (called by solvers) -----------------------------------
  void derivs(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &wallDerivs, Matrix &vertexDerivs);
  void derivsWithAbs(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                     Matrix &cellDerivs, Matrix &wallDerivs, Matrix &vertexDerivs,
                     Matrix &sdydtCell, Matrix &sdydtWall, Matrix &sdydtVertex);
  void initiateReactions(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                         Matrix &cellDerivs, Matrix &wallDerivs,
                         Matrix &vertexDerivs);
  void updateReactions(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                       double step);
  // Direction machinery: not ported yet (models with a direction block are
  // rejected at read time), so these are structural no-ops.
  void initiateDirection(Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
                         Matrix &) {}
  void updateDirection(double, Matrix &, Matrix &, Matrix &, Matrix &,
                       Matrix &, Matrix &) {}

  void checkCompartmentChange(Matrix &cellData, Matrix &wallData,
                              Matrix &vertexData, Matrix &cellDerivs,
                              Matrix &wallDerivs, Matrix &vertexDerivs);

  // --- topology surgery --------------------------------------------------------
  // Divides cell cellI along the line between v1Pos (on cell-local wall wI)
  // and v2Pos (on cell-local wall w3I). Appends: 1 cell, 2 vertices, 3 walls.
  // Cell variables listed in volumeChangeList are split in proportion to the
  // daughters' areas. threshold >= 0 nudges cut points off near-vertices.
  void divideCell(size_t cellI, size_t wI, size_t w3I, Vec3 v1Pos, Vec3 v2Pos,
                  Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                  Matrix &cellDerivs, Matrix &wallDerivs, Matrix &vertexDerivs,
                  const std::vector<size_t> &volumeChangeList, double threshold);
  // Removes cell cellI (swap-with-last), turning shared walls into boundary
  // walls and deleting orphaned walls/vertices.
  void removeCell(size_t cellI, Matrix &cellData, Matrix &wallData,
                  Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
                  Matrix &vertexDerivs);

  // Establishes the cyclic ordering invariant (wall k joins vertex k and
  // vertex (k+1)%n) with tissue-consistent orientation flags.
  void sortCellWallAndCellVertex();
  void sortWallAndVertex(size_t cellI);

  // Validates the whole topology; on any error prints and exits (verbose!=0)
  // — run after every step like the legacy simulator.
  void checkConnectivity(int verbose) const;

private:
  void sortCellRecursive(size_t cellI, std::vector<size_t> &sortedFlag,
                         size_t &numSorted);
  // swap-with-last removal helpers; fix all index references to the moved
  // element (legacy pointer re-pointing equivalent)
  void removeCellTopo(size_t index);
  void removeWallTopo(size_t index);
  void removeVertexTopo(size_t index);

  std::vector<CellTopo> cells_;
  std::vector<Wall> walls_;
  std::vector<VertexTopo> vertices_;
  std::vector<std::array<size_t, 2>> sisterVertexIndex_;

  Matrix cellState_, wallState_, vertexState_;
  size_t numCellVariable_ = 0;

  std::vector<std::unique_ptr<Reaction>> reactions_;
  std::vector<std::unique_ptr<CompartmentChange>> compartmentChanges_;
};

} // namespace tissue

#endif
