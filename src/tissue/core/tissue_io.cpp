//
// Tissue file I/O: model files (reaction/compartment-change/direction blocks)
// and init files (topology + state), matching the legacy grammar exactly.
// Both formats are whitespace-delimited token streams; '#' starts a comment
// to end of line (legacy myFiles::openFile behavior).
//
#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "tissue/compartment/compartment_change.h"
#include "tissue/core/tissue.h"
#include "tissue/io/file_utils.h"
#include "tissue/reactions/reaction.h"

namespace tissue {

namespace {

[[noreturn]] void parseError(const std::string &what) {
  throw std::runtime_error(what);
}

// Reads one reaction-shaped block: id, numParameter, numLevel, per-level
// counts, parameters, indices (row-major). Shared by reactions, compartment
// changes and direction rules (identical legacy grammar).
struct RuleBlock {
  std::string id;
  std::vector<double> parameters;
  std::vector<std::vector<size_t>> indices;
};

RuleBlock readRuleBlock(std::istream &in, const char *context) {
  RuleBlock block;
  size_t numParameter = 0, numLevel = 0;
  in >> block.id >> numParameter >> numLevel;
  if (!in)
    parseError(std::string(context) + ": failed to read rule header (id, "
                                      "numParameter, numVariableIndexLevel).");
  std::vector<size_t> counts(numLevel);
  for (size_t l = 0; l < numLevel; ++l)
    in >> counts[l];
  block.parameters.resize(numParameter);
  for (size_t p = 0; p < numParameter; ++p)
    in >> block.parameters[p];
  block.indices.resize(numLevel);
  for (size_t l = 0; l < numLevel; ++l) {
    block.indices[l].resize(counts[l]);
    for (size_t j = 0; j < counts[l]; ++j)
      in >> block.indices[l][j];
  }
  if (!in)
    parseError(std::string(context) + ": failed while reading rule '" +
               block.id + "' (parameters/indices).");
  return block;
}

} // namespace

// --- model file ---------------------------------------------------------------

void Tissue::readModel(const std::string &file, int verbose) {
  auto in = io::openCommentFiltered(file);
  if (!in)
    parseError("Tissue::readModel: cannot open file " + file);
  readModel(*in, verbose);
}

void Tissue::readModel(std::istream &in, int verbose) {
  size_t numReactionVal = 0, numCompartmentVal = 0, numDirectionVal = 0;
  in >> numReactionVal >> numCompartmentVal >> numDirectionVal;
  if (!in)
    parseError("Tissue::readModel: failed to read "
               "'numReaction numCompartmentChange numDirection' header.");
  if (numDirectionVal > 1)
    parseError("Tissue::readModel: numDirection must be 0 or 1.");

  reactions_.clear();
  for (size_t i = 0; i < numReactionVal; ++i) {
    RuleBlock block = readRuleBlock(in, "Tissue::readModel(reaction)");
    if (verbose)
      std::cerr << "Reaction " << block.id << " added." << std::endl;
    reactions_.push_back(Reaction::create(block.id, block.parameters, block.indices));
  }
  compartmentChanges_.clear();
  for (size_t i = 0; i < numCompartmentVal; ++i) {
    RuleBlock block = readRuleBlock(in, "Tissue::readModel(compartmentChange)");
    if (verbose)
      std::cerr << "CompartmentChange " << block.id << " added." << std::endl;
    compartmentChanges_.push_back(
        CompartmentChange::create(block.id, block.parameters, block.indices));
  }
  if (numDirectionVal == 1) {
    RuleBlock update = readRuleBlock(in, "Tissue::readModel(directionUpdate)");
    RuleBlock division = readRuleBlock(in, "Tissue::readModel(directionDivision)");
    setDirection(
        DirectionUpdate::create(update.id, update.parameters, update.indices),
        DirectionDivision::create(division.id, division.parameters,
                                  division.indices));
    if (verbose)
      std::cerr << "Direction " << update.id << " / " << division.id
                << " added." << std::endl;
  }
}

// --- init file -----------------------------------------------------------------

void Tissue::readInit(const std::string &file, int verbose) {
  auto in = io::openCommentFiltered(file);
  if (!in)
    parseError("Tissue::readInit: cannot open file " + file);
  readInit(*in, verbose);
}

namespace {
// Shared topology + vertex/wall sections of readInit and readInitCenterTri.
void readTopologyAndWallSections(std::istream &in, Tissue &T,
                                 std::vector<CellTopo> &cells,
                                 std::vector<Wall> &walls,
                                 std::vector<VertexTopo> &vertices,
                                 Matrix &vertexState, Matrix &wallState) {
  size_t numCellVal = 0, numWallVal = 0, numVertexVal = 0;
  in >> numCellVal >> numWallVal >> numVertexVal;
  if (!in || !numCellVal || !numWallVal || !numVertexVal)
    parseError("Tissue::readInit: bad 'numCell numWall numVertex' header.");
  cells.assign(numCellVal, {});
  walls.assign(numWallVal, {});
  vertices.assign(numVertexVal, {});

  // Wall connectivity: w c1 c2 v1 v2 (ints; -1 = background).
  for (size_t nW = 0; nW < numWallVal; ++nW) {
    long long w, c1, c2, v1, v2;
    in >> w >> c1 >> c2 >> v1 >> v2;
    if (!in)
      parseError("Tissue::readInit: failed reading wall connectivity row " +
                 std::to_string(nW) + ".");
    if (static_cast<size_t>(w) != nW)
      parseError("Tissue::readInit: wall rows must be in index order.");
    if (v1 < 0 || static_cast<size_t>(v1) >= numVertexVal || v2 < 0 ||
        static_cast<size_t>(v2) >= numVertexVal)
      parseError("Tissue::readInit: bad vertex index in wall row " +
                 std::to_string(nW) + ".");
    auto toCell = [&](long long c) -> size_t {
      if (c == -1)
        return kBackground;
      if (c < 0 || static_cast<size_t>(c) >= numCellVal)
        parseError("Tissue::readInit: bad cell index in wall row " +
                   std::to_string(nW) + ".");
      return static_cast<size_t>(c);
    };
    Wall &wall = walls[nW];
    wall.vertex1 = static_cast<size_t>(v1);
    wall.vertex2 = static_cast<size_t>(v2);
    wall.cell1 = toCell(c1);
    wall.cell2 = toCell(c2);
    vertices[wall.vertex1].walls.push_back(nW);
    vertices[wall.vertex2].walls.push_back(nW);
    for (size_t c : {wall.cell1, wall.cell2}) {
      if (Tissue::isBackground(c))
        continue;
      cells[c].walls.push_back(nW);
      for (size_t v : {wall.vertex1, wall.vertex2}) {
        if (!cells[c].hasVertex(v)) {
          cells[c].vertices.push_back(v);
          vertices[v].cells.push_back(c);
        }
      }
    }
  }

  // Vertex positions.
  size_t numVertexTmp = 0, dimensionVal = 0;
  in >> numVertexTmp >> dimensionVal;
  if (!in || numVertexTmp != numVertexVal)
    parseError("Tissue::readInit: vertex section header mismatch.");
  if (dimensionVal != 2 && dimensionVal != 3)
    parseError("Tissue::readInit: dimension must be 2 or 3.");
  vertexState.assign(numVertexVal, dimensionVal);
  for (size_t i = 0; i < numVertexVal; ++i)
    for (size_t d = 0; d < dimensionVal; ++d)
      in >> vertexState[i][d];
  if (!in)
    parseError("Tissue::readInit: failed reading vertex positions.");

  // Wall data: numWall numLength numVariable, then (length, vars...) rows.
  size_t numWallTmp = 0, numLength = 0, numWallVar = 0;
  in >> numWallTmp >> numLength >> numWallVar;
  if (!in || numWallTmp != numWallVal)
    parseError("Tissue::readInit: wall data section header mismatch.");
  if (numLength != 1)
    parseError("Tissue::readInit: numLength must be 1 (legacy supports only "
               "one resting length per wall).");
  wallState.assign(numWallVal, 1 + numWallVar);
  for (size_t i = 0; i < numWallVal; ++i)
    for (size_t j = 0; j < 1 + numWallVar; ++j)
      in >> wallState[i][j];
  if (!in)
    parseError("Tissue::readInit: failed reading wall data.");
  (void)T;
}
} // namespace

void Tissue::readInit(std::istream &in, int verbose) {
  readTopologyAndWallSections(in, *this, cells_, walls_, vertices_,
                              vertexState_, wallState_);
  // Cell data: numCell numVariable then rows (section absent when 0 vars).
  size_t numCellTmp = 0, numCellVar = 0;
  in >> numCellTmp >> numCellVar;
  if (!in || numCellTmp != numCell())
    parseError("Tissue::readInit: cell data section header mismatch.");
  numCellVariable_ = numCellVar;
  cellState_.assign(numCell(), numCellVar);
  if (numCellVar) {
    for (size_t i = 0; i < numCell(); ++i)
      for (size_t j = 0; j < numCellVar; ++j)
        in >> cellState_[i][j];
    if (!in)
      parseError("Tissue::readInit: failed reading cell variables.");
  }
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}

void Tissue::readInitCenterTri(const std::string &file, int verbose) {
  auto in = io::openCommentFiltered(file);
  if (!in)
    parseError("Tissue::readInitCenterTri: cannot open file " + file);
  readInitCenterTri(*in, verbose);
}

void Tissue::readInitCenterTri(std::istream &in, int verbose) {
  readTopologyAndWallSections(in, *this, cells_, walls_, vertices_,
                              vertexState_, wallState_);
  if (dimension() != 3)
    parseError("Tissue::readInitCenterTri: requires a 3D init file.");
  size_t numCellTmp = 0, numCellVar = 0;
  in >> numCellTmp >> numCellVar;
  if (!in || numCellTmp != numCell())
    parseError("Tissue::readInitCenterTri: cell data section header mismatch.");
  numCellVariable_ = numCellVar;
  cellState_.assign(numCell(), numCellVar);
  for (size_t i = 0; i < numCell(); ++i) {
    for (size_t j = 0; j < numCellVar; ++j)
      in >> cellState_[i][j];
    CellTopo &c = cells_[i];
    c.centerPosition.resize(dimension());
    for (size_t d = 0; d < dimension(); ++d)
      in >> c.centerPosition[d];
    c.edgeLength.resize(c.numVertex());
    for (size_t k = 0; k < c.numVertex(); ++k)
      in >> c.edgeLength[k];
  }
  if (!in)
    parseError("Tissue::readInitCenterTri: failed reading cell section.");
  sortCellWallAndCellVertex();
  checkConnectivity(verbose);
}

// --- init printing ----------------------------------------------------------------

void Tissue::printInitCenterTri(const Matrix &cellData, const Matrix &wallData,
                                const Matrix &vertexData,
                                std::ostream &os) const {
  const size_t dim = dimension();
  auto oldPrecision = os.precision();
  os.precision(20);
  os << numCell() << " " << numWall() << " " << numVertex() << std::endl;
  for (size_t i = 0; i < numWall(); ++i) {
    const Wall &w = walls_[i];
    auto cellOut = [&](size_t c) -> long long {
      return isBackground(c) ? -1 : static_cast<long long>(c);
    };
    os << i << " " << cellOut(w.cell1) << " " << cellOut(w.cell2) << " "
       << w.vertex1 << " " << w.vertex2 << std::endl;
  }
  os << std::endl << numVertex() << " " << dim << std::endl;
  for (size_t i = 0; i < numVertex(); ++i) {
    for (size_t d = 0; d < dim; ++d)
      os << vertexData[i][d] << " ";
    os << std::endl;
  }
  os << std::endl << numWall() << " 1 " << (wallData.cols() - 1) << std::endl;
  for (size_t i = 0; i < numWall(); ++i) {
    for (size_t j = 0; j < wallData.rowSize(i); ++j)
      os << wallData[i][j] << " ";
    os << std::endl;
  }
  os << std::endl << numCell() << " " << numCellVariable_ << std::endl;
  for (size_t i = 0; i < numCell(); ++i) {
    const size_t n = cells_[i].numVertex();
    if (cellData.rowSize(i) < numCellVariable_ + dim + n)
      parseError("Tissue::printInitCenterTri: cell " + std::to_string(i) +
                 " has no center-triangulation data in its row.");
    for (size_t j = 0; j < numCellVariable_ + dim + n; ++j)
      os << cellData[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  os.precision(oldPrecision);
}

void Tissue::printInit(const Matrix &cellData, const Matrix &wallData,
                       const Matrix &vertexData, std::ostream &os) const {
  assert(numCell() == cellData.rows());
  assert(numWall() == wallData.rows());
  assert(numVertex() == vertexData.rows());
  auto oldPrecision = os.precision();
  os.precision(20);

  os << numCell() << " " << numWall() << " " << numVertex() << std::endl;
  // Wall connectivity (background printed as -1).
  for (size_t i = 0; i < numWall(); ++i) {
    const Wall &w = walls_[i];
    auto cellOut = [&](size_t c) -> long long {
      return isBackground(c) ? -1 : static_cast<long long>(c);
    };
    os << i << " " << cellOut(w.cell1) << " " << cellOut(w.cell2) << " "
       << w.vertex1 << " " << w.vertex2 << std::endl;
  }
  os << std::endl;
  // Vertex positions.
  os << numVertex() << " " << dimension() << std::endl;
  for (size_t i = 0; i < numVertex(); ++i) {
    for (size_t d = 0; d < dimension(); ++d)
      os << vertexData[i][d] << " ";
    os << std::endl;
  }
  os << std::endl;
  // Wall data.
  os << numWall() << " 1 " << (wallData.cols() - 1) << std::endl;
  for (size_t i = 0; i < numWall(); ++i) {
    for (size_t j = 0; j < wallData.rowSize(i); ++j)
      os << wallData[i][j] << " ";
    os << std::endl;
  }
  os << std::endl;
  // Cell data, truncated to the declared cell variables (center-triangulation
  // extras appended past numCellVariable are not part of the init format).
  os << numCell() << " " << numCellVariable_ << std::endl;
  if (numCellVariable_) {
    for (size_t i = 0; i < numCell(); ++i) {
      for (size_t j = 0; j < numCellVariable_; ++j)
        os << cellData[i][j] << " ";
      os << std::endl;
    }
    os << std::endl;
  }
  os.precision(oldPrecision);
}

} // namespace tissue
