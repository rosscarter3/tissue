#include "tissue/core/tissue.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <set>

#include "tissue/compartment/compartment_change.h"
#include "tissue/core/random.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {

namespace {
inline int mathSign(double x) { return x >= 0 ? 1 : -1; } // legacy myMath::sign
} // namespace

Tissue::Tissue() = default;
Tissue::~Tissue() = default;

// --- geometry ---------------------------------------------------------------

double Tissue::cellVolume(size_t c, const Matrix &vertexData,
                          bool signedArea) const {
  const CellTopo &cell = cells_[c];
  const size_t n = cell.numVertex();
  assert(n > 0);
  const size_t dim = vertexData.cols();
  if (dim == 2) {
    // Shoelace over the sorted vertex cycle.
    double tmpVolume = 0.0;
    for (size_t k = 0; k < n; ++k) {
      auto p = vertexData[cell.vertices[k]];
      auto pNext = vertexData[cell.vertices[(k + 1) % n]];
      tmpVolume += p[0] * pNext[1] - p[1] * pNext[0];
    }
    tmpVolume *= 0.5;
    return signedArea ? tmpVolume : std::fabs(tmpVolume);
  }
  if (dim == 3) {
    // Triangle fan about the cell center; unsigned surface area.
    Vec3 center = cellPosition(c, vertexData);
    double volume = 0.0;
    for (size_t w = 0; w < cell.numWall(); ++w) {
      auto p1 = vertexData[walls_[cell.walls[w]].vertex1];
      auto p2 = vertexData[walls_[cell.walls[w]].vertex2];
      double r1[3], r2[3];
      for (size_t d = 0; d < 3; ++d) {
        r1[d] = p1[d] - center[d];
        r2[d] = p2[d] - center[d];
      }
      double cross = 0.0;
      for (size_t d = 0; d < 3; ++d) {
        double component = r1[(d + 1) % 3] * r2[(d + 2) % 3] -
                           r1[(d + 2) % 3] * r2[(d + 1) % 3];
        cross += component * component;
      }
      volume += 0.5 * std::sqrt(cross);
    }
    return volume;
  }
  std::cerr << "Tissue::cellVolume() dimension " << dim << " not supported."
            << std::endl;
  std::exit(-1);
}

Vec3 Tissue::cellPosition(size_t c, const Matrix &vertexData) const {
  const CellTopo &cell = cells_[c];
  const size_t dim = vertexData.cols();
  Vec3 pos{0.0, 0.0, 0.0};
  if (dim == 2) {
    // Area-weighted polygon centroid (signed area cancels orientation).
    double area = cellVolume(c, vertexData, true);
    const size_t n = cell.numVertex();
    for (size_t i = 0; i < n; ++i) {
      auto p = vertexData[cell.vertices[i]];
      auto pNext = vertexData[cell.vertices[(i + 1) % n]];
      double factor = p[0] * pNext[1] - pNext[0] * p[1];
      for (size_t d = 0; d < 2; ++d)
        pos[d] += factor * (p[d] + pNext[d]);
    }
    for (size_t d = 0; d < 2; ++d)
      pos[d] /= 6.0 * area;
    return pos;
  }
  if (dim == 3) {
    // Wall-length-weighted mean of wall midpoints.
    double sumLength = 0.0;
    for (size_t w = 0; w < cell.numWall(); ++w) {
      auto p1 = vertexData[walls_[cell.walls[w]].vertex1];
      auto p2 = vertexData[walls_[cell.walls[w]].vertex2];
      double length = 0.0;
      for (size_t d = 0; d < 3; ++d)
        length += (p1[d] - p2[d]) * (p1[d] - p2[d]);
      length = std::sqrt(length);
      sumLength += length;
      for (size_t d = 0; d < 3; ++d)
        pos[d] += 0.5 * (p1[d] + p2[d]) * length;
    }
    for (size_t d = 0; d < 3; ++d)
      pos[d] /= sumLength;
    return pos;
  }
  std::cerr << "Tissue::cellPosition() dimension " << dim << " not supported."
            << std::endl;
  std::exit(-1);
}

double Tissue::wallLengthFromVertices(size_t w, const Matrix &vertexData) const {
  return distance(vertexData[walls_[w].vertex1], vertexData[walls_[w].vertex2]);
}

std::optional<std::array<double, 2>>
Tissue::randomPositionInCell(size_t c, const Matrix &vertexData,
                             int numberOfTries) const {
  if (vertexData.cols() != 2) {
    std::cerr << "Tissue::randomPositionInCell only supports two dimensions.\n";
    std::exit(EXIT_FAILURE);
  }
  const CellTopo &cell = cells_[c];
  double xmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::min(); // legacy quirk: min(), not lowest()
  double ymin = std::numeric_limits<double>::max();
  double ymax = std::numeric_limits<double>::min();
  for (size_t v : cell.vertices) {
    double x = vertexData[v][0], y = vertexData[v][1];
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
  const size_t n = cell.numVertex();
  for (int tryCounter = 0; tryCounter < numberOfTries; ++tryCounter) {
    double rx = xmin + (xmax - xmin) * random::Rnd();
    double ry = ymin + (ymax - ymin) * random::Rnd();
    int sign = 0;
    bool success = true;
    for (size_t k = 0; k < n; ++k) {
      auto p1 = vertexData[cell.vertices[k]];
      auto p2 = vertexData[cell.vertices[(k + 1) % n]];
      double vx = p2[0] - p1[0], vy = p2[1] - p1[1];
      double dx = rx - p1[0], dy = ry - p1[1];
      int s = mathSign(vx * dy - vy * dx);
      if (!sign)
        sign = s;
      else if (sign != s) {
        success = false;
        break;
      }
    }
    if (success)
      return std::array<double, 2>{rx, ry};
  }
  return std::nullopt;
}

// --- derivative dispatch -----------------------------------------------------

void Tissue::derivs(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                    Matrix &cellDerivs, Matrix &wallDerivs,
                    Matrix &vertexDerivs) {
  cellDerivs.fill(0.0);
  wallDerivs.fill(0.0);
  vertexDerivs.fill(0.0);
  for (auto &r : reactions_)
    r->derivs(*this, cellData, wallData, vertexData, cellDerivs, wallDerivs,
              vertexDerivs);
}

// Force-balance solvers need the two kinds of vertex contribution apart:
// forces (relaxable, vanish at equilibrium) and prescribed velocities
// (not relaxable, applied as motion over the growth step).
void Tissue::derivsSplit(Matrix &cellData, Matrix &wallData, Matrix &vertexData,
                         Matrix &cellDerivs, Matrix &wallDerivs,
                         Matrix &vertexDerivs, Matrix &vertexVel) {
  cellDerivs.fill(0.0);
  wallDerivs.fill(0.0);
  vertexDerivs.fill(0.0);
  vertexVel.fill(0.0);
  for (auto &r : reactions_) {
    if (r->prescribesVelocity())
      r->velocityDerivs(*this, cellData, wallData, vertexData, vertexVel);
    else
      r->derivs(*this, cellData, wallData, vertexData, cellDerivs, wallDerivs,
                vertexDerivs);
  }
}

bool Tissue::hasPrescribedVelocity() const {
  for (const auto &r : reactions_)
    if (r->prescribesVelocity())
      return true;
  return false;
}

std::vector<size_t> Tissue::positionalCellVariables() const {
  std::vector<size_t> v;
  for (const auto &r : reactions_)
    r->positionalCellVariables(v);
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  return v;
}

void Tissue::derivsWithAbs(Matrix &cellData, Matrix &wallData,
                           Matrix &vertexData, Matrix &cellDerivs,
                           Matrix &wallDerivs, Matrix &vertexDerivs,
                           Matrix &sdydtCell, Matrix &sdydtWall,
                           Matrix &sdydtVertex) {
  cellDerivs.fill(0.0);
  wallDerivs.fill(0.0);
  vertexDerivs.fill(0.0);
  sdydtCell.fill(0.0);
  sdydtWall.fill(0.0);
  sdydtVertex.fill(0.0);
  for (auto &r : reactions_)
    r->derivsWithAbs(*this, cellData, wallData, vertexData, cellDerivs,
                     wallDerivs, vertexDerivs, sdydtCell, sdydtWall,
                     sdydtVertex);
}

void Tissue::initiateReactions(Matrix &cellData, Matrix &wallData,
                               Matrix &vertexData, Matrix &cellDerivs,
                               Matrix &wallDerivs, Matrix &vertexDerivs) {
  for (auto &r : reactions_)
    r->initiate(*this, cellData, wallData, vertexData, cellDerivs, wallDerivs,
                vertexDerivs);
}

void Tissue::updateReactions(Matrix &cellData, Matrix &wallData,
                             Matrix &vertexData, double step) {
  for (auto &r : reactions_)
    r->update(*this, cellData, wallData, vertexData, step);
}

// --- compartment change sweep --------------------------------------------------

void Tissue::checkCompartmentChange(Matrix &cellData, Matrix &wallData,
                                    Matrix &vertexData, Matrix &cellDerivs,
                                    Matrix &wallDerivs, Matrix &vertexDerivs) {
  unsigned int guardCounter = 0;
  for (size_t l = 0; l < numCompartmentChange(); ++l) {
    for (size_t i = 0; i < numCell(); ++i) {
      if (++guardCounter > 1000000) {
        std::cerr << "Tissue::checkCompartmentChange() more than a million "
                     "compartment checks in one step; aborting." << std::endl;
        std::exit(EXIT_FAILURE);
      }
      CompartmentChange &rule = *compartmentChanges_[l];
      if (rule.flag(*this, i, cellData, wallData, vertexData, cellDerivs,
                    wallDerivs, vertexDerivs)) {
        rule.update(*this, i, cellData, wallData, vertexData, cellDerivs,
                    wallDerivs, vertexDerivs);
        if (rule.numChange() == 1) {
          // Division: locally re-sort mother, daughter, and their neighbors.
          // Cells with a single real neighbor are sorted in a second pass so
          // their orientation seed exists.
          std::set<size_t> sortCell;
          sortCell.insert(i);
          sortCell.insert(numCell() - 1);
          for (size_t w : cells_[i].walls)
            sortCell.insert(walls_[w].otherCell(i));
          size_t daughter = numCell() - 1;
          for (size_t w : cells_[daughter].walls)
            sortCell.insert(walls_[w].otherCell(daughter));
          sortCell.erase(kBackground);
          std::vector<size_t> oneNeighborCells;
          for (size_t cellToSort : sortCell) {
            size_t numNeighbors = 0;
            for (size_t w : cells_[cellToSort].walls)
              if (!isBackground(walls_[w].otherCell(cellToSort)))
                ++numNeighbors;
            if (numNeighbors == 1)
              oneNeighborCells.push_back(cellToSort);
            else
              sortWallAndVertex(cellToSort);
          }
          for (size_t cellToSort : oneNeighborCells)
            sortWallAndVertex(cellToSort);
        } else if (rule.numChange() == -1) {
          --i; // swapped-in cell must be re-tested
        } else if (rule.numChange() < -1) {
          i = numCell() + 1; // many removed: abort this rule's sweep
        }
      }
    }
  }
}

// --- sorting -------------------------------------------------------------------

void Tissue::sortCellWallAndCellVertex() {
  std::vector<size_t> sortedFlag(numCell(), 0);
  size_t numSorted = 0;
  size_t startSortIndex = 0;
  if (numCell()) {
    do {
      if (!sortedFlag[startSortIndex])
        sortCellRecursive(startSortIndex, sortedFlag, numSorted);
      startSortIndex++;
    } while (numSorted < numCell() && startSortIndex < numCell());
  }
  std::cerr << "Tissue::sortCellWallAndCellVertex() " << numSorted << " of "
            << numCell() << " faces (cells) sorted recursively." << std::endl;
}

void Tissue::sortCellRecursive(size_t cellI, std::vector<size_t> &sortedFlag,
                               size_t &numSorted) {
  if (sortedFlag[cellI])
    return;
  sortWallAndVertex(cellI);
  sortedFlag[cellI]++;
  numSorted++;
  for (size_t w : cells_[cellI].walls) {
    size_t next = walls_[w].otherCell(cellI);
    if (!isBackground(next))
      sortCellRecursive(next, sortedFlag, numSorted);
  }
}

void Tissue::sortWallAndVertex(size_t cellI) {
  CellTopo &cell = cells_[cellI];
  assert(cell.numWall() == cell.numVertex());
  const size_t n = cell.numWall();

  std::vector<size_t> tmpWall(n, kBackground);
  std::vector<size_t> tmpVertex(n, kBackground);

  // Choose a starting wall: the first with a set cellSort1 flag (legacy seeds
  // from previously sorted neighbors via this flag), else wall 0.
  size_t wallIndex = 0;
  bool foundCellSortFlag = false;
  while (wallIndex < n && !walls_[cell.walls[wallIndex]].cellSort1)
    ++wallIndex;
  size_t wallIndexStart;
  size_t vertexIndex;
  if (wallIndex < n) {
    foundCellSortFlag = true;
    wallIndexStart = vertexIndex = wallIndex;
  } else {
    wallIndex = wallIndexStart = vertexIndex = 0;
  }

  size_t numWallSorted = 0, numVertexSorted = 0;
  Wall &firstWall = walls_[cell.walls[wallIndex]];
  auto sortError = [&](const char *msg) {
    std::cerr << "Tissue::sortWallAndVertex() " << msg << std::endl;
    std::exit(EXIT_FAILURE);
  };

  if (!foundCellSortFlag) {
    // No previous orientation information: define it from this wall.
    tmpWall[wallIndex] = cell.walls[wallIndex];
    ++numWallSorted;
    tmpVertex[vertexIndex] = firstWall.vertex1;
    ++numVertexSorted;
    ++vertexIndex;
    tmpVertex[vertexIndex] = firstWall.vertex2;
    ++numVertexSorted;
    if (firstWall.cell1 == cellI) {
      firstWall.cellSort1 = 1;
      if (!isBackground(firstWall.cell2))
        firstWall.cellSort2 = -1;
    } else if (firstWall.cell2 == cellI) {
      firstWall.cellSort2 = 1;
      if (!isBackground(firstWall.cell1))
        firstWall.cellSort1 = -1;
    } else {
      sortError("Wall to be sorted does not belong to the cell");
    }
  } else {
    // Seed orientation from earlier sorting information.
    tmpWall[wallIndex] = cell.walls[wallIndex];
    ++numWallSorted;
    int sort;
    if (firstWall.cell1 == cellI)
      sort = firstWall.cellSort1;
    else if (firstWall.cell2 == cellI)
      sort = firstWall.cellSort2;
    else {
      sortError("Wall to be sorted does not belong to the cell");
      return;
    }
    size_t vFirst, vSecond;
    if (sort == 1) {
      vFirst = firstWall.vertex1;
      vSecond = firstWall.vertex2;
    } else if (sort == -1) {
      vFirst = firstWall.vertex2;
      vSecond = firstWall.vertex1;
    } else {
      sortError("Wall to be sorted marked as but is not sorted.");
      return;
    }
    tmpVertex[vertexIndex] = vFirst;
    ++numVertexSorted;
    vertexIndex = (vertexIndex + 1) % n;
    tmpVertex[vertexIndex] = vSecond;
    ++numVertexSorted;
  }

  // Sort the remaining walls by walking the cycle.
  while (wallIndex < wallIndexStart + n - 1) {
    for (size_t wI = 0; wI < n; ++wI) {
      Wall &w = walls_[cell.walls[wI]];
      if (cell.walls[wI] != tmpWall[wallIndex % n] &&
          w.hasVertex(tmpVertex[vertexIndex])) {
        ++wallIndex;
        tmpWall[wallIndex % n] = cell.walls[wI];
        ++numWallSorted;
        if (wallIndex < wallIndexStart + n - 1) {
          if (tmpVertex[vertexIndex] == w.vertex1) {
            vertexIndex = (vertexIndex + 1) % n;
            tmpVertex[vertexIndex] = w.vertex2;
            ++numVertexSorted;
          } else if (tmpVertex[vertexIndex] == w.vertex2) {
            vertexIndex = (vertexIndex + 1) % n;
            tmpVertex[vertexIndex] = w.vertex1;
            ++numVertexSorted;
          } else {
            sortError("Wrong vertex index in wall.");
          }
        } else {
          // Point vertexIndex at the initial vertex (not added again).
          vertexIndex = (vertexIndex + 1) % n;
        }
        // Verify or record this wall's orientation flag for this cell.
        int *ownSort, *otherSort;
        size_t otherCell;
        if (w.cell1 == cellI) {
          ownSort = &w.cellSort1;
          otherSort = &w.cellSort2;
          otherCell = w.cell2;
        } else if (w.cell2 == cellI) {
          ownSort = &w.cellSort2;
          otherSort = &w.cellSort1;
          otherCell = w.cell1;
        } else {
          sortError("Wall to be sorted does not belong to the cell");
          return;
        }
        if (*ownSort == 1) {
          if (w.vertex2 != tmpVertex[vertexIndex])
            sortError("Current sorting direction for wall does not comply "
                      "with previous.");
        } else if (*ownSort == -1) {
          if (w.vertex1 != tmpVertex[vertexIndex])
            sortError("Current sorting direction for wall does not comply "
                      "with previous.");
        } else { // not sorted before
          if (w.vertex2 == tmpVertex[vertexIndex]) {
            *ownSort = 1;
            if (!isBackground(otherCell))
              *otherSort = -1;
          } else if (w.vertex1 == tmpVertex[vertexIndex]) {
            *ownSort = -1;
            if (!isBackground(otherCell))
              *otherSort = 1;
          } else {
            sortError("tmpVertex not in wall sorted!");
          }
        }
        break;
      }
    }
  }

  if (numWallSorted != n || numVertexSorted != n) {
    std::cerr << "Tissue::sortWallAndVertex() Cell " << cellI << " has sorted "
              << numWallSorted << " walls out of " << n << " and "
              << numVertexSorted << " vertices out of " << n << std::endl;
    std::exit(-1);
  }
  cell.walls = tmpWall;
  cell.vertices = tmpVertex;
}

// --- connectivity check ----------------------------------------------------------

namespace {
// Error sink for checkConnectivity's parallel loops. Each partition appends
// to its own list, tagged with the index that produced the message, so the
// report can be replayed in the order a serial pass would have produced it
// however the loop was split. A clean tissue appends nothing, so the happy
// path allocates nothing.
struct ConnErrors {
  std::vector<std::vector<std::pair<size_t, std::string>>> parts;
  explicit ConnErrors(size_t n) : parts(n) {}
  void add(size_t part, size_t index, std::string msg) {
    parts[part].emplace_back(index, std::move(msg));
  }
  // Partition order first, then a stable sort by index: within one partition
  // indices already ascend, so this is exactly the serial order.
  std::vector<std::pair<size_t, std::string>> merged() const {
    std::vector<std::pair<size_t, std::string>> all;
    for (const auto &p : parts)
      all.insert(all.end(), p.begin(), p.end());
    std::stable_sort(all.begin(), all.end(),
                     [](const auto &a, const auto &b) { return a.first < b.first; });
    return all;
  }
};

// The per-cell checks below are quadratic in the cell's own size (hasVertex
// and the touches-exactly-2-walls count are both linear scans over the cell),
// so a cell costs microseconds, not nanoseconds, on a finely resampled mesh.
// The pool's default grain assumes the opposite and would run these serially.
constexpr size_t kCellGrain = 16;
constexpr size_t kFlatGrain = 4096;
} // namespace

// Runs after every accepted step, so it is on the hot path even though it is
// pure validation: on an 11980-wall per-face mesh the serial version was 14%
// of total wall time, almost all of it in the quadratic per-cell loop. It is
// const and reads only topology, so the loops parallelise directly; the
// checks themselves are unchanged.
void Tissue::checkConnectivity(int verbose) const {
  const size_t parts = ThreadPool::instance().numThreads();
  ConnErrors cellErr(parts), wallErr(parts), vertErr(parts), sortErr(parts);

  // Cells reference valid walls/vertices; no duplicated vertices.
  ThreadPool::instance().parallelFor(
      numCell(), kCellGrain, [&](size_t b, size_t e, size_t p) {
        for (size_t i = b; i < e; ++i) {
          const CellTopo &c = cells_[i];
          for (size_t w : c.walls)
            if (w >= numWall())
              cellErr.add(p, i, "cell " + std::to_string(i) +
                                    " references bad wall index");
          for (size_t k = 0; k < c.vertices.size(); ++k) {
            if (c.vertices[k] >= numVertex())
              cellErr.add(p, i, "cell " + std::to_string(i) +
                                    " references bad vertex index");
            for (size_t l = k + 1; l < c.vertices.size(); ++l)
              if (c.vertices[k] == c.vertices[l])
                cellErr.add(p, i, "cell " + std::to_string(i) +
                                      " lists duplicate vertex");
          }
        }
      });
  // Walls: cells in range or background, distinct; vertices in range, distinct.
  ThreadPool::instance().parallelFor(
      numWall(), kFlatGrain, [&](size_t b, size_t e, size_t p) {
        for (size_t i = b; i < e; ++i) {
          const Wall &w = walls_[i];
          if (!isBackground(w.cell1) && w.cell1 >= numCell())
            wallErr.add(p, i, "wall " + std::to_string(i) + " has bad cell1");
          if (!isBackground(w.cell2) && w.cell2 >= numCell())
            wallErr.add(p, i, "wall " + std::to_string(i) + " has bad cell2");
          if (w.cell1 == w.cell2)
            wallErr.add(p, i,
                        "wall " + std::to_string(i) + " connects a cell to itself");
          if (w.vertex1 >= numVertex() || w.vertex2 >= numVertex())
            wallErr.add(p, i, "wall " + std::to_string(i) + " has bad vertex");
          if (w.vertex1 == w.vertex2)
            wallErr.add(p, i, "wall " + std::to_string(i) + " has equal vertices");
        }
      });
  // Vertices: cells in range, never background, no duplicates; same for walls.
  ThreadPool::instance().parallelFor(
      numVertex(), kFlatGrain, [&](size_t b, size_t e, size_t p) {
        for (size_t i = b; i < e; ++i) {
          const VertexTopo &v = vertices_[i];
          for (size_t k = 0; k < v.cells.size(); ++k) {
            if (isBackground(v.cells[k]))
              vertErr.add(p, i, "vertex " + std::to_string(i) +
                                    " lists background cell");
            else if (v.cells[k] >= numCell())
              vertErr.add(p, i, "vertex " + std::to_string(i) +
                                    " references bad cell");
            for (size_t l = k + 1; l < v.cells.size(); ++l)
              if (v.cells[k] == v.cells[l])
                vertErr.add(p, i, "vertex " + std::to_string(i) +
                                      " lists duplicate cell");
          }
          for (size_t k = 0; k < v.walls.size(); ++k) {
            if (v.walls[k] >= numWall())
              vertErr.add(p, i, "vertex " + std::to_string(i) +
                                    " references bad wall");
            for (size_t l = k + 1; l < v.walls.size(); ++l)
              if (v.walls[k] == v.walls[l])
                vertErr.add(p, i, "vertex " + std::to_string(i) +
                                      " lists duplicate wall");
          }
        }
      });
  // Per cell: numWall==numVertex; wall endpoints are cell vertices; every
  // cell vertex touches exactly 2 of the cell's walls; sorting invariant.
  ThreadPool::instance().parallelFor(
      numCell(), kCellGrain, [&](size_t b, size_t e, size_t p) {
        for (size_t i = b; i < e; ++i) {
          const CellTopo &c = cells_[i];
          if (c.numWall() != c.numVertex())
            sortErr.add(p, i,
                        "cell " + std::to_string(i) + " has numWall != numVertex");
          for (size_t w : c.walls) {
            if (!c.hasVertex(walls_[w].vertex1) || !c.hasVertex(walls_[w].vertex2))
              sortErr.add(p, i, "cell " + std::to_string(i) +
                                    " wall endpoint is not a cell vertex");
          }
          for (size_t v : c.vertices) {
            size_t count = 0;
            for (size_t w : c.walls)
              if (walls_[w].hasVertex(v))
                ++count;
            if (count != 2)
              sortErr.add(p, i, "cell " + std::to_string(i) + " vertex " +
                                    std::to_string(v) +
                                    " does not touch exactly 2 cell walls");
          }
          const size_t n = c.numWall();
          for (size_t k = 0; k < n; ++k) {
            const Wall &w = walls_[c.walls[k]];
            int sort =
                (w.cell1 == i) ? w.cellSort1 : (w.cell2 == i) ? w.cellSort2 : 0;
            size_t vk = c.vertices[k];
            size_t vk1 = c.vertices[(k + 1) % n];
            if (sort == -1) {
              if (!(vk1 == w.vertex1 && vk == w.vertex2))
                sortErr.add(p, i, "cell " + std::to_string(i) +
                                      " sorting invariant broken");
            } else if (sort == 1) {
              if (!(vk == w.vertex1 && vk1 == w.vertex2))
                sortErr.add(p, i, "cell " + std::to_string(i) +
                                      " sorting invariant broken");
            }
          }
        }
      });

  size_t errors = 0;
  for (const ConnErrors *e : {&cellErr, &wallErr, &vertErr, &sortErr})
    for (const auto &[index, msg] : e->merged()) {
      ++errors;
      if (verbose)
        std::cerr << "Tissue::checkConnectivity() " << msg << std::endl;
    }
  if (errors) {
    std::cerr << "Tissue::checkConnectivity() " << errors
              << " errors found; aborting." << std::endl;
    std::exit(-1);
  }
}

// --- swap-with-last topology removal helpers -----------------------------------

namespace {
inline void replaceIndex(std::vector<size_t> &list, size_t from, size_t to) {
  for (size_t &x : list)
    if (x == from)
      x = to;
}
inline bool eraseIndex(std::vector<size_t> &list, size_t value) {
  auto it = std::find(list.begin(), list.end(), value);
  if (it == list.end())
    return false;
  list.erase(it);
  return true;
}
} // namespace

void Tissue::removeCellTopo(size_t index) {
  size_t last = numCell() - 1;
  if (index != last) {
    cells_[index] = std::move(cells_[last]);
    // Re-point references to the moved cell.
    for (size_t w : cells_[index].walls) {
      if (walls_[w].cell1 == last)
        walls_[w].cell1 = index;
      else if (walls_[w].cell2 == last)
        walls_[w].cell2 = index;
    }
    for (size_t v : cells_[index].vertices)
      replaceIndex(vertices_[v].cells, last, index);
  }
  cells_.pop_back();
}

void Tissue::removeWallTopo(size_t index) {
  size_t last = numWall() - 1;
  if (index != last) {
    walls_[index] = walls_[last];
    const Wall &w = walls_[index];
    if (!isBackground(w.cell1))
      replaceIndex(cells_[w.cell1].walls, last, index);
    if (!isBackground(w.cell2))
      replaceIndex(cells_[w.cell2].walls, last, index);
    replaceIndex(vertices_[w.vertex1].walls, last, index);
    replaceIndex(vertices_[w.vertex2].walls, last, index);
  }
  walls_.pop_back();
}

void Tissue::removeVertexTopo(size_t index) {
  size_t last = numVertex() - 1;
  if (index != last) {
    vertices_[index] = std::move(vertices_[last]);
    for (size_t c : vertices_[index].cells)
      replaceIndex(cells_[c].vertices, last, index);
    for (size_t w : vertices_[index].walls) {
      if (walls_[w].vertex1 == last)
        walls_[w].vertex1 = index;
      else if (walls_[w].vertex2 == last)
        walls_[w].vertex2 = index;
    }
  }
  vertices_.pop_back();
}

// --- cell removal -----------------------------------------------------------------

void Tissue::removeCell(size_t cellI, Matrix &cellData, Matrix &wallData,
                        Matrix &vertexData, Matrix &cellDerivs,
                        Matrix &wallDerivs, Matrix &vertexDerivs) {
  assert(cellI < numCell());
  // Mark boundary walls (other side background) for removal; interior walls
  // become boundary walls of the neighbor.
  std::vector<size_t> wallRemove;
  for (size_t w : cells_[cellI].walls) {
    Wall &wl = walls_[w];
    if ((wl.cell1 == cellI && isBackground(wl.cell2)) ||
        (wl.cell2 == cellI && isBackground(wl.cell1))) {
      wallRemove.push_back(w);
    } else if (wl.cell1 == cellI) {
      wl.cell1 = kBackground;
    } else if (wl.cell2 == cellI) {
      wl.cell2 = kBackground;
    } else {
      std::cerr << "Tissue::removeCell() wall not connected to cell" << std::endl;
      std::exit(-1);
    }
  }
  // Detach the cell (and its dying walls) from its vertices; drop vertices
  // left with no cells and no walls (swap-with-last, mirrored in the data).
  // Iterate positions over the live vertex list: swap fixes propagate into it.
  for (size_t k = 0; k < cells_[cellI].numVertex(); ++k) {
    size_t vI = cells_[cellI].vertices[k];
    VertexTopo &v = vertices_[vI];
    eraseIndex(v.cells, cellI);
    for (size_t w : wallRemove)
      eraseIndex(v.walls, w);
    if (v.cells.empty() && v.walls.empty()) {
      assert(vI < vertexData.rows());
      vertexData.removeRowSwap(vI);
      vertexDerivs.removeRowSwap(vI);
      removeVertexTopo(vI);
      std::cerr << "Vertex " << vI << " removed" << std::endl;
    } else if (v.cells.empty() || v.walls.empty()) {
      std::cerr << "Tissue::removeCell() strange vertex." << std::endl;
      std::exit(-1);
    }
  }
  // Remove walls now bounded by the dying cell and background. Re-check the
  // live list each iteration (swap-with-last may reshuffle it).
  for (size_t k = 0; k < cells_[cellI].numWall(); ++k) {
    size_t wI = cells_[cellI].walls[k];
    Wall &wl = walls_[wI];
    if ((isBackground(wl.cell1) && wl.cell2 == cellI) ||
        (isBackground(wl.cell2) && wl.cell1 == cellI)) {
      assert(wI < wallData.rows());
      wallData.removeRowSwap(wI);
      wallDerivs.removeRowSwap(wI);
      removeWallTopo(wI);
      // A removed wall's stale entry may now alias the moved wall; the
      // condition re-check above keeps this loop safe (legacy behavior).
    }
  }
  // Remove the cell itself.
  assert(cellI < cellData.rows());
  cellData.removeRowSwap(cellI);
  cellDerivs.removeRowSwap(cellI);
  removeCellTopo(cellI);

  assert(cellData.rows() == numCell());
  assert(wallData.rows() == numWall());
  assert(vertexData.rows() == numVertex());
}

// --- cell division ------------------------------------------------------------------

void Tissue::divideCell(size_t cellI, size_t wI, size_t w3I, Vec3 v1Pos,
                        Vec3 v2Pos, Matrix &cellData, Matrix &wallData,
                        Matrix &vertexData, Matrix &cellDerivs,
                        Matrix &wallDerivs, Matrix &vertexDerivs,
                        const std::vector<size_t> &volumeChangeList,
                        double threshold) {
  const size_t Nc = numCell(), Nw = numWall(), Nv = numVertex();
  const size_t i = cellI;
  const size_t dimension = vertexData.cols();
  // Global indices of the two cut walls (cell-local wI/w3I refer to the
  // mother's wall list, which stays unchanged until the end).
  const size_t wIg = cells_[i].walls[wI];
  const size_t w3Ig = cells_[i].walls[w3I];

  // Nudge cut points away from existing vertices if within threshold.
  if (threshold >= 0.0) {
    const Wall &w1 = walls_[wIg];
    const Wall &w2 = walls_[w3Ig];
    double w1L = wallLengthFromVertices(wIg, vertexData);
    double w2L = wallLengthFromVertices(w3Ig, vertexData);
    double t1 = 0.0, t2 = 0.0;
    for (size_t d = 0; d < dimension; ++d) {
      t1 += (v1Pos[d] - vertexData[w1.vertex2][d]) *
            (v1Pos[d] - vertexData[w1.vertex2][d]);
      t2 += (v2Pos[d] - vertexData[w2.vertex2][d]) *
            (v2Pos[d] - vertexData[w2.vertex2][d]);
    }
    t1 = std::sqrt(t1) / w1L;
    t2 = std::sqrt(t2) / w2L;
    assert(t1 >= 0.0 && t1 <= 1.0);
    assert(t2 >= 0.0 && t2 <= 1.0);
    if (t1 < threshold) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from " << t1 << " to "
                << threshold << std::endl;
      for (size_t d = 0; d < dimension; ++d)
        v1Pos[d] = vertexData[w1.vertex2][d] +
                   threshold * (vertexData[w1.vertex1][d] - vertexData[w1.vertex2][d]);
    } else if (t1 > 1.0 - threshold) {
      std::cerr << "Tissue::divideCell() Moving vertex 1 from " << t1 << " to "
                << 1.0 - threshold << std::endl;
      for (size_t d = 0; d < dimension; ++d)
        v1Pos[d] = vertexData[w1.vertex1][d] +
                   threshold * (vertexData[w1.vertex2][d] - vertexData[w1.vertex1][d]);
    }
    if (t2 < threshold) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from " << t2 << " to "
                << threshold << std::endl;
      for (size_t d = 0; d < dimension; ++d)
        v2Pos[d] = vertexData[w2.vertex2][d] +
                   threshold * (vertexData[w2.vertex1][d] - vertexData[w2.vertex2][d]);
    } else if (t2 > 1.0 - threshold) {
      std::cerr << "Tissue::divideCell() Moving vertex 2 from " << t2 << " to "
                << 1.0 - threshold << std::endl;
      for (size_t d = 0; d < dimension; ++d)
        v2Pos[d] = vertexData[w2.vertex1][d] +
                   threshold * (vertexData[w2.vertex2][d] - vertexData[w2.vertex1][d]);
    }
  }

  // New cell Nc: full copy of the mother's variables (rescaled at the end
  // when volumeChangeList is non-empty).
  cells_.push_back(cells_[i]); // wall/vertex lists overwritten below
  cellData.appendRowCopy(i);
  cellDerivs.appendRowCopy(i); // legacy copies row 0; content is irrelevant

  // New vertices Nv (at v1Pos) and Nv+1 (at v2Pos).
  vertices_.emplace_back();
  vertices_.emplace_back();
  vertexData.appendRow(std::span<const double>(v1Pos.data(), dimension));
  vertexData.appendRow(std::span<const double>(v2Pos.data(), dimension));
  vertexDerivs.appendRow(dimension);
  vertexDerivs.appendRow(dimension);

  // New wall Nw between the daughters; its non-length variables are a copy of
  // wall 0's row (legacy resize-fill quirk, preserved).
  walls_.emplace_back();
  wallData.appendRowCopy(0);
  {
    double tmpLength = 0.0;
    for (size_t d = 0; d < dimension; ++d)
      tmpLength += (v1Pos[d] - v2Pos[d]) * (v1Pos[d] - v2Pos[d]);
    wallData[Nw][0] = std::sqrt(tmpLength);
  }

  // Wall Nw+1 continues wI; split resting length by the geometric fraction.
  walls_.push_back(walls_[wIg]);
  {
    double oldL = wallData[wIg][0];
    size_t v1w = walls_[wIg].vertex1, v2w = walls_[wIg].vertex2;
    double tmpLength = 0.0, tmpLengthFrac = 0.0;
    for (size_t d = 0; d < dimension; ++d) {
      tmpLengthFrac += (v1Pos[d] - vertexData[v1w][d]) * (v1Pos[d] - vertexData[v1w][d]);
      tmpLength += (vertexData[v2w][d] - vertexData[v1w][d]) *
                   (vertexData[v2w][d] - vertexData[v1w][d]);
    }
    double lengthFrac = std::sqrt(tmpLengthFrac) / std::sqrt(tmpLength);
    wallData.appendRowCopy(wIg);
    wallData[wIg][0] = lengthFrac * oldL;
    wallData[Nw + 1][0] = oldL - wallData[wIg][0];
  }

  // Wall Nw+2 continues w3I.
  walls_.push_back(walls_[w3Ig]);
  {
    double oldL = wallData[w3Ig][0];
    size_t v1w = walls_[w3Ig].vertex1, v2w = walls_[w3Ig].vertex2;
    double tmpLength = 0.0, tmpLengthFrac = 0.0;
    for (size_t d = 0; d < dimension; ++d) {
      tmpLengthFrac += (v2Pos[d] - vertexData[v1w][d]) * (v2Pos[d] - vertexData[v1w][d]);
      tmpLength += (vertexData[v2w][d] - vertexData[v1w][d]) *
                   (vertexData[v2w][d] - vertexData[v1w][d]);
    }
    double lengthFrac = std::sqrt(tmpLengthFrac) / std::sqrt(tmpLength);
    wallData.appendRowCopy(w3Ig);
    wallData[w3Ig][0] = lengthFrac * oldL;
    wallData[Nw + 2][0] = oldL - wallData[w3Ig][0];
  }
  wallDerivs.appendRow(wallData.rowSize(Nw));
  wallDerivs.appendRow(wallData.rowSize(Nw + 1));
  wallDerivs.appendRow(wallData.rowSize(Nw + 2));

  // Partition the mother's boundary cycle into the 'old' (keeps wI side from
  // vertex1) and 'new' halves by walking wall-to-wall from wI to w3I.
  std::vector<size_t> oldVIndex, newVIndex, oldWIndex, newWIndex;
  std::vector<char> usedWIndex(cells_[i].numWall(), 0);
  const CellTopo &mother = cells_[i];
  auto walkHalf = [&](size_t startLocal, size_t startGlobal, size_t startVertex,
                      std::vector<size_t> &wOut, std::vector<size_t> &vOut) {
    size_t tmpWIndex = startGlobal;
    size_t tmpVIndex = startVertex;
    size_t nextW = startLocal;
    do {
      wOut.push_back(tmpWIndex);
      vOut.push_back(tmpVIndex);
      usedWIndex[nextW] = 1;
      nextW = mother.numWall();
      size_t flag = 0;
      for (size_t w = 0; w < mother.numWall(); ++w) {
        if (!usedWIndex[w] && walls_[mother.walls[w]].hasVertex(tmpVIndex)) {
          nextW = w;
          flag++;
        }
      }
      if (flag != 1) {
        std::cerr << "Tissue::divideCell() " << flag
                  << " walls marked for next wall..." << std::endl;
        std::exit(-1);
      }
      tmpWIndex = mother.walls[nextW];
      tmpVIndex = walls_[mother.walls[nextW]].otherVertex(tmpVIndex);
    } while (nextW != w3I);
    // Which half of the second cut wall closes this cycle?
    if (walls_[w3Ig].vertex1 == vOut.back())
      wOut.push_back(w3Ig);
    else if (walls_[w3Ig].vertex2 == vOut.back())
      wOut.push_back(Nw + 2);
    else {
      std::cerr << "Tissue::divideCell() wrong last index (not in w3I)"
                << std::endl;
      std::exit(-1);
    }
  };
  // Old half: from wI through its vertex1 side.
  walkHalf(wI, wIg, walls_[wIg].vertex1, oldWIndex, oldVIndex);
  // New half: from the wI copy (Nw+1) through the vertex2 side.
  usedWIndex[wI] = 0;
  usedWIndex[w3I] = 0;
  walkHalf(wI, Nw + 1, walls_[Nw + 1].vertex2, newWIndex, newVIndex);

  oldVIndex.push_back(Nv);
  oldVIndex.push_back(Nv + 1);
  newVIndex.push_back(Nv);
  newVIndex.push_back(Nv + 1);
  oldWIndex.push_back(Nw);
  newWIndex.push_back(Nw);

  auto contains = [](const std::vector<size_t> &v, size_t x) {
    return std::find(v.begin(), v.end(), x) != v.end();
  };

  // Wire up the dividing wall.
  walls_[Nw].vertex1 = Nv;
  walls_[Nw].vertex2 = Nv + 1;
  walls_[Nw].cell1 = i;
  walls_[Nw].cell2 = Nc;
  walls_[Nw].cellSort1 = walls_[Nw].cellSort2 = 0;

  // Split wI: mother keeps [vertex1, Nv]; Nw+1 gets [Nv, old vertex2].
  walls_[wIg].vertex2 = Nv;
  walls_[Nw + 1].vertex1 = Nv;
  // Re-point Nw+1's mother side to the daughter, and tell the outside
  // neighbor (if any) about the new sub-wall and vertex.
  if (walls_[Nw + 1].cell1 == i) {
    walls_[Nw + 1].cell1 = Nc;
    if (walls_[Nw + 1].cell2 < Nc) {
      cells_[walls_[Nw + 1].cell2].walls.push_back(Nw + 1);
      cells_[walls_[Nw + 1].cell2].vertices.push_back(Nv);
    }
  } else if (walls_[Nw + 1].cell2 == i) {
    walls_[Nw + 1].cell2 = Nc;
    if (walls_[Nw + 1].cell1 < Nc) {
      cells_[walls_[Nw + 1].cell1].walls.push_back(Nw + 1);
      cells_[walls_[Nw + 1].cell1].vertices.push_back(Nv);
    }
  } else {
    std::cerr << "Tissue::divideCell() First wall not connected to dividing cell"
              << std::endl;
    std::exit(-1);
  }

  // Split w3I between whichever halves ended up with it.
  if (contains(oldWIndex, w3Ig)) {
    if (contains(oldVIndex, walls_[w3Ig].vertex1)) {
      walls_[w3Ig].vertex2 = Nv + 1;
      walls_[Nw + 2].vertex1 = Nv + 1;
    } else {
      walls_[w3Ig].vertex1 = Nv + 1;
      walls_[Nw + 2].vertex2 = Nv + 1;
    }
  } else {
    if (contains(oldVIndex, walls_[Nw + 2].vertex1)) {
      walls_[Nw + 2].vertex2 = Nv + 1;
      walls_[w3Ig].vertex1 = Nv + 1;
    } else {
      walls_[Nw + 2].vertex1 = Nv + 1;
      walls_[w3Ig].vertex2 = Nv + 1;
    }
  }

  // Whichever w3I half belongs to the new cell gets re-pointed; the outside
  // neighbor learns about Nw+2 and Nv+1.
  size_t newWallIndex = contains(newWIndex, w3Ig) ? w3Ig : Nw + 2;
  if (walls_[newWallIndex].cell1 == i) {
    walls_[newWallIndex].cell1 = Nc;
    if (walls_[newWallIndex].cell2 < Nc) {
      cells_[walls_[newWallIndex].cell2].walls.push_back(Nw + 2);
      cells_[walls_[newWallIndex].cell2].vertices.push_back(Nv + 1);
    }
  } else if (walls_[newWallIndex].cell2 == i) {
    walls_[newWallIndex].cell2 = Nc;
    if (walls_[newWallIndex].cell1 < Nc) {
      cells_[walls_[newWallIndex].cell1].walls.push_back(Nw + 2);
      cells_[walls_[newWallIndex].cell1].vertices.push_back(Nv + 1);
    }
  } else {
    std::cerr << "Tissue::divideCell() Second wall not connected to dividing cell"
              << std::endl;
    std::exit(-1);
  }

  // Vertex adjacency for the new-half vertices.
  for (size_t v = 0; v < newVIndex.size(); ++v) {
    size_t vIdx = newVIndex[v];
    if (vIdx == Nv) {
      VertexTopo &vt = vertices_[Nv];
      vt.cells = {i, Nc};
      // wI's outside neighbor also borders this vertex.
      if (walls_[wIg].cell1 == i || walls_[wIg].cell1 == Nc) {
        if (!isBackground(walls_[wIg].cell2))
          vt.cells.push_back(walls_[wIg].cell2);
      } else if (walls_[wIg].cell2 == i || walls_[wIg].cell2 == Nc) {
        if (!isBackground(walls_[wIg].cell1))
          vt.cells.push_back(walls_[wIg].cell1);
      } else {
        std::cerr << "Tissue::divideCell() Wall wI not connected to dividing cell"
                  << std::endl;
        std::exit(-1);
      }
      vt.walls = {wIg, Nw, Nw + 1};
    } else if (vIdx == Nv + 1) {
      VertexTopo &vt = vertices_[Nv + 1];
      vt.cells = {i, Nc};
      if (walls_[w3Ig].cell1 == i || walls_[w3Ig].cell1 == Nc) {
        if (!isBackground(walls_[w3Ig].cell2))
          vt.cells.push_back(walls_[w3Ig].cell2);
      } else if (walls_[w3Ig].cell2 == i || walls_[w3Ig].cell2 == Nc) {
        if (!isBackground(walls_[w3Ig].cell1))
          vt.cells.push_back(walls_[w3Ig].cell1);
      } else {
        std::cerr << "Tissue::divideCell() Wall w3I not connected to dividing cell"
                  << std::endl;
        std::exit(-1);
      }
      vt.walls = {w3Ig, Nw, Nw + 2};
    } else {
      // Interior vertices of the new half: mother -> daughter, and cut-wall
      // references remapped to the new halves where appropriate.
      replaceIndex(vertices_[vIdx].cells, i, Nc);
      for (size_t &w : vertices_[vIdx].walls) {
        if (w == wIg)
          w = Nw + 1;
        else if (w == w3Ig && !contains(newWIndex, w3Ig))
          w = Nw + 2;
      }
    }
  }
  // Old-half vertices connected to w3I switch to Nw+2 if the old cell did not
  // keep w3I.
  for (size_t v = 0; v < oldVIndex.size(); ++v) {
    if (oldVIndex[v] < Nv && !contains(oldWIndex, w3Ig))
      replaceIndex(vertices_[oldVIndex[v]].walls, w3Ig, Nw + 2);
  }

  // Final cell wall/vertex lists; daughter-half walls still pointing at the
  // mother are re-pointed.
  cells_[i].walls = oldWIndex;
  cells_[i].vertices = oldVIndex;
  for (size_t w : newWIndex) {
    if (w != Nw && walls_[w].cell1 == i)
      walls_[w].cell1 = Nc;
    else if (w != Nw && walls_[w].cell2 == i)
      walls_[w].cell2 = Nc;
  }
  cells_[Nc].walls = newWIndex;
  cells_[Nc].vertices = newVIndex;
  assert(walls_[Nw].cell1 == i);
  assert(walls_[Nw].cell2 == Nc);

  // Split size-dependent variables in proportion to daughter areas.
  if (!volumeChangeList.empty()) {
    sortWallAndVertex(i);
    sortWallAndVertex(Nc);
    double Vi = cellVolume(i, vertexData);
    double Vn = cellVolume(Nc, vertexData);
    double fi = Vi / (Vi + Vn);
    double fn = Vn / (Vi + Vn);
    for (size_t k : volumeChangeList) {
      cellData[i][k] *= fi;
      cellData[Nc][k] *= fn;
    }
  }
}

} // namespace tissue
