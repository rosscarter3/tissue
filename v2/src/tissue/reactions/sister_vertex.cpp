//
// Sister-vertex reactions: pairs of vertices (e.g. duplicated vertices of
// per-face meshes) that are constrained to move together. Ported from legacy
// sisterVertex.cc.
//
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// Registers vertex pairs closer than d_max as sisters (all-pairs sweep).
class SisterVertexInitiateFromDistance : public Reaction {
public:
  SisterVertexInitiateFromDistance(const ParameterList &p, const IndexLevels &i) {
    configure("SisterVertex::InitiateFromDistance", p, i, 1, {}, {"d_max"});
  }
  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
                Matrix &, Matrix &) override {
    const size_t n = T.numVertex();
    const double dMax = parameter(0);
    size_t count = 0;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) {
        if (distance(vertexData[i], vertexData[j]) <= dMax) {
          T.addSisterVertex(i, j);
          ++count;
        }
      }
    std::cerr << "SisterVertex::InitiateFromDistance::initiate() added "
              << count << " sisters by distance rule." << std::endl;
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(SisterVertexInitiateFromDistance,
                         "SisterVertex::InitiateFromDistance")

// Sums the vertex-derivative contributions over each connected sister group
// and assigns the sum to every member (so they move in concert). Must be
// placed after all reactions that update vertex positions: this overwrites.
class SisterVertexCombineDerivatives : public Reaction {
public:
  SisterVertexCombineDerivatives(const ParameterList &p, const IndexLevels &i) {
    configure("SisterVertex::CombineDerivatives", p, i, 0, {}, {});
  }
  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
                Matrix &) override {
    // Transitive-closure grouping of sister pairs, ported verbatim from the
    // legacy algorithm (including the rescan-when-grown logic).
    groups_.clear();
    const size_t n = T.numSisterVertex();
    std::vector<std::array<size_t, 3>> tmp(n);
    for (size_t i = 0; i < n; ++i)
      tmp[i] = {T.sisterVertex(i, 0), T.sisterVertex(i, 1), 0};
    size_t counter = 1;
    for (size_t i = 0; i < n; ++i) {
      if (tmp[i][2] != 0)
        continue;
      tmp[i][2] = counter;
      groups_.push_back({tmp[i][0], tmp[i][1]});
      auto &group = groups_.back();
      size_t j = i + 1;
      size_t m = group.size();
      while (j < n) {
        if (tmp[j][2] == 0) {
          if (m < group.size()) { // group grew: restart the scan
            m = group.size();
            j = i + 1;
          }
          size_t ww = 0;
          for (size_t k = 0; k < m; ++k) {
            if (tmp[j][0] == group[k])
              ww += 1;
            if (tmp[j][1] == group[k])
              ww += 2;
          }
          if (ww == 1)
            group.push_back(tmp[j][1]);
          if (ww == 2)
            group.push_back(tmp[j][0]);
          if (ww != 0)
            tmp[j][2] = counter;
        }
        ++j;
      }
      ++counter;
    }
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const size_t dimension = vertexData.cols();
    for (const auto &group : groups_) {
      for (size_t d = 0; d < dimension; ++d) {
        double sum = 0.0;
        for (size_t v : group)
          sum += vertexDerivs[v][d];
        for (size_t v : group)
          vertexDerivs[v][d] = sum; // assignment, not accumulation
      }
    }
  }

private:
  std::vector<std::vector<size_t>> groups_;
};
TISSUE_REGISTER_REACTION(SisterVertexCombineDerivatives,
                         "SisterVertex::CombineDerivatives")

// Spring pulling each sister pair together; optionally breaks pairs whose
// distance exceeds BreakLength during update.
class SisterVertexSpring : public Reaction {
public:
  SisterVertexSpring(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1 && p.size() != 2)
      throw std::runtime_error("SisterVertex::Spring: uses one or two "
                               "parameters (K_spring, [BreakLength]).");
    configure("SisterVertex::Spring", p, i, p.size(), {},
              p.size() == 2
                  ? std::vector<std::string>{"K_spring", "BreakLength"}
                  : std::vector<std::string>{"K_spring"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &vertexData, Matrix &,
              Matrix &, Matrix &vertexDerivs) override {
    const double k = parameter(0);
    const size_t dimension = vertexData.cols();
    for (size_t s = 0; s < T.numSisterVertex(); ++s) {
      size_t v0 = T.sisterVertex(s, 0);
      size_t v1 = T.sisterVertex(s, 1);
      for (size_t d = 0; d < dimension; ++d) {
        double f = -k * (vertexData[v0][d] - vertexData[v1][d]);
        vertexDerivs[v0][d] += f;
        vertexDerivs[v1][d] -= f;
      }
    }
  }
  void update(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              double) override {
    if (numParameter() != 2)
      return;
    auto &sisters = T.sisterVertices();
    for (size_t s = sisters.size(); s-- > 0;) {
      if (distance(vertexData[sisters[s][0]], vertexData[sisters[s][1]]) >
          parameter(1)) {
        sisters[s] = sisters.back();
        sisters.pop_back();
      }
    }
  }
};
TISSUE_REGISTER_REACTION(SisterVertexSpring, "SisterVertex::Spring")


// Reads sister pairs from a file literally named "sister" in the working
// directory: a count followed by that many vertex-index pairs. The filename is
// hard-coded in legacy and kept so existing model directories still work.
class SisterVertexInitiateFromFile : public Reaction {
public:
  SisterVertexInitiateFromFile(const ParameterList &p, const IndexLevels &i) {
    configure("SisterVertex::InitiateFromFile", p, i, 0, {}, {});
  }
  void initiate(Tissue &T, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
                Matrix &) override {
    std::ifstream in("sister");
    if (!in)
      throw std::runtime_error(
          "SisterVertex::InitiateFromFile: cannot open file 'sister' (the "
          "filename is fixed; it is read from the working directory).");
    size_t n = 0;
    in >> n;
    auto &sisters = T.sisterVertices();
    sisters.resize(n);
    for (size_t i = 0; i < n; ++i) {
      size_t v1, v2;
      if (!(in >> v1 >> v2))
        throw std::runtime_error("SisterVertex::InitiateFromFile: file 'sister' "
                                 "ended before " + std::to_string(n) +
                                 " vertex pairs were read.");
      sisters[i] = {v1, v2};
    }
    std::cerr << "SisterVertex::InitiateFromFile::initiate() read " << n
              << " sisters from file sister." << std::endl;
  }
  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(SisterVertexInitiateFromFile,
                         "SisterVertex::InitiateFromFile")

// As SisterVertex::Spring, but the spring constant is scaled by the mean of a
// cell variable over the two cells the pair belongs to - so an adhesion
// molecule expressed per cell can switch the constraint on and off. Each
// sister vertex is assumed to border exactly one cell (the rest being
// background), which is what the per-face meshes this is written for give.
//
// FIX: legacy's gate reads `cellData[cell1] > 0 || cellData[cell1] > 0`,
// testing the first cell twice. The pair was therefore silently inert whenever
// the first vertex's cell was zero, however strongly the second expressed,
// which made the force depend on the order the pair happened to be listed in.
// The second test reads cell2 here, matching the comment in legacy ("check if
// either cell has non-zero value") and the symmetry of the force itself.
class SisterVertexSpringCellConc : public Reaction {
public:
  SisterVertexSpringCellConc(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 1 && p.size() != 2)
      throw std::runtime_error("SisterVertex::SpringCellConc: uses one or two "
                               "parameters (k_spring, [lengthFractionBreak]).");
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "SisterVertex::SpringCellConc: one cell variable index.");
    configure("SisterVertex::SpringCellConc", p, i, p.size(), {1},
              p.size() == 2
                  ? std::vector<std::string>{"K_spring", "BreakLength"}
                  : std::vector<std::string>{"K_spring"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    const size_t cIndex = variableIndex(0, 0);
    const size_t dimension = vertexData.cols();
    for (size_t s = 0; s < T.numSisterVertex(); ++s) {
      const size_t v0 = T.sisterVertex(s, 0);
      const size_t v1 = T.sisterVertex(s, 1);
      if (T.vertex(v0).cells.empty() || T.vertex(v1).cells.empty())
        continue; // no cell to read the concentration from
      const double c0 = cellData[T.vertex(v0).cells[0]][cIndex];
      const double c1 = cellData[T.vertex(v1).cells[0]][cIndex];
      if (c0 <= 0.0 && c1 <= 0.0)
        continue;
      const double factor = 0.5 * (c0 + c1) * parameter(0);
      for (size_t d = 0; d < dimension; ++d) {
        const double f = -factor * (vertexData[v0][d] - vertexData[v1][d]);
        vertexDerivs[v0][d] += f;
        vertexDerivs[v1][d] -= f;
      }
    }
  }
  void update(Tissue &T, Matrix &, Matrix &, Matrix &vertexData,
              double) override {
    if (numParameter() != 2)
      return;
    auto &sisters = T.sisterVertices();
    std::vector<size_t> remove;
    for (size_t i = 0; i < sisters.size(); ++i)
      if (distance(vertexData[sisters[i][0]], vertexData[sisters[i][1]]) >
          parameter(1))
        remove.push_back(i);
    // Reverse swap-with-last compaction, matching legacy's survivor ordering.
    for (size_t k = remove.size(); k-- > 0;) {
      const size_t i = remove[k];
      const size_t last = sisters.size() - 1;
      if (i != last)
        sisters[i] = sisters[last];
      sisters.pop_back();
    }
  }
};
TISSUE_REGISTER_REACTION(SisterVertexSpringCellConc,
                         "SisterVertex::SpringCellConc")

} // namespace
} // namespace tissue
