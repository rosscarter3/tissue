//
// CenterTriangulation::Initiate: extends each cell's variable row with a
// central point (3D) and internal edge resting lengths, used by the CT
// mechanics/growth reactions. Ported from legacy centerTriangulation.cc.
//
// Cell row layout after initiation (base = variableIndex(0,0) = numVariable):
//   single edge: [vars..., cx, cy, cz, L_0..L_{n-1}]
//   double edge: [vars..., cx, cy, cz, (La_0, Lb_0)...(La_{n-1}, Lb_{n-1}),
//                 growth components g_0..g_{n-1}]
//
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class CTInitiate : public Reaction {
public:
  CTInitiate(const ParameterList &p, const IndexLevels &i) {
    if (p.size() > 2)
      throw std::runtime_error(
          "CenterTriangulation::Initiate: uses zero, one or two parameters "
          "(overRide_flag, [doubleEdge_flag]).");
    if (p.size() == 2 && p[1] == 2.0)
      throw std::runtime_error(
          "CenterTriangulation::Initiate: doubleEdge_flag=2 is obsolete; use "
          "overRide_flag=0 doubleEdge_flag=1.");
    std::vector<std::string> ids;
    if (!p.empty())
      ids.push_back("overRide_flag");
    if (p.size() == 2)
      ids.push_back("doubleEdge_flag");
    configure("CenterTriangulation::Initiate", p, i, p.size(), {1},
              std::move(ids));
  }

  void initiate(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &cellDerivs, Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "CenterTriangulation::Initiate requires a 3D tissue.");
    const size_t numVariable = T.numCellVariable();
    if (variableIndex(0, 0) != numVariable)
      throw std::runtime_error(
          "CenterTriangulation::Initiate: wrong index given as start index "
          "for additional variables (must equal the number of cell "
          "variables).");
    const bool doubleEdge = numParameter() == 2 && parameter(1) == 1.0;
    const bool fromScratch =
        T.cell(0).centerPosition.empty() ||
        (numParameter() && parameter(0) == 1.0);

    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numVertex();
      const size_t newSize =
          doubleEdge ? numVariable + 3 + 3 * n : numVariable + 3 + n;
      cellData.resizeRow(c, newSize);
      cellDerivs.resizeRow(c, newSize);

      Vec3 com;
      std::vector<double> edge(n);
      if (fromScratch) {
        com = T.cellPosition(c, vertexData);
        for (size_t k = 0; k < n; ++k) {
          auto v = vertexData[cell.vertices[k]];
          double d2 = 0.0;
          for (size_t d = 0; d < 3; ++d)
            d2 += (com[d] - v[d]) * (com[d] - v[d]);
          edge[k] = std::sqrt(d2);
        }
      } else {
        if (cell.centerPosition.empty())
          throw std::runtime_error(
              "CenterTriangulation::Initiate: init file has no center "
              "triangulation for cell " + std::to_string(c) +
              ". Did you forget to provide -centerTri_init to simulator?");
        for (size_t d = 0; d < 3; ++d)
          com[d] = cell.centerPosition[d];
        // Note: edgeLength was read in file vertex order and the cell was
        // sorted afterwards, matching legacy behavior.
        for (size_t k = 0; k < n; ++k)
          edge[k] = cell.edgeLength[k];
      }
      for (size_t d = 0; d < 3; ++d)
        cellData[c][numVariable + d] = com[d];
      if (doubleEdge) {
        for (size_t k = 0; k < n; ++k) {
          cellData[c][numVariable + 3 + 2 * k] = edge[k];
          cellData[c][numVariable + 3 + 2 * k + 1] = edge[k];
        }
        for (size_t k = 0; k < n; ++k)
          cellData[c][numVariable + 3 + 2 * n + k] = 0.0;
      } else {
        for (size_t k = 0; k < n; ++k)
          cellData[c][numVariable + 3 + k] = edge[k];
      }
    }
    std::cerr << "CenterTriangulation::Initiate: initiating "
              << (doubleEdge ? "double" : "single") << " edge CT from "
              << (fromScratch ? "scratch." : "init file.") << std::endl;
  }

  void derivs(Tissue &, Matrix &, Matrix &, Matrix &, Matrix &, Matrix &,
              Matrix &) override {}
};
TISSUE_REGISTER_REACTION(CTInitiate, "CenterTriangulation::Initiate")

} // namespace
} // namespace tissue
