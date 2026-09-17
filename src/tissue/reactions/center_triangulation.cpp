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
    const size_t dim = vertexData.cols(); // 3D in legacy; 2D also supported
    const size_t numVariable = T.numCellVariable();
    if (variableIndex(0, 0) != numVariable)
      throw std::runtime_error(
          "CenterTriangulation::Initiate: wrong index given as start index "
          "for additional variables (must equal the number of cell "
          "variables).");
    const bool doubleEdge = numParameter() == 2 && parameter(1) == 1.0;
    if (doubleEdge && dim != 3)
      throw std::runtime_error(
          "CenterTriangulation::Initiate: double edge layout is 3D only.");
    const bool fromScratch =
        T.cell(0).centerPosition.empty() ||
        (numParameter() && parameter(0) == 1.0);

    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numVertex();
      const size_t newSize =
          doubleEdge ? numVariable + 3 + 3 * n : numVariable + dim + n;
      cellData.resizeRow(c, newSize);
      cellDerivs.resizeRow(c, newSize);

      Vec3 com;
      std::vector<double> edge(n);
      if (fromScratch) {
        com = T.cellPosition(c, vertexData);
        for (size_t k = 0; k < n; ++k) {
          auto v = vertexData[cell.vertices[k]];
          double d2 = 0.0;
          for (size_t d = 0; d < dim; ++d)
            d2 += (com[d] - v[d]) * (com[d] - v[d]);
          edge[k] = std::sqrt(d2);
        }
      } else {
        if (cell.centerPosition.empty())
          throw std::runtime_error(
              "CenterTriangulation::Initiate: init file has no center "
              "triangulation for cell " + std::to_string(c) +
              ". Did you forget to provide -centerTri_init to simulator?");
        for (size_t d = 0; d < dim; ++d)
          com[d] = cell.centerPosition[d];
        // Note: edgeLength was read in file vertex order and the cell was
        // sorted afterwards, matching legacy behavior.
        for (size_t k = 0; k < n; ++k)
          edge[k] = cell.edgeLength[k];
      }
      for (size_t d = 0; d < dim; ++d)
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
          cellData[c][numVariable + dim + k] = edge[k];
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

// Maxwell-type relaxation of the center-triangulation edge resting lengths
// toward their current geometric lengths:
//   dL_k/dt = k_relax (d_k - L_k)
// Combined with CenterTriangulation::EdgeSpring this makes the cell interior
// a viscoelastic body: instantaneously stiff (suppressing the shear/buckling
// mode of quad cells), fluid on timescales longer than 1/k_relax so slow
// growth is not resisted. (New in v2; used by the apical-hook model.)
class CTEdgeRelaxation : public Reaction {
public:
  CTEdgeRelaxation(const ParameterList &p, const IndexLevels &i) {
    configure("CenterTriangulation::EdgeRelaxation", p, i, 1, {1},
              {"k_relax"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t posIndex = variableIndex(0, 0);
    const size_t dim = vertexData.cols();
    const double kRelax = parameter(0);
    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t lengthIndex = posIndex + dim;
      for (size_t k = 0; k < cell.numVertex(); ++k) {
        size_t v = cell.vertices[k];
        double d2 = 0.0;
        for (size_t dd = 0; dd < dim; ++dd) {
          double diff = vertexData[v][dd] - cellData[c][posIndex + dd];
          d2 += diff * diff;
        }
        cellDerivs[c][lengthIndex + k] +=
            kRelax * (std::sqrt(d2) - cellData[c][lengthIndex + k]);
      }
    }
  }
};
TISSUE_REGISTER_REACTION(CTEdgeRelaxation, "CenterTriangulation::EdgeRelaxation")

} // namespace
} // namespace tissue
