//
// Small bookkeeping reactions ported from legacy cellTime.cc, dilution.cc and
// turgorGrowth.cc.
//
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

// dT/dt = 1, giving a per-cell clock that division rules can read. Legacy
// *assigns* rather than accumulates, so any earlier contribution to this
// variable in the same step is discarded; kept, since that is what models
// using it were written against.
class CellTimeDerivative : public Reaction {
public:
  CellTimeDerivative(const ParameterList &p, const IndexLevels &i) {
    configure("CellTimeDerivative", p, i, 0, {1}, {});
  }
  void derivs(Tissue &T, Matrix &, Matrix &, Matrix &, Matrix &cellDerivs,
              Matrix &, Matrix &) override {
    const size_t cIndex = variableIndex(0, 0);
    for (size_t cell = 0; cell < T.numCell(); ++cell)
      cellDerivs[cell][cIndex] = 1.0;
  }
};
TISSUE_REGISTER_REACTION(CellTimeDerivative, "CellTimeDerivative")

// Dilution of cell concentrations by area change: dc/dt -= c (dA/dt)/A, with
// dA/dt from the shoelace formula differentiated through the vertex
// velocities. 2D only, as in legacy.
class DilutionFromVertexDerivs : public Reaction {
public:
  DilutionFromVertexDerivs(const ParameterList &p, const IndexLevels &i) {
    if (!p.empty())
      throw std::runtime_error("Dilution::FromVertexDerivs: uses no parameters.");
    if (i.size() != 1 || i[0].empty())
      throw std::runtime_error(
          "Dilution::FromVertexDerivs: level 0 lists the cell variables to "
          "dilute.");
    configure("Dilution::FromVertexDerivs", p, i, 0, {i[0].size()}, {});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
              Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2)
      throw std::runtime_error(
          "Dilution::FromVertexDerivs: requires a 2D tissue.");
    for (size_t n = 0; n < T.numCell(); ++n) {
      const CellTopo &cell = T.cell(n);
      const double area = T.cellVolume(n, vertexData, /*signedArea=*/true);
      double areaDerivs = 0.0;
      const size_t nv = cell.numVertex();
      for (size_t k = 0; k < nv; ++k) {
        const size_t vI = cell.vertices[k];
        const size_t vP = cell.vertices[(k + 1) % nv];
        areaDerivs += vertexData[vP][1] * vertexDerivs[vI][0] -
                      vertexData[vI][1] * vertexDerivs[vP][0] -
                      vertexData[vP][0] * vertexDerivs[vI][1] +
                      vertexData[vI][0] * vertexDerivs[vP][1];
      }
      const double fac = 0.5 * areaDerivs / area;
      for (size_t k = 0; k < numVariableIndex(0); ++k)
        cellDerivs[n][variableIndex(0, k)] -=
            cellData[n][variableIndex(0, k)] * fac;
    }
  }
};
TISSUE_REGISTER_REACTION(DilutionFromVertexDerivs, "Dilution::FromVertexDerivs",
                         "DilutionFromVertexDerivs")

// Water uptake driven by the difference between a maximal turgor and the
// current one, the latter inferred from how far the stored water volume
// exceeds the geometric cell volume. Scaled by the cell's perimeter.
class TurgorGrowthWaterVolume : public Reaction {
public:
  TurgorGrowthWaterVolume(const ParameterList &p, const IndexLevels &i) {
    if (i.empty() || i.size() > 2 || i[0].size() != 1 ||
        (i.size() == 2 && i[1].size() != 1))
      throw std::runtime_error(
          "TurgorGrowth::WaterVolume: level 0 = the stored water volume; "
          "optional level 1 = a cell variable to write the turgor into.");
    std::vector<size_t> shape{1};
    if (i.size() == 2)
      shape.push_back(1);
    configure("TurgorGrowth::WaterVolume", p, i, 5, shape,
              {"k_p", "P_max", "k_pp", "denyShrink_flag",
               "allowNegTurgor_flag"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &) override {
    const size_t vIndex = variableIndex(0, 0);
    const bool storeP = numVariableIndexLevel() == 2;
    for (size_t n = 0; n < T.numCell(); ++n) {
      const CellTopo &cell = T.cell(n);
      double totalLength = 0.0;
      for (size_t w : cell.walls)
        totalLength += T.wallLengthFromVertices(w, vertexData);
      const double volume = T.cellVolume(n, vertexData);
      double P = (cellData[n][vIndex] - volume) / volume;
      if (P < 0.0 && !parameter(4))
        P = 0.0;
      P *= parameter(2);
      if (storeP)
        cellData[n][variableIndex(1, 0)] = P;
      if (!parameter(3) || parameter(1) - P > 0.0)
        cellDerivs[n][vIndex] +=
            parameter(0) * (parameter(1) - P) * totalLength;
    }
  }
};
TISSUE_REGISTER_REACTION(TurgorGrowthWaterVolume, "TurgorGrowth::WaterVolume",
                         "WaterVolumeFromTurgor")

} // namespace
} // namespace tissue
