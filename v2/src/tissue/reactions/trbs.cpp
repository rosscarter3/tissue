//
// Triangular biquadratic springs (TRBS) on center-triangulated cells: proper
// 2D elasticity (Young's modulus + Poisson ratio) for polygonal cells
// embedded in 3D, following Delingette (2008) and the legacy
// VertexFromTRBScenterTriangulation (mechanicalTRBS.cc). Each cell is a fan
// of triangles (center, vertex k, vertex k+1); per triangle the biquadratic
// edge strains give tensile and angular stiffnesses from the Lame
// coefficients, and forces on the center (a cell variable) and the two
// vertices.
//
// The optional legacy levels for storing strain/stress eigendirections (MT
// feedback) are not ported; only the isotropic 2-level form is accepted.
//
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class VertexFromTRBScenterTriangulation : public Reaction {
public:
  VertexFromTRBScenterTriangulation(const ParameterList &p,
                                    const IndexLevels &i) {
    if (p.size() != 2)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulation: uses two parameters, Young "
          "modulus and Poisson coefficient.");
    if (i.size() != 2 || i[0].size() != 1 || i[1].size() != 1)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulation: wall length index in level 0; "
          "start of center-triangulation cell variables in level 1. (The "
          "legacy 4-level strain/stress-direction storage is not ported.)");
    configure("VertexFromTRBScenterTriangulation", p, i, 2, {1, 1},
              {"Y_mod", "P_ratio"});
  }

  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &,
              Matrix &vertexDerivs) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulation requires a 3D tissue.");
    const size_t wallLengthIndex = variableIndex(0, 0);
    const size_t comIndex = variableIndex(1, 0);
    const size_t lengthInternalIndex = comIndex + 3;
    const double young = parameter(0);
    const double poisson = parameter(1);
    const double lambda = young * poisson / (1 - poisson * poisson);
    const double mio = young / (1 + poisson);

    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numWall();
            for (size_t k = 0; k < n; ++k) {
              const size_t k1 = (k + 1) % n;
              const size_t v2 = cell.vertices[k];
              const size_t v3 = cell.vertices[k1];
              const size_t w2 = cell.walls[k];

              // Triangle nodes: 0 = center (com), 1 = vertex(k), 2 = vertex(k+1)
              double pos[3][3];
              for (size_t d = 0; d < 3; ++d) {
                pos[0][d] = cellData[c][comIndex + d];
                pos[1][d] = vertexData[v2][d];
                pos[2][d] = vertexData[v3][d];
              }
              double restingLength[3] = {
                  cellData[c][lengthInternalIndex + k],
                  wallData[w2][wallLengthIndex],
                  cellData[c][lengthInternalIndex + k1]};
              auto dist = [&](int a2, int b2) {
                double s = 0;
                for (size_t d = 0; d < 3; ++d) {
                  double diff = pos[a2][d] - pos[b2][d];
                  s += diff * diff;
                }
                return std::sqrt(s);
              };
              double length[3] = {dist(0, 1), dist(1, 2), dist(0, 2)};

              // Resting area: numerically stable Heron (edges sorted).
              double lhe[3] = {restingLength[0], restingLength[1],
                               restingLength[2]};
              std::sort(lhe, lhe + 3);
              const double aHe = lhe[2], bHe = lhe[1], cHe = lhe[0];
              const double restingArea =
                  0.25 * std::sqrt(((bHe + cHe) + aHe) *
                                   (-(aHe - bHe) + cHe) * ((aHe - bHe) + cHe) *
                                   ((bHe - cHe) + aHe));

              // Resting angles (law of cosines) and stiffnesses.
              auto angleOf = [&](double la, double lb, double opposite) {
                return std::acos((la * la + lb * lb - opposite * opposite) /
                                 (2.0 * la * lb));
              };
              const double angle0 =
                  angleOf(restingLength[0], restingLength[2], restingLength[1]);
              const double angle1 =
                  angleOf(restingLength[0], restingLength[1], restingLength[2]);
              const double angle2 =
                  angleOf(restingLength[1], restingLength[2], restingLength[0]);
              const double temp = 1.0 / (restingArea * 16.0);
              const double cot0 = 1.0 / std::tan(angle0);
              const double cot1 = 1.0 / std::tan(angle1);
              const double cot2 = 1.0 / std::tan(angle2);
              const double tensile[3] = {
                  (2 * cot2 * cot2 * (lambda + mio) + mio) * temp,
                  (2 * cot0 * cot0 * (lambda + mio) + mio) * temp,
                  (2 * cot1 * cot1 * (lambda + mio) + mio) * temp};
              const double angular[3] = {
                  (2 * cot1 * cot2 * (lambda + mio) - mio) * temp,
                  (2 * cot0 * cot2 * (lambda + mio) - mio) * temp,
                  (2 * cot0 * cot1 * (lambda + mio) - mio) * temp};

              // Biquadratic strains.
              const double delta[3] = {
                  length[0] * length[0] - restingLength[0] * restingLength[0],
                  length[1] * length[1] - restingLength[1] * restingLength[1],
                  length[2] * length[2] - restingLength[2] * restingLength[2]};

              // Forces (exact legacy expressions, mechanicalTRBS.cc:1089-1169).
              const double f0c =
                  tensile[0] * delta[0] + angular[1] * delta[1] +
                  angular[0] * delta[2];
              const double f0v =
                  tensile[2] * delta[2] + angular[2] * delta[1] +
                  angular[0] * delta[0];
              const double f1c =
                  tensile[0] * delta[0] + angular[0] * delta[2] +
                  angular[1] * delta[1];
              const double f1v =
                  tensile[1] * delta[1] + angular[2] * delta[2] +
                  angular[1] * delta[0];
              const double f2c =
                  tensile[2] * delta[2] + angular[0] * delta[0] +
                  angular[2] * delta[1];
              const double f2v =
                  tensile[1] * delta[1] + angular[1] * delta[0] +
                  angular[2] * delta[2];
              for (size_t d = 0; d < 3; ++d) {
                const double forceCom =
                    f0c * (pos[1][d] - pos[0][d]) + f0v * (pos[2][d] - pos[0][d]);
                const double forceV2 =
                    f1c * (pos[0][d] - pos[1][d]) + f1v * (pos[2][d] - pos[1][d]);
                const double forceV3 =
                    f2c * (pos[0][d] - pos[2][d]) + f2v * (pos[1][d] - pos[2][d]);
                cellDerivs[c][comIndex + d] += forceCom;
                vOut[v2][d] += forceV2;
                vOut[v3][d] += forceV3;
              }
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(VertexFromTRBScenterTriangulation,
                         "VertexFromTRBScenterTriangulation")

// Turgor pressure on a closed center-triangulated shell: per CT triangle
// (center, vertex k, vertex k+1) the force P * A * n_hat is distributed
// equally to its three nodes, with the triangle's area vector from the
// sorted cycle orientation. The global sign is fixed each call from the
// mesh's signed volume so that P > 0 always inflates. (New in v2, replacing
// the PCA-normal Pressure3D reactions for CT shells.)
class Pressure3DCenterTriangulation : public Reaction {
public:
  Pressure3DCenterTriangulation(const ParameterList &p, const IndexLevels &i) {
    configure("Pressure3D::CenterTriangulation", p, i, 1, {1}, {"P_force"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "Pressure3D::CenterTriangulation requires a 3D tissue.");
    const size_t comIndex = variableIndex(0, 0);
    // Signed volume of the closed shell (divergence theorem over the fan).
    double signedVolume = 0.0;
    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numVertex();
      const double *com = &cellData[c][comIndex];
      for (size_t k = 0; k < n; ++k) {
        auto a = vertexData[cell.vertices[k]];
        auto b2 = vertexData[cell.vertices[(k + 1) % n]];
        signedVolume += com[0] * (a[1] * b2[2] - a[2] * b2[1]) +
                        com[1] * (a[2] * b2[0] - a[0] * b2[2]) +
                        com[2] * (a[0] * b2[1] - a[1] * b2[0]);
      }
    }
    const double pForce =
        signedVolume >= 0.0 ? parameter(0) : -parameter(0);

    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numVertex();
            const double *com = &cellData[c][comIndex];
            for (size_t k = 0; k < n; ++k) {
              const size_t va = cell.vertices[k];
              const size_t vb = cell.vertices[(k + 1) % n];
              double e1[3], e2[3];
              for (size_t d = 0; d < 3; ++d) {
                e1[d] = vertexData[va][d] - com[d];
                e2[d] = vertexData[vb][d] - com[d];
              }
              // area vector = 0.5 e1 x e2; per-node force = P * areaVec / 3
              const double fx = pForce * 0.5 / 3.0 * (e1[1] * e2[2] - e1[2] * e2[1]);
              const double fy = pForce * 0.5 / 3.0 * (e1[2] * e2[0] - e1[0] * e2[2]);
              const double fz = pForce * 0.5 / 3.0 * (e1[0] * e2[1] - e1[1] * e2[0]);
              cellDerivs[c][comIndex] += fx;
              cellDerivs[c][comIndex + 1] += fy;
              cellDerivs[c][comIndex + 2] += fz;
              vOut[va][0] += fx;
              vOut[va][1] += fy;
              vOut[va][2] += fz;
              vOut[vb][0] += fx;
              vOut[vb][1] += fy;
              vOut[vb][2] += fz;
            }
          }
        });
  }
};
TISSUE_REGISTER_REACTION(Pressure3DCenterTriangulation,
                         "Pressure3D::CenterTriangulation")

} // namespace
} // namespace tissue
