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
// An optional third index level enables the cell stress state to be computed
// and stored (per cell: [principal stress direction (x,y,z), anisotropy
// a = 1 - s2/s1, s1]), the input to the dynamic anisotropic material of
// Walia, Carter et al. (2024) Eq. 3 (see WallMechanics::FiberSpring).
// It is evaluated in update(), i.e. between solver steps rather than inside
// every derivative evaluation: CMT reorientation is a slow (~hour) process,
// and holding the material constant within a step keeps the adaptive solver
// from chasing its own feedback between RK stages. An optional third
// parameter sets the minimum simulated time between refreshes (0 = every
// step); since the stress pass costs about as much as ten derivative
// evaluations, a refresh interval well below the CMT response time makes it
// negligible. The
// stress tensor computation is the exact legacy code path
// (mechanicalTRBS.cc:843-1360): per-triangle Cauchy stress from the
// deformation gradient, rotated to the global frame, area-averaged over the
// cell, and diagonalized with the legacy Jacobi iteration.
//
#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"

namespace tissue {
namespace {

class VertexFromTRBScenterTriangulation : public Reaction {
public:
  VertexFromTRBScenterTriangulation(const ParameterList &p,
                                    const IndexLevels &i) {
    if (p.size() != 2 && p.size() != 3)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulation: uses two or three parameters, "
          "Young modulus, Poisson coefficient and an optional stress-state "
          "refresh interval (simulated time; 0 = every step).");
    bool ok = (i.size() == 2 || i.size() == 3) && i[0].size() == 1 &&
              i[1].size() == 1 && (i.size() == 2 || i[2].size() == 1);
    if (!ok)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulation: wall length index in level 0; "
          "start of center-triangulation cell variables in level 1; optional "
          "level 2 = start index for storing the cell stress state "
          "[dir x, dir y, dir z, anisotropy, sigma1] (5 cell variables). "
          "Storing at index 0 puts the direction on VTK's native cell-vector "
          "slot and the anisotropy on its length slot.");
    std::vector<std::string> ids{"Y_mod", "P_ratio"};
    if (p.size() == 3)
      ids.push_back("stress_interval");
    configure("VertexFromTRBScenterTriangulation", p, i, p.size(),
              std::vector<size_t>(i.size(), 1), std::move(ids));
  }

  // Forces only; the stress state is refreshed once per step in update().
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &,
              Matrix &vertexDerivs) override {
    evaluate(T, cellData, wallData, vertexData, cellDerivs, vertexDerivs,
             /*wantForces=*/true, /*wantStress=*/false);
  }

  void initiate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &cellDerivs, Matrix &,
                Matrix &vertexDerivs) override {
    // Prime the stress state so the first derivative evaluation already has a
    // material orientation (otherwise the first step is isotropic).
    if (numVariableIndexLevel() == 3)
      evaluate(T, cellData, wallData, vertexData, cellDerivs, vertexDerivs,
               false, true);
  }


  // The center-triangulation vertex lives in the cell row; its "derivative"
  // is a force, so it relaxes with the vertices rather than growing.
  void positionalCellVariables(std::vector<size_t> &out) const override {
    const size_t com = variableIndex(1, 0);
    out.push_back(com);
    out.push_back(com + 1);
    out.push_back(com + 2);
  }

  void update(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, double h) override {
    if (numVariableIndexLevel() != 3)
      return;
    if (numParameter() == 3 && parameter(2) > 0.0) {
      sinceStressUpdate_ += h;
      if (sinceStressUpdate_ < parameter(2))
        return;
      sinceStressUpdate_ = 0.0;
    }
    // Scratch sinks kept as members: the stress pass writes no forces, but
    // evaluate() shares one code path with derivs().
    if (!scratchCell_.sameShape(cellData))
      scratchCell_.reshapeLike(cellData);
    if (!scratchVertex_.sameShape(vertexData))
      scratchVertex_.reshapeLike(vertexData);
    evaluate(T, cellData, wallData, vertexData, scratchCell_, scratchVertex_,
             false, true);
  }

private:
  Matrix scratchCell_, scratchVertex_;
  double sinceStressUpdate_ = 0.0;

  void evaluate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &cellDerivs, Matrix &vertexDerivs,
                bool wantForces, bool wantStress) {
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
    const bool storeStress = wantStress && numVariableIndexLevel() == 3;
    const size_t stressIndex = storeStress ? variableIndex(2, 0) : 0;

    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numWall();
            double stressCellGlobal[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
            double totalRestingArea = 0.0;
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

              // Resting area: numerically stable Heron (edges sorted),
              // guarded against degenerate/inverted trial configurations so a
              // bad adaptive-solver trial step is rejected by error control
              // instead of poisoning the state with NaNs.
              double lhe[3] = {restingLength[0], restingLength[1],
                               restingLength[2]};
              std::sort(lhe, lhe + 3);
              const double aHe = lhe[2], bHe = lhe[1], cHe = lhe[0];
              const double heron = ((bHe + cHe) + aHe) * (-(aHe - bHe) + cHe) *
                                   ((aHe - bHe) + cHe) * ((bHe - cHe) + aHe);
              const double restingArea =
                  0.25 * std::sqrt(std::max(heron, 1e-12));

              // Resting angles (law of cosines) and stiffnesses.
              auto angleOf = [&](double la, double lb, double opposite) {
                double arg = (la * la + lb * lb - opposite * opposite) /
                             (2.0 * la * lb);
                arg = std::max(-1.0 + 1e-9, std::min(1.0 - 1e-9, arg));
                return std::acos(arg);
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
              if (wantForces) {
                for (size_t d = 0; d < 3; ++d) {
                  const double forceCom = f0c * (pos[1][d] - pos[0][d]) +
                                          f0v * (pos[2][d] - pos[0][d]);
                  const double forceV2 = f1c * (pos[0][d] - pos[1][d]) +
                                         f1v * (pos[2][d] - pos[1][d]);
                  const double forceV3 = f2c * (pos[0][d] - pos[2][d]) +
                                         f2v * (pos[1][d] - pos[2][d]);
                  cellDerivs[c][comIndex + d] += forceCom;
                  vOut[v2][d] += forceV2;
                  vOut[v3][d] += forceV3;
                }
              }

              if (storeStress) {
                // Per-triangle Cauchy stress in the element plane, rotated to
                // the global frame (legacy mechanicalTRBS.cc:843-1055).
                const double trE = (delta[1] * cot0 + delta[2] * cot1 +
                                    delta[0] * cot2) /
                                   (4.0 * restingArea);
                const double curAngle1 = angleOf(length[0], length[1], length[2]);
                const double Qa = std::cos(curAngle1) * length[0];
                const double Qc = std::sin(curAngle1) * length[0];
                const double Qb = length[1];
                const double restAngle1 =
                    angleOf(restingLength[0], restingLength[1], restingLength[2]);
                const double Pa = std::cos(restAngle1) * restingLength[0];
                const double Pc = std::sin(restAngle1) * restingLength[0];
                const double Pb = restingLength[1];
                const double shapeResting[3][2] = {
                    {0.0, 1.0 / Pc},
                    {-1.0 / Pb, (Pa - Pb) / (Pb * Pc)},
                    {1.0 / Pb, -Pa / (Pb * Pc)}};
                const double posLocal[3][2] = {{Qa, Qc}, {0, 0}, {Qb, 0}};
                double F[2][2] = {{0, 0}, {0, 0}};
                for (int ii = 0; ii < 3; ++ii) {
                  F[0][0] += posLocal[ii][0] * shapeResting[ii][0];
                  F[1][0] += posLocal[ii][1] * shapeResting[ii][0];
                  F[0][1] += posLocal[ii][0] * shapeResting[ii][1];
                  F[1][1] += posLocal[ii][1] * shapeResting[ii][1];
                }
                double Bc[2][2]; // left Cauchy-Green B = F F^T
                Bc[0][0] = F[0][0] * F[0][0] + F[0][1] * F[0][1];
                Bc[1][0] = F[1][0] * F[0][0] + F[1][1] * F[0][1];
                Bc[0][1] = F[0][0] * F[1][0] + F[0][1] * F[1][1];
                Bc[1][1] = F[1][0] * F[1][0] + F[1][1] * F[1][1];
                double B2[2][2];
                B2[0][0] = Bc[0][0] * Bc[0][0] + Bc[0][1] * Bc[1][0];
                B2[1][0] = Bc[1][0] * Bc[0][0] + Bc[1][1] * Bc[1][0];
                B2[0][1] = Bc[0][0] * Bc[0][1] + Bc[0][1] * Bc[1][1];
                B2[1][1] = Bc[1][0] * Bc[0][1] + Bc[1][1] * Bc[1][1];
                const double curArea = 0.25 * std::sqrt(std::max(
                    (length[0] + length[1] + length[2]) *
                        (-length[0] + length[1] + length[2]) *
                        (length[0] - length[1] + length[2]) *
                        (length[0] + length[1] - length[2]),
                    1e-12));
                double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                const double fac = curArea / restingArea;
                S[0][0] = fac * ((lambda * trE - mio / 2) * Bc[0][0] +
                                 (mio / 2) * B2[0][0]);
                S[1][0] = fac * ((lambda * trE - mio / 2) * Bc[1][0] +
                                 (mio / 2) * B2[1][0]);
                S[0][1] = fac * ((lambda * trE - mio / 2) * Bc[0][1] +
                                 (mio / 2) * B2[0][1]);
                S[1][1] = fac * ((lambda * trE - mio / 2) * Bc[1][1] +
                                 (mio / 2) * B2[1][1]);
                // Rotation local->global from the triangle's frame.
                double X[3], Bv[3], Z[3], Yv[3];
                double tA = 0, tB = 0;
                for (size_t d = 0; d < 3; ++d) {
                  X[d] = pos[2][d] - pos[1][d];
                  Bv[d] = pos[0][d] - pos[1][d];
                  tA += X[d] * X[d];
                  tB += Bv[d] * Bv[d];
                }
                tA = std::sqrt(tA);
                tB = std::sqrt(tB);
                for (size_t d = 0; d < 3; ++d) {
                  X[d] /= tA;
                  Bv[d] /= tB;
                }
                Z[0] = X[1] * Bv[2] - X[2] * Bv[1];
                Z[1] = X[2] * Bv[0] - X[0] * Bv[2];
                Z[2] = X[0] * Bv[1] - X[1] * Bv[0];
                double zn = std::sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
                for (size_t d = 0; d < 3; ++d)
                  Z[d] /= zn;
                Yv[0] = Z[1] * X[2] - Z[2] * X[1];
                Yv[1] = Z[2] * X[0] - Z[0] * X[2];
                Yv[2] = Z[0] * X[1] - Z[1] * X[0];
                double R[3][3];
                for (size_t d = 0; d < 3; ++d) {
                  R[d][0] = X[d];
                  R[d][1] = Yv[d];
                  R[d][2] = Z[d];
                }
                // S_global = R S R^T, accumulated area-weighted.
                for (int r = 0; r < 3; ++r)
                  for (int t = 0; t < 3; ++t) {
                    double v = 0.0;
                    for (int u = 0; u < 3; ++u)
                      for (int w = 0; w < 3; ++w)
                        v += R[r][u] * S[u][w] * R[t][w];
                    stressCellGlobal[r][t] += restingArea * v;
                  }
                totalRestingArea += restingArea;
              }
            }

            if (storeStress && totalRestingArea > 0.0) {
              for (int r = 0; r < 3; ++r)
                for (int t = 0; t < 3; ++t)
                  stressCellGlobal[r][t] /= totalRestingArea;
              // Jacobi diagonalization (legacy mechanicalTRBS.cc:1269-1334).
              double eig[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
              double pivot = 1.0;
              const double pi = 3.1415;
              int iterations = 0;
              // A symmetric 3x3 needs only a handful of Jacobi sweeps; the
              // cap bounds the cost when an eigenvalue pair is degenerate.
              while (pivot > 0.00001 && ++iterations < 20) {
                int I = 1, J = 0;
                pivot = std::fabs(stressCellGlobal[1][0]);
                if (std::fabs(stressCellGlobal[2][0]) > pivot) {
                  pivot = std::fabs(stressCellGlobal[2][0]);
                  I = 2;
                  J = 0;
                }
                if (std::fabs(stressCellGlobal[2][1]) > pivot) {
                  pivot = std::fabs(stressCellGlobal[2][1]);
                  I = 2;
                  J = 1;
                }
                double rotAngle;
                if (std::fabs(stressCellGlobal[I][I] - stressCellGlobal[J][J]) <
                    0.00001)
                  rotAngle = pi / 4;
                else
                  rotAngle = 0.5 * std::atan((2 * stressCellGlobal[I][J]) /
                                             (stressCellGlobal[J][J] -
                                              stressCellGlobal[I][I]));
                const double Si = std::sin(rotAngle);
                const double Co = std::cos(rotAngle);
                double rot[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
                rot[I][I] = Co;
                rot[J][J] = Co;
                rot[I][J] = Si;
                rot[J][I] = -Si;
                double tmp[3][3];
                for (int r = 0; r < 3; ++r)
                  for (int t = 0; t < 3; ++t) {
                    tmp[r][t] = 0.0;
                    for (int w = 0; w < 3; ++w)
                      tmp[r][t] += stressCellGlobal[r][w] * rot[w][t];
                  }
                for (int r = 0; r < 3; ++r)
                  for (int t = 0; t < 3; ++t) {
                    stressCellGlobal[r][t] = 0.0;
                    for (int w = 0; w < 3; ++w)
                      stressCellGlobal[r][t] += rot[w][r] * tmp[w][t];
                  }
                for (int r = 0; r < 3; ++r)
                  for (int t = 0; t < 3; ++t)
                    tmp[r][t] = eig[r][t];
                for (int r = 0; r < 3; ++r)
                  for (int t = 0; t < 3; ++t) {
                    eig[r][t] = 0.0;
                    for (int w = 0; w < 3; ++w)
                      eig[r][t] += tmp[r][w] * rot[w][t];
                  }
              }
              // In-plane principal pair: two largest eigenvalues (the shell
              // normal eigenvalue is ~0 in a membrane under tension).
              int order[3] = {0, 1, 2};
              double ev[3] = {stressCellGlobal[0][0], stressCellGlobal[1][1],
                              stressCellGlobal[2][2]};
              for (int r = 0; r < 3; ++r)
                for (int t = r + 1; t < 3; ++t)
                  if (ev[order[t]] > ev[order[r]])
                    std::swap(order[r], order[t]);
              const double s1 = ev[order[0]];
              const double s2 = ev[order[1]];
              double a = 0.0;
              if (std::fabs(s1) > 1e-12)
                a = 1.0 - std::min(s1, s2) / std::max(s1, s2);
              if (a < 0.0)
                a = 0.0;
              if (a > 1.0)
                a = 1.0;
              for (size_t d = 0; d < 3; ++d)
                cellData[c][stressIndex + d] = eig[d][order[0]];
              cellData[c][stressIndex + 3] = a;
              cellData[c][stressIndex + 4] = s1;
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
  // The center vertex is a position, not a concentration; see the same
  // override on VertexFromTRBScenterTriangulation.
  void positionalCellVariables(std::vector<size_t> &out) const override {
    const size_t com = variableIndex(0, 0);
    out.push_back(com);
    out.push_back(com + 1);
    out.push_back(com + 2);
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

// Dynamic anisotropic wall material of Walia, Carter et al. (2024), Eq. 3:
// the fiber (cellulose/CMT) part of the wall stiffness, Y_f, redistributes
// between the principal stress directions according to the current stress
// anisotropy a (stored per cell by VertexFromTRBScenterTriangulation):
//   g(a)        = a^n / ((1-a)^n k^n + a^n)
//   Y_principal = 0.5 (1 + g) Y_f     (along maximal stress: CMTs align with
//   Y_second    = 0.5 (1 - g) Y_f      stress, cellulose follows CMTs)
// Implemented as oriented reinforcement springs on the walls: each wall gets
// stiffness K = Y_f * [0.5(1+g) cos^2(t) + 0.5(1-g) sin^2(t)] * w, averaged
// over its adjacent cells, where t is the angle between the wall and the
// cell's principal stress direction and w = A_cell/(2 l_wall) the transverse
// width the wall represents. At a=0 this is the isotropic contribution
// Y_f/2 in every direction, so the total wall modulus is the matrix TRBS
// Y_m plus this term, exactly the paper's decomposition.
class WallMechanicsFiberSpring : public Reaction {
public:
  WallMechanicsFiberSpring(const ParameterList &p, const IndexLevels &i) {
    if (p.size() != 3)
      throw std::runtime_error(
          "WallMechanics::FiberSpring: uses three parameters (Y_fiber, "
          "K_hill, n_hill).");
    configure("WallMechanics::FiberSpring", p, i, 3, {1, 1},
              {"Y_fiber", "K_hill", "n_hill"});
    // level 0: wall resting length index; level 1: cell stress-state start
    // index ([dir x, dir y, dir z, a, s1], written by the TRBS reaction).
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &, Matrix &,
              Matrix &vertexDerivs) override {
    const size_t lengthIndex = variableIndex(0, 0);
    const size_t sIndex = variableIndex(1, 0);
    const double yFiber = parameter(0);
    const double kHill = parameter(1);
    const double nHill = parameter(2);
    const size_t dim = vertexData.cols();
    // Per-cell caches: area (cellVolume is O(cell walls)) and the Eq. 3
    // redistribution factor g(a), which depends only on the cell's stress
    // anisotropy — computing it here keeps the two pow() calls out of the
    // per-wall loop.
    areas_.resize(T.numCell());
    gFactor_.resize(T.numCell());
    parallelFor(T.numCell(), [&](size_t b, size_t e) {
      for (size_t c = b; c < e; ++c) {
        areas_[c] = T.cellVolume(c, vertexData);
        const double a =
            std::max(0.0, std::min(1.0 - 1e-9, cellData[c][sIndex]));
        const double an = std::pow(a, nHill);
        gFactor_[c] =
            an / (std::pow(1.0 - a, nHill) * std::pow(kHill, nHill) + an);
      }
    });
    parallelScatter1(
        T.numWall(), vertexDerivs, [&](size_t b, size_t e, Matrix &out) {
          for (size_t w = b; w < e; ++w) {
            const Wall &wall = T.wall(w);
            const size_t v1 = wall.vertex1;
            const size_t v2 = wall.vertex2;
            double u[3] = {0, 0, 0};
            double d = 0.0;
            for (size_t dd = 0; dd < dim; ++dd) {
              u[dd] = vertexData[v1][dd] - vertexData[v2][dd];
              d += u[dd] * u[dd];
            }
            d = std::sqrt(d);
            if (d <= 0.0)
              continue;
            for (size_t dd = 0; dd < dim; ++dd)
              u[dd] /= d;
            // Stiffness contributions from the adjacent cells.
            double K = 0.0;
            for (size_t cIdx : {wall.cell1, wall.cell2}) {
              if (Tissue::isBackground(cIdx))
                continue;
              const double g = gFactor_[cIdx];
              double nvec[3] = {cellData[cIdx][sIndex], cellData[cIdx][sIndex + 1],
                                cellData[cIdx][sIndex + 2]};
              double nn = std::sqrt(nvec[0] * nvec[0] + nvec[1] * nvec[1] +
                                    nvec[2] * nvec[2]);
              double cos2 = 0.0;
              if (nn > 1e-12) {
                double dot = (u[0] * nvec[0] + u[1] * nvec[1] + u[2] * nvec[2]) / nn;
                cos2 = dot * dot;
              } else {
                cos2 = 0.5; // no direction yet (first evaluation): isotropic
              }
              const double orient =
                  0.5 * (1.0 + g) * cos2 + 0.5 * (1.0 - g) * (1.0 - cos2);
              const double width = areas_[cIdx] / (2.0 * d);
              K += yFiber * orient * width;
            }
            const double L = wallData[w][lengthIndex];
            double coeff = K * (1.0 / L - 1.0 / d);
            if (d <= 0.0 && L <= 0.0)
              coeff = 0.0;
            for (size_t dd = 0; dd < dim; ++dd) {
              double div = (vertexData[v1][dd] - vertexData[v2][dd]) * coeff;
              out[v1][dd] -= div;
              out[v2][dd] += div;
            }
          }
        });
  }

private:
  std::vector<double> areas_;   // per-call cell area cache
  std::vector<double> gFactor_; // per-call Eq. 3 redistribution cache
};
TISSUE_REGISTER_REACTION(WallMechanicsFiberSpring, "WallMechanics::FiberSpring")

} // namespace
} // namespace tissue
