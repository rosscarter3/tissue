//
// Transversely isotropic TRBS: the fibre-reinforced material of legacy's
// VertexFromTRBScenterTriangulationMT and its relatives (mechanicalTRBS.cc).
//
// The material has two Young's moduli - one along the fibre (microtubule)
// direction, one across it. The isotropic part of the element force is the
// ordinary TRBS force evaluated with the *transverse* moduli; the fibre
// contributes an additive correction built from the extra Lame pair
// (lambdaL - lambdaT, mioL - mioT) and the fibre direction pulled back into
// the element's rest frame. Both pieces live in trbs_core.h and are shared
// with the isotropic reactions.
//
// Legacy carries five variants of this, ~5700 live lines between them, each a
// copy of the same routine. What they actually vary is small: where the two
// moduli come from, whether the cell is center-triangulated, and what gets
// stored afterwards.
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/reactions/reaction.h"
#include "tissue/reactions/trbs_core.h"

namespace tissue {
namespace {

// Level-0 index layout, 11 entries (12 with the loosening compound).
enum MtIdx {
  kWallLength = 0,
  kMTDirection = 1,
  kStrainAniso = 2,
  kStressAniso = 3,
  kAreaRatio = 4,
  kYoungT = 5,
  kAnisoEnergy = 6,
  kYoungL = 7,
  kMisesStress = 8,
  kStressTensor = 9,
  kNormalVector = 10,
  kLoosening = 11
};

// Cell variables legacy reads by hard-coded number in some MF-flag modes.
// They are not configurable and not documented anywhere but the code.
constexpr size_t kAdHocConcIndex = 13; // MF flags 6, 7, 8, 9
constexpr size_t kAdHocFlagIndex = 40; // MF flag 0

struct Moduli {
  double youngL;
  double youngT;
};

struct Lame {
  double lambdaL, mioL, lambdaT, mioT;
};

// Plane strain (flag 0) or plane stress (flag 1). Note the shear modulus is
// mio = Y/(2(1+p)) here, half the convention the isotropic reactions use, so
// the stiffness coefficients below read (lambdaT + 2 mioT, 2 mioT) rather
// than (lambda + mio, mio) - the same quantity either way.
inline Lame lameOf(const Moduli &m, double poissonL, double poissonT,
                   bool planeStress) {
  Lame l;
  if (!planeStress) {
    l.lambdaL = m.youngL * poissonL / ((1 + poissonL) * (1 - 2 * poissonL));
    l.lambdaT = m.youngT * poissonT / ((1 + poissonT) * (1 - 2 * poissonT));
  } else {
    l.lambdaL = m.youngL * poissonL / (1 - poissonL * poissonL);
    l.lambdaT = m.youngT * poissonT / (1 - poissonT * poissonT);
  }
  l.mioL = m.youngL / (2 * (1 + poissonL));
  l.mioT = m.youngT / (2 * (1 + poissonT));
  return l;
}

class VertexFromTRBScenterTriangulationMT : public Reaction {
public:
  VertexFromTRBScenterTriangulationMT(const ParameterList &p,
                                      const IndexLevels &i) {
    if (p.size() != 11 && p.size() != 13)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationMT: uses 11 or 13 parameters "
          "(Y_matrix, Y_fibre, poisson_L, poisson_T, MF flag, neighbour "
          "weight, max stress/strain, plane-stress flag, MT angle, MT update "
          "flag, double-resting-length flag, [Hill K, Hill n]).");
    const bool ok =
        (i.size() == 2 || i.size() == 4) &&
        (i[0].size() == 11 || i[0].size() == 12) && i[1].size() == 1 &&
        (i.size() == 2 || (i[2].size() <= 3 && i[3].size() <= 2));
    if (!ok)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationMT: level 0 holds 11 or 12 cell "
          "variable indices (wall length, MT direction, strain anisotropy, "
          "stress anisotropy, area ratio, transverse modulus, anisotropic "
          "energy, longitudinal modulus, MT stress/strain, stress tensor, "
          "cell normal, [loosening compound]); level 1 the start of the "
          "center-triangulation cell variables; optional levels 2 and 3 store "
          "strain and stress directions.");
    if (p[2] < 0 || p[2] >= 0.5 || p[3] < 0 || p[3] >= 0.5)
      throw std::runtime_error("VertexFromTRBScenterTriangulationMT: Poisson "
                               "ratios must satisfy 0 <= p < 0.5.");
    if (p[7] != 0 && p[7] != 1)
      throw std::runtime_error("VertexFromTRBScenterTriangulationMT: "
                               "parameter 7 must be 0 (plane strain) or 1 "
                               "(plane stress).");
    if (p[9] < 0 || p[9] > 4 || p[9] != std::floor(p[9]))
      throw std::runtime_error("VertexFromTRBScenterTriangulationMT: the MT "
                               "update flag (parameter 9) must be 0-4.");
    std::vector<size_t> shape;
    for (const auto &lvl : i)
      shape.push_back(lvl.size());
    configure("VertexFromTRBScenterTriangulationMT", p, i, p.size(), shape,
              {"Y_mod_M", "Y_mod_F", "P_ratio_L", "P_ratio_T", "MF_flag",
               "neighbourweight", "stressmax", "plane_stress_flag",
               "TETA_anisotropy", "MT_update_flag", "double_length_flag",
               "Hill_K", "Hill_n"});
  }

  // The centre is a cell variable driven by a force, so it relaxes with the
  // vertices rather than being integrated as a concentration.
  void positionalCellVariables(std::vector<size_t> &out) const override {
    const size_t com = variableIndex(1, 0);
    out.push_back(com);
    out.push_back(com + 1);
    out.push_back(com + 2);
  }

  void initiate(Tissue &, Matrix &cellData, Matrix &, Matrix &vertexData,
                Matrix &, Matrix &, Matrix &) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationMT requires a 3D tissue.");
    // Legacy reads two cell variables by hard-coded number in some MF modes
    // and would run off the end of the row without saying so.
    const double mf = parameter(4);
    if (mf == 0 && cellData.cols() <= kAdHocFlagIndex)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationMT: MF flag 0 reads cell variable "
          "40 (a hard-coded legacy switch that halves the fibre modulus when "
          "it equals 100), but this tissue has fewer cell variables.");
    if ((mf >= 6 && mf <= 9) && cellData.cols() <= kAdHocConcIndex)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationMT: MF flags 6-9 read cell "
          "variable 13 (a hard-coded legacy concentration), but this tissue "
          "has fewer cell variables.");
  }

  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &,
              Matrix &vertexDerivs) override {
    const size_t comIndex = variableIndex(1, 0);
    const size_t lengthInternalIndex = comIndex + 3;
    const size_t wallLengthIndex = variableIndex(0, kWallLength);
    const size_t mtIndex = variableIndex(0, kMTDirection);
    const bool planeStress = parameter(7) == 1.0;
    const bool doubleRest = parameter(10) == 1.0;

    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numWall();
      if (cell.numVertex() != n)
        throw std::runtime_error("VertexFromTRBScenterTriangulationMT: needs "
                                 "the same number of vertices and walls.");

      const Moduli mod = moduliFor(cellData, c);
      if (parameter(4) == 1.0) // this mode reports the transverse modulus back
        cellData[c][variableIndex(0, kYoungT)] = mod.youngT;
      const Lame lame = lameOf(mod, parameter(2), parameter(3), planeStress);
      const double dLambda = lame.lambdaL - lame.lambdaT;
      const double dMio = lame.mioL - lame.mioT;

      if (parameter(9) == 1.0)
        setDirectionFromAngle(cellData, c);
      const double dir[3] = {cellData[c][mtIndex], cellData[c][mtIndex + 1],
                             cellData[c][mtIndex + 2]};

      for (size_t k = 0; k < n; ++k) {
        const size_t k1 = (k + 1) % n;
        const size_t v2 = cell.vertices[k];
        const size_t v3 = cell.vertices[k1];

        trbs::Element el;
        for (size_t d = 0; d < 3; ++d) {
          el.pos[0][d] = cellData[c][comIndex + d];
          el.pos[1][d] = vertexData[v2][d];
          el.pos[2][d] = vertexData[v3][d];
        }
        if (doubleRest) {
          // Two independent resting lengths per internal edge, and the wall's
          // own share added to the cell's.
          el.rest[0] = cellData[c][lengthInternalIndex + 2 * k + 1];
          el.rest[1] = wallData[cell.walls[k]][wallLengthIndex] +
                       cellData[c][lengthInternalIndex + 2 * n + k];
          el.rest[2] = cellData[c][lengthInternalIndex + 2 * k1];
        } else {
          el.rest[0] = cellData[c][lengthInternalIndex + k];
          el.rest[1] = wallData[cell.walls[k]][wallLengthIndex];
          el.rest[2] = cellData[c][lengthInternalIndex + k1];
        }
        el.cur[0] = trbs::nodeDistance(el, 0, 1);
        el.cur[1] = trbs::nodeDistance(el, 1, 2);
        el.cur[2] = trbs::nodeDistance(el, 0, 2);
        trbs::completeElement(el);

        double f[3][3];
        trbs::elementForces(
            el, trbs::stiffnessFrom(el, lame.lambdaT + 2 * lame.mioT,
                                    2 * lame.mioT),
            f);
        const trbs::LocalFrame lf = trbs::localFrameOf(el);
        double deltaF[3][3];
        trbs::anisotropicDeltaForce(el, lf, dLambda, dMio, dir, deltaF);
        for (int node = 0; node < 3; ++node)
          for (int d = 0; d < 3; ++d)
            f[node][d] += deltaF[node][d];

        // A near-degenerate ("sliver") element. Legacy reacts by *doubling*
        // the derivatives already accumulated on these three nodes and
        // dropping this element's own force, which is not a meaningful
        // remedy, but it is reachable and is reproduced rather than
        // reinterpreted. It warns, as legacy does.
        const bool sliver = std::fabs(1 - el.cos[0]) < 0.0001 ||
                            std::fabs(1 - el.cos[1]) < 0.0001 ||
                            std::fabs(1 - el.cos[2]) < 0.0001;
        if (sliver) {
          double t[3];
          for (int d = 0; d < 3; ++d)
            t[d] = cellDerivs[c][comIndex + d] + vertexDerivs[v2][d] +
                   vertexDerivs[v3][d];
          for (int d = 0; d < 3; ++d) {
            cellDerivs[c][comIndex + d] += t[d];
            vertexDerivs[v2][d] += t[d];
            vertexDerivs[v3][d] += t[d];
          }
          std::cerr << "VertexFromTRBScenterTriangulationMT::derivs() WARNING!"
                    << std::endl
                    << " there is a sliver in cell: " << c << std::endl;
          continue;
        }
        for (int d = 0; d < 3; ++d) {
          cellDerivs[c][comIndex + d] += f[0][d];
          vertexDerivs[v2][d] += f[1][d];
          vertexDerivs[v3][d] += f[2][d];
        }
      }
    }
  }

  // Everything the forces do not need: the cell's area-averaged true strain
  // and true stress, their principal values and directions, the anisotropy
  // measures those imply, the stored energies, and the cell normal. Other
  // reactions read these - FiberModel drives the longitudinal modulus from
  // the stress anisotropy stored here, which is what closes the microtubule
  // feedback loop.
  //
  // It repeats the element loop rather than sharing one with derivs(), as
  // legacy does: the two need different quantities, and the update runs once
  // per step while derivs runs once per solver stage.
  void update(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, double) override {
    const size_t comIndex = variableIndex(1, 0);
    const size_t lengthInternalIndex = comIndex + 3;
    const size_t wallLengthIndex = variableIndex(0, kWallLength);
    const size_t mtIndex = variableIndex(0, kMTDirection);
    const bool planeStress = parameter(7) == 1.0;
    const bool doubleRest = parameter(10) == 1.0;
    const double neighbourWeight = parameter(5);
    const double eps = 0.000001;
    const bool storeDirs = numVariableIndexLevel() == 4;

    for (size_t c = 0; c < T.numCell(); ++c) {
      const CellTopo &cell = T.cell(c);
      const size_t n = cell.numWall();
      const Moduli mod = moduliFor(cellData, c, /*inUpdate=*/true);
      const Lame lame = lameOf(mod, parameter(2), parameter(3), planeStress);
      const double dLambda = lame.lambdaL - lame.lambdaT;
      const double dMio = lame.mioL - lame.mioT;

      // Unlike derivs(), this branch has no four-cell special case.
      if (parameter(9) == 1.0) {
        cellData[c][mtIndex] = std::cos(parameter(8));
        cellData[c][mtIndex + 1] = std::sin(parameter(8));
        cellData[c][mtIndex + 2] = 0.0;
      }
      const double dir[3] = {cellData[c][mtIndex], cellData[c][mtIndex + 1],
                             cellData[c][mtIndex + 2]};

      double strainCell[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      double stressCell[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      double normalGlob[3] = {0, 0, 0};
      double totalRestArea = 0, totalArea = 0;
      double energyIso = 0, energyAniso = 0;

      for (size_t k = 0; k < n; ++k) {
        const size_t k1 = (k + 1) % n;
        trbs::Element el;
        for (size_t d = 0; d < 3; ++d) {
          el.pos[0][d] = cellData[c][comIndex + d];
          el.pos[1][d] = vertexData[cell.vertices[k]][d];
          el.pos[2][d] = vertexData[cell.vertices[k1]][d];
        }
        if (doubleRest) {
          el.rest[0] = cellData[c][lengthInternalIndex + 2 * k + 1];
          el.rest[1] = wallData[cell.walls[k]][wallLengthIndex] +
                       cellData[c][lengthInternalIndex + 2 * n + k];
          el.rest[2] = cellData[c][lengthInternalIndex + 2 * k1];
        } else {
          el.rest[0] = cellData[c][lengthInternalIndex + k];
          el.rest[1] = wallData[cell.walls[k]][wallLengthIndex];
          el.rest[2] = cellData[c][lengthInternalIndex + k1];
        }
        el.cur[0] = trbs::nodeDistance(el, 0, 1);
        el.cur[1] = trbs::nodeDistance(el, 1, 2);
        el.cur[2] = trbs::nodeDistance(el, 0, 2);
        trbs::completeElement(el);

        const trbs::LocalFrame lf = trbs::localFrameOf(el);
        double unusedForce[3][3];
        const trbs::AnisoTerms an = trbs::anisotropicDeltaForce(
            el, lf, dLambda, dMio, dir, unusedForce);

        const double area =
            0.25 * std::sqrt(std::max((el.cur[0] + el.cur[1] + el.cur[2]) *
                                          (-el.cur[0] + el.cur[1] + el.cur[2]) *
                                          (el.cur[0] - el.cur[1] + el.cur[2]) *
                                          (el.cur[0] + el.cur[1] - el.cur[2]),
                                      1e-12));
        // Left Cauchy-Green B and its square.
        double B[2][2];
        B[0][0] = lf.F[0][0] * lf.F[0][0] + lf.F[0][1] * lf.F[0][1];
        B[0][1] = lf.F[0][0] * lf.F[1][0] + lf.F[0][1] * lf.F[1][1];
        B[1][0] = lf.F[1][0] * lf.F[0][0] + lf.F[1][1] * lf.F[0][1];
        B[1][1] = lf.F[1][0] * lf.F[1][0] + lf.F[1][1] * lf.F[1][1];
        const double detB = B[0][0] * B[1][1] - B[1][0] * B[0][1];
        double B2[2][2];
        B2[0][0] = B[0][0] * B[0][0] + B[0][1] * B[1][0];
        B2[0][1] = B[0][0] * B[0][1] + B[0][1] * B[1][1];
        B2[1][0] = B[1][0] * B[0][0] + B[1][1] * B[1][0];
        B2[1][1] = B[1][0] * B[0][1] + B[1][1] * B[1][1];

        // Almansi (true) strain e = (I - B^-1)/2.
        double strainT[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        strainT[0][0] = 0.5 * (1 - B[1][1] / detB);
        strainT[0][1] = 0.5 * B[0][1] / detB;
        strainT[1][0] = 0.5 * B[1][0] / detB;
        strainT[1][1] = 0.5 * (1 - B[0][0] / detB);

        // True stress: isotropic part plus the fibre correction. Note the
        // area factor here is 1/detF, the inverse of the one the isotropic
        // reaction's Cauchy stress uses - these are two different code paths
        // in legacy and each is reproduced as it stands.
        const double areaFactor = el.restArea / area;
        double stressT[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int r = 0; r < 2; ++r)
          for (int t = 0; t < 2; ++t)
            stressT[r][t] =
                areaFactor * ((lame.lambdaT * an.trE - lame.mioT) * B[r][t] +
                              lame.mioT * B2[r][t]);
        double dSFt[2][2];
        for (int r = 0; r < 2; ++r)
          for (int t = 0; t < 2; ++t)
            dSFt[r][t] =
                an.deltaS[r][0] * lf.F[t][0] + an.deltaS[r][1] * lf.F[t][1];
        for (int r = 0; r < 2; ++r)
          for (int t = 0; t < 2; ++t)
            stressT[r][t] +=
                areaFactor * (lf.F[r][0] * dSFt[0][t] + lf.F[r][1] * dSFt[1][t]);

        // Both tensors rotated out to the global frame and area weighted.
        for (double (*M)[3] : {strainT, stressT}) {
          double tmp[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
          for (int r = 0; r < 3; ++r)
            for (int t = 0; t < 3; ++t)
              for (int w = 0; w < 3; ++w)
                tmp[r][t] += lf.R[r][w] * M[w][t];
          for (int r = 0; r < 3; ++r)
            for (int t = 0; t < 3; ++t) {
              M[r][t] = 0;
              for (int w = 0; w < 3; ++w)
                M[r][t] += tmp[r][w] * lf.R[t][w];
            }
        }
        for (int r = 0; r < 3; ++r)
          for (int t = 0; t < 3; ++t) {
            strainCell[r][t] += area * strainT[r][t];
            stressCell[r][t] += area * stressT[r][t];
          }
        for (int d = 0; d < 3; ++d)
          normalGlob[d] += area * lf.R[d][2]; // the element normal
        totalRestArea += el.restArea;
        totalArea += area;
        energyIso += ((lame.lambdaT / 2) * an.trE * an.trE + lame.mioT * an.I2) *
                     el.restArea;
        energyAniso +=
            ((dLambda / 2) * an.I4 * an.trE + dMio * an.I5) * el.restArea;
      }

      for (int r = 0; r < 3; ++r)
        for (int t = 0; t < 3; ++t) {
          strainCell[r][t] /= totalArea;
          stressCell[r][t] /= totalArea;
        }
      for (int d = 0; d < 3; ++d)
        normalGlob[d] /= totalArea;
      const double nrm = std::sqrt(normalGlob[0] * normalGlob[0] +
                                   normalGlob[1] * normalGlob[1] +
                                   normalGlob[2] * normalGlob[2]);
      if (nrm > 0)
        for (int d = 0; d < 3; ++d)
          normalGlob[d] /= nrm;

      const size_t stressTensorIndex = variableIndex(0, kStressTensor);
      const size_t normalIndex = variableIndex(0, kNormalVector);
      if (neighbourWeight > 0) {
        // The stress pass is redone below over the neighbour-averaged tensor,
        // so only the raw tensor is stored here.
        cellData[c][stressTensorIndex] = stressCell[0][0];
        cellData[c][stressTensorIndex + 1] = stressCell[1][1];
        cellData[c][stressTensorIndex + 2] = stressCell[2][2];
        cellData[c][stressTensorIndex + 3] = stressCell[0][1];
        cellData[c][stressTensorIndex + 4] = stressCell[0][2];
        cellData[c][stressTensorIndex + 5] = stressCell[1][2];
      }
      for (int d = 0; d < 3; ++d)
        cellData[c][normalIndex + d] = normalGlob[d];

      if (neighbourWeight == 0)
        storeStressState(cellData, c, stressCell, normalGlob, eps, storeDirs,
                         /*byMagnitude=*/true, /*writeMises=*/true);

      // Strain: diagonalize, rank by |value|, and store.
      double eigStrain[3][3];
      trbs::jacobiStrainStyle(strainCell, eigStrain, eps);
      int i1, i2, i3;
      rankPrincipal(strainCell, i1, i2, i3);
      const double s1 = strainCell[i1][i1];
      const double s2 = strainCell[i2][i2];
      cellData[c][variableIndex(0, kStrainAniso)] =
          std::fabs(s1) < eps ? 0.0 : 1 - std::fabs(s2 / s1);

      double perp[3] = {
          normalGlob[1] * eigStrain[2][i1] - normalGlob[2] * eigStrain[1][i1],
          normalGlob[2] * eigStrain[0][i1] - normalGlob[0] * eigStrain[2][i1],
          normalGlob[0] * eigStrain[1][i1] - normalGlob[1] * eigStrain[0][i1]};
      const double pn = std::sqrt(perp[0] * perp[0] + perp[1] * perp[1] +
                                  perp[2] * perp[2]);
      if (std::fabs(pn) < 0.0001) // maximal strain along the cell normal
        for (int d = 0; d < 3; ++d)
          perp[d] = eigStrain[d][i1];

      if (storeDirs) {
        const size_t nStrain = numVariableIndex(2);
        if (nStrain >= 1)
          storeDirection(cellData, c, variableIndex(2, 0), eigStrain, i1, s1);
        if (nStrain >= 2) {
          const size_t at = variableIndex(2, 1);
          for (int d = 0; d < 3; ++d)
            cellData[c][at + d] = perp[d];
          cellData[c][at + 3] = s2;
        }
        if (nStrain == 3)
          storeDirection(cellData, c, variableIndex(2, 2), eigStrain, i2, s2);
      }

      // Legacy normalizes the stored direction every step. It cancels in the
      // force (the pullback is normalized anyway), but it keeps the variable
      // readable by anything else that looks at it.
      const double mn =
          std::sqrt(cellData[c][mtIndex] * cellData[c][mtIndex] +
                    cellData[c][mtIndex + 1] * cellData[c][mtIndex + 1] +
                    cellData[c][mtIndex + 2] * cellData[c][mtIndex + 2]);
      if (mn > 0)
        for (int d = 0; d < 3; ++d)
          cellData[c][mtIndex + d] /= mn;

      cellData[c][variableIndex(0, kAreaRatio)] = totalRestArea;
      cellData[c][variableIndex(0, kYoungT)] = energyIso / totalRestArea;
      cellData[c][variableIndex(0, kAnisoEnergy)] = energyAniso / totalRestArea;
    }

    if (neighbourWeight > 0)
      neighbourStressPass(T, cellData, neighbourWeight, eps, storeDirs);
  }

private:
  // Young's moduli along and across the fibre. The MF flag (parameter 4)
  // selects between a constant anisotropic material, a modulus supplied by a
  // FiberModel reaction, and a set of ad-hoc variants the original author
  // added for particular figures - several of which read cell variables by
  // hard-coded number. They are ported as they stand.
  // `inUpdate` selects legacy's *other* material chain. The two differ in
  // two places and nowhere else: the update pass has no MF-flag-9 branch at
  // all (the moduli stay at 1 there), and its MF-flag-0 branch omits the
  // ad-hoc cell-variable-40 switch. Neither difference looks intentional, but
  // both change results, so both are reproduced.
  Moduli moduliFor(Matrix &cellData, size_t c, bool inUpdate = false) const {
    const double youngMatrix = parameter(0);
    double youngFiber = parameter(1);
    const double mf = parameter(4);
    Moduli m{1.0, 1.0};

    if (mf == 1) { // anisotropy supplied by FiberModel
      m.youngL = cellData[c][variableIndex(0, kYoungL)];
      m.youngT = 2 * youngMatrix + youngFiber - m.youngL;
      return m;
    }
    if (mf == 0) { // constant anisotropic material
      m.youngL = youngMatrix + youngFiber;
      m.youngT = youngMatrix;
      if (!inUpdate && cellData[c][kAdHocFlagIndex] == 100) { // ad-hoc switch
        m.youngL = youngMatrix + youngFiber / 2;
        m.youngT = youngMatrix + youngFiber / 2;
      }
      return m;
    }
    if (mf == 2) { // Y_matrix is the total stiffness here
      m.youngL = youngFiber;
      m.youngT = youngMatrix - m.youngL;
      return m;
    }
    if (mf == 3) { // fixed total stiffness, Y_fibre is the anisotropy
      const double total = youngMatrix, aniso = youngFiber;
      m.youngT = ((1 - aniso) / (2 - aniso)) * total;
      m.youngL = total - m.youngT;
      return m;
    }
    if (mf == 5) {
      youngFiber = cellData[c][variableIndex(0, kAnisoEnergy)];
      const double fiberL = cellData[c][variableIndex(0, kYoungL)];
      m.youngL = youngMatrix + fiberL;
      m.youngT = youngMatrix + youngFiber - fiberL;
      return m;
    }
    if (mf == 9 && inUpdate)
      return m; // legacy's update chain has no branch for this flag
    if (mf >= 6 && mf <= 9) {
      const double conc = cellData[c][kAdHocConcIndex];
      const double kConc = 0.005, nConc = 2.0;
      const double frac =
          1 - std::pow(conc, nConc) /
                  (std::pow(kConc, nConc) + std::pow(conc, nConc));
      if (mf == 6) {
        m.youngL = cellData[c][variableIndex(0, kYoungL)];
        m.youngT = 2 * youngMatrix + youngFiber - m.youngL;
        m.youngL *= 0.25 + 0.75 * frac;
        m.youngT *= 0.25 + 0.75 * frac;
        return m;
      }
      m.youngL = cellData[c][variableIndex(0, kYoungL)] - youngMatrix;
      m.youngT = 2 * youngMatrix + youngFiber - m.youngL - 2 * youngMatrix;
      if (mf == 9) {
        if (conc == 1.0 || conc == 0.0) {
          m.youngL *= 2;
          m.youngT *= 2;
        } else if (conc == 2.0) {
          m.youngL *= 0.5;
          m.youngT *= 0.5;
        }
        return m;
      }
      m.youngL *= frac;
      m.youngT *= frac;
      const double add =
          (mf == 8) ? youngMatrix * (0.2 + 0.8 * frac) : youngMatrix;
      m.youngL += add;
      m.youngT += add;
      return m;
    }
    if (mf == 10) { // FiberModel modulus, locally loosened by a compound
      m.youngL = cellData[c][variableIndex(0, kYoungL)];
      m.youngT = 2 * youngMatrix + youngFiber - m.youngL;
      const double A = cellData[c][variableIndex(0, kLoosening)];
      const double frac = 1 - std::pow(A, hillN()) /
                                  (std::pow(hillK(), hillN()) +
                                   std::pow(A, hillN()));
      m.youngL *= frac;
      m.youngT *= frac;
      return m;
    }
    // MF flag 4 has no branch in legacy either; both moduli stay at 1.
    return m;
  }

  // Rank the three diagonal entries by magnitude: i1 largest, then i2, i3.
  // Legacy writes this out longhand in four places and, in one of them,
  // with the i2/i3 assignment for the middle case transposed relative to the
  // others; the ordering it produces is the same either way because the pair
  // is immediately sorted by magnitude straight after.
  static void rankPrincipal(const double A[3][3], int &i1, int &i2, int &i3) {
    i1 = 0;
    if (std::fabs(A[1][1]) > std::fabs(A[i1][i1]))
      i1 = 1;
    if (std::fabs(A[2][2]) > std::fabs(A[i1][i1]))
      i1 = 2;
    i2 = (i1 + 1) % 3;
    i3 = (i1 + 2) % 3;
    if (std::fabs(A[i3][i3]) > std::fabs(A[i2][i2]))
      std::swap(i2, i3);
  }

  // As rankPrincipal, but on the signed values.
  static void rankPrincipalSigned(const double A[3][3], int &i1, int &i2,
                                  int &i3) {
    i1 = 0;
    if (A[1][1] > A[i1][i1])
      i1 = 1;
    if (A[2][2] > A[i1][i1])
      i1 = 2;
    i2 = (i1 + 1) % 3;
    i3 = (i1 + 2) % 3;
    if (A[i3][i3] > A[i2][i2])
      std::swap(i2, i3);
  }

  // Direction (3 components) then its principal value, at `at`.
  static void storeDirection(Matrix &cellData, size_t c, size_t at,
                             const double eig[3][3], int col, double value) {
    for (int d = 0; d < 3; ++d)
      cellData[c][at + d] = eig[d][col];
    cellData[c][at + 3] = value;
  }

  // Diagonalize the cell stress, pick the in-plane principal pair, and store
  // the Mises stress, the stress anisotropy and (optionally) the directions.
  // `byMagnitude` and `writeMises` are not options, they are the difference
  // between legacy's two stress passes: the direct one ranks the principal
  // values by magnitude and records the Mises stress, the neighbour-averaged
  // one ranks by signed value and records no Mises stress.
  void storeStressState(Matrix &cellData, size_t c, double stressCell[3][3],
                        const double normalGlob[3], double eps, bool storeDirs,
                        bool byMagnitude, bool writeMises) const {
    double eig[3][3];
    trbs::jacobiStressStyle(stressCell, eig, eps);
    int i1, i2, i3;
    if (byMagnitude)
      rankPrincipal(stressCell, i1, i2, i3);
    else
      rankPrincipalSigned(stressCell, i1, i2, i3);
    double v1 = stressCell[i1][i1];
    double v2 = stressCell[i2][i2];
    // If the largest principal direction is close to the cell normal it is
    // not an in-plane stress; step down to the next pair.
    if (std::fabs(normalGlob[0] * eig[0][i1] + normalGlob[1] * eig[1][i1] +
                  normalGlob[2] * eig[2][i1]) > 0.7) {
      i1 = i2;
      i2 = i3;
      v1 = v2;
      v2 = stressCell[i3][i3];
    }
    if (writeMises)
      cellData[c][variableIndex(0, kMisesStress)] =
          std::sqrt(v1 * v1 + v2 * v2 - v1 * v2);
    cellData[c][variableIndex(0, kStressAniso)] =
        std::fabs(v1) < eps
            ? 0.0
            : (parameter(6) == 0 ? 1 - std::fabs(v2 / v1)
                                 : (1 - std::fabs(v2 / v1)) * (v1 / parameter(6)));
    if (storeDirs && numVariableIndex(3) >= 1)
      storeDirection(cellData, c, variableIndex(3, 0), eig, i1, v1);
    if (storeDirs && numVariableIndex(3) == 2)
      storeDirection(cellData, c, variableIndex(3, 1), eig, i2, v2);
  }

  // With a positive neighbour weight the stress each cell reports is a
  // weighted mean of its own and its neighbours', which smooths the field the
  // microtubules respond to. Done in a second sweep so every cell reads its
  // neighbours' *own* stress rather than a partially updated one.
  void neighbourStressPass(Tissue &T, Matrix &cellData, double weight,
                           double eps, bool storeDirs) const {
    const size_t sti = variableIndex(0, kStressTensor);
    const size_t normalIndex = variableIndex(0, kNormalVector);
    for (size_t c = 0; c < T.numCell(); ++c) {
      double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      size_t counter = 0;
      for (size_t w : T.cell(c).walls) {
        const size_t nb = T.wall(w).otherCell(c);
        if (Tissue::isBackground(nb))
          continue;
        S[0][0] += weight * cellData[nb][sti];
        S[1][1] += weight * cellData[nb][sti + 1];
        S[2][2] += weight * cellData[nb][sti + 2];
        S[0][1] += weight * cellData[nb][sti + 3];
        S[2][0] += weight * cellData[nb][sti + 4];
        S[1][2] += weight * cellData[nb][sti + 5];
        ++counter;
      }
      if (counter != 0) {
        S[0][0] /= counter;
        S[1][1] /= counter;
        S[2][2] /= counter;
        S[0][1] /= counter;
        S[2][0] /= counter;
        S[1][2] /= counter;
      }
      S[0][0] += (1 - weight) * cellData[c][sti];
      S[1][1] += (1 - weight) * cellData[c][sti + 1];
      S[2][2] += (1 - weight) * cellData[c][sti + 2];
      S[0][1] += (1 - weight) * cellData[c][sti + 3];
      S[2][0] += (1 - weight) * cellData[c][sti + 4];
      S[1][2] += (1 - weight) * cellData[c][sti + 5];
      S[0][2] = S[2][0];
      S[1][0] = S[0][1];
      S[2][1] = S[1][2];
      const double normalGlob[3] = {cellData[c][normalIndex],
                                    cellData[c][normalIndex + 1],
                                    cellData[c][normalIndex + 2]};
      storeStressState(cellData, c, S, normalGlob, eps, storeDirs,
                       /*byMagnitude=*/false, /*writeMises=*/false);
    }
  }

  // Legacy reads the loosening Hill constants from parameters 9 and 10 - the
  // MT update flag and the double-resting-length flag - even though the
  // 13-parameter form exists precisely to carry them at 11 and 12, where
  // nothing reads them. Mode 10 is therefore unusable as written: the Hill
  // exponent can only be 0 or 1 and the half-max only 0-4, and neither of
  // those two flags can then be set independently. The 11-parameter form
  // keeps legacy's reading so existing models reproduce exactly; the
  // 13-parameter form uses the documented positions.
  double hillK() const { return numParameter() == 13 ? parameter(11) : parameter(9); }
  double hillN() const { return numParameter() == 13 ? parameter(12) : parameter(10); }

  // MT update flag 1 overwrites the direction from the TETA parameter. The
  // nested branch is a hard-coded four-cell test setup from legacy: unless
  // parameter 10 is exactly 100 it gives cells 0 and 2 the angle TETA and
  // cells 1 and 3 the angle parameter(10) - which is also the
  // double-resting-length flag, so the two collide.
  void setDirectionFromAngle(Matrix &cellData, size_t c) const {
    const size_t mt = variableIndex(0, kMTDirection);
    double angle = parameter(8);
    if (parameter(10) != 100 && (c == 1 || c == 3))
      angle = parameter(10);
    cellData[c][mt] = std::cos(angle);
    cellData[c][mt + 1] = std::sin(angle);
    cellData[c][mt + 2] = 0.0;
  }
};
TISSUE_REGISTER_REACTION(VertexFromTRBScenterTriangulationMT,
                         "VertexFromTRBScenterTriangulationMT")

} // namespace
} // namespace tissue
