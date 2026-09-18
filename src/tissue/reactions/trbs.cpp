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
#include <string>

#include "tissue/core/tissue.h"
#include "tissue/parallel/scatter.h"
#include "tissue/parallel/thread_pool.h"
#include "tissue/reactions/reaction.h"
#include "tissue/reactions/trbs_core.h"

namespace tissue {
namespace {

// Shared body of the center-triangulated TRBS reactions. Legacy repeats this
// loop verbatim in each variant; the only thing they actually vary is where
// Young's modulus comes from, so that arrives as a per-cell callback.
//
// Forces go to the cell centre (a cell variable, since the centre is not a
// mesh vertex) and to the two cycle vertices. When storeStress is set the
// cell's area-averaged Cauchy stress is diagonalized and written as
// [dir x, dir y, dir z, anisotropy, sigma1].
template <class YoungFn>
void ctTrbsEvaluate(Tissue &T, Matrix &cellData, Matrix &wallData,
                    Matrix &vertexData, Matrix &cellDerivs,
                    Matrix &vertexDerivs, size_t wallLengthIndex,
                    size_t comIndex, double poisson, YoungFn youngOf,
                    bool wantForces, bool storeStress, size_t stressIndex,
                    const char *who) {
  if (vertexData.cols() != 3)
    throw std::runtime_error(std::string(who) + " requires a 3D tissue.");
  const size_t lengthInternalIndex = comIndex + 3;

  parallelScatter1(
      T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
        for (size_t c = b; c < e; ++c) {
          const CellTopo &cell = T.cell(c);
          const size_t n = cell.numWall();
          if (cell.numVertex() != n)
            throw std::runtime_error(
                std::string(who) +
                ": needs the same number of vertices and walls per cell.");
          const double young = youngOf(c);
          const double lambda = young * poisson / (1 - poisson * poisson);
          const double mio = young / (1 + poisson);
          double stressCellGlobal[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
          double totalRestingArea = 0.0;
          for (size_t k = 0; k < n; ++k) {
            const size_t k1 = (k + 1) % n;
            const size_t v2 = cell.vertices[k];
            const size_t v3 = cell.vertices[k1];

            // Triangle nodes: 0 = cell centre, 1 = vertex(k), 2 = vertex(k+1).
            trbs::Element el;
            for (size_t d = 0; d < 3; ++d) {
              el.pos[0][d] = cellData[c][comIndex + d];
              el.pos[1][d] = vertexData[v2][d];
              el.pos[2][d] = vertexData[v3][d];
            }
            el.rest[0] = cellData[c][lengthInternalIndex + k];
            el.rest[1] = wallData[cell.walls[k]][wallLengthIndex];
            el.rest[2] = cellData[c][lengthInternalIndex + k1];
            el.cur[0] = trbs::nodeDistance(el, 0, 1);
            el.cur[1] = trbs::nodeDistance(el, 1, 2);
            el.cur[2] = trbs::nodeDistance(el, 0, 2);
            trbs::completeElement(el);

            if (wantForces) {
              double f[3][3];
              trbs::elementForces(
                  el, trbs::stiffnessFrom(el, lambda + mio, mio), f);
              for (size_t d = 0; d < 3; ++d) {
                cellDerivs[c][comIndex + d] += f[0][d];
                vOut[v2][d] += f[1][d];
                vOut[v3][d] += f[2][d];
              }
            }
            if (storeStress) {
              trbs::addCauchyStress(el, lambda, mio, stressCellGlobal);
              totalRestingArea += el.restArea;
            }
          }

          if (storeStress && totalRestingArea > 0.0) {
            for (int r = 0; r < 3; ++r)
              for (int t = 0; t < 3; ++t)
                stressCellGlobal[r][t] /= totalRestingArea;
            const trbs::Principal pr = trbs::principalOf(stressCellGlobal);
            for (size_t d = 0; d < 3; ++d)
              cellData[c][stressIndex + d] = pr.dir[d];
            cellData[c][stressIndex + 3] = pr.anisotropy;
            cellData[c][stressIndex + 4] = pr.s1;
          }
        }
      });
}

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
    const double young = parameter(0);
    ctTrbsEvaluate(T, cellData, wallData, vertexData, cellDerivs, vertexDerivs,
                   variableIndex(0, 0), variableIndex(1, 0), parameter(1),
                   [&](size_t) { return young; }, wantForces,
                   wantStress && numVariableIndexLevel() == 3,
                   numVariableIndexLevel() == 3 ? variableIndex(2, 0) : 0,
                   "VertexFromTRBScenterTriangulation");
  }
};
TISSUE_REGISTER_REACTION(VertexFromTRBScenterTriangulation,
                         "VertexFromTRBScenterTriangulation")

// Same material, but Young's modulus is set per cell by an *inhibitory* Hill
// function of a cell concentration:
//     Y = Y_min + Y_max K^n / (K^n + c^n)
// so the species softens the wall as it accumulates. Legacy ships this as a
// 250-line copy of the reaction above with that one line changed.
class VertexFromTRBScenterTriangulationConcentrationHill : public Reaction {
public:
  VertexFromTRBScenterTriangulationConcentrationHill(const ParameterList &p,
                                                     const IndexLevels &i) {
    const bool ok = (i.size() == 2 || i.size() == 3) && i[0].size() == 2 &&
                    i[1].size() == 1 && (i.size() == 2 || i[2].size() == 1);
    if (!ok)
      throw std::runtime_error(
          "VertexFromTRBScenterTriangulationConcentrationHill: wall length "
          "and concentration indices in level 0; start of the "
          "center-triangulation cell variables in level 1; optional level 2 = "
          "start index for storing the cell stress state (5 cell variables).");
    configure("VertexFromTRBScenterTriangulationConcentrationHill", p, i, 5,
              i.size() == 3 ? std::vector<size_t>{2, 1, 1}
                            : std::vector<size_t>{2, 1},
              {"Y_mod_min", "Y_mod_max", "P_ratio", "K_hill", "n_hill"});
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &,
              Matrix &vertexDerivs) override {
    evaluate(T, cellData, wallData, vertexData, cellDerivs, vertexDerivs, true,
             false);
  }
  void initiate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &cellDerivs, Matrix &,
                Matrix &vertexDerivs) override {
    if (numVariableIndexLevel() == 3)
      evaluate(T, cellData, wallData, vertexData, cellDerivs, vertexDerivs,
               false, true);
  }
  void positionalCellVariables(std::vector<size_t> &out) const override {
    const size_t com = variableIndex(1, 0);
    out.push_back(com);
    out.push_back(com + 1);
    out.push_back(com + 2);
  }
  void update(Tissue &T, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, double) override {
    if (numVariableIndexLevel() != 3)
      return;
    if (!scratchCell_.sameShape(cellData))
      scratchCell_.reshapeLike(cellData);
    if (!scratchVertex_.sameShape(vertexData))
      scratchVertex_.reshapeLike(vertexData);
    evaluate(T, cellData, wallData, vertexData, scratchCell_, scratchVertex_,
             false, true);
  }

private:
  Matrix scratchCell_, scratchVertex_;

  void evaluate(Tissue &T, Matrix &cellData, Matrix &wallData,
                Matrix &vertexData, Matrix &cellDerivs, Matrix &vertexDerivs,
                bool wantForces, bool wantStress) {
    const size_t concIndex = variableIndex(0, 1);
    const double kPow = std::pow(parameter(3), parameter(4));
    ctTrbsEvaluate(
        T, cellData, wallData, vertexData, cellDerivs, vertexDerivs,
        variableIndex(0, 0), variableIndex(1, 0), parameter(2),
        [&](size_t c) {
          return parameter(0) +
                 parameter(1) * kPow /
                     (kPow + std::pow(cellData[c][concIndex], parameter(4)));
        },
        wantForces, wantStress && numVariableIndexLevel() == 3,
        numVariableIndexLevel() == 3 ? variableIndex(2, 0) : 0,
        "VertexFromTRBScenterTriangulationConcentrationHill");
  }
};
TISSUE_REGISTER_REACTION(VertexFromTRBScenterTriangulationConcentrationHill,
                         "VertexFromTRBScenterTriangulationConcentrationHill")

// The same elasticity without center triangulation: the cell *is* one
// triangle, its three walls are the three edges, and there is no centre node.
// Only defined for triangular cells, as in legacy.
class VertexFromTRBS : public Reaction {
public:
  VertexFromTRBS(const ParameterList &p, const IndexLevels &i) {
    if (i.size() != 1 || i[0].size() != 1)
      throw std::runtime_error(
          "VertexFromTRBS: one index, the wall resting-length index.");
    configure("VertexFromTRBS", p, i, 2, {1}, {"Y_mod", "P_ratio"});
  }
  void derivs(Tissue &T, Matrix &, Matrix &wallData, Matrix &vertexData,
              Matrix &, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error("VertexFromTRBS requires a 3D tissue.");
    const size_t wallLengthIndex = variableIndex(0, 0);
    const double young = parameter(0), poisson = parameter(1);
    const double lambda = young * poisson / (1 - poisson * poisson);
    const double mio = young / (1 + poisson);
    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            if (cell.numWall() != 3)
              throw std::runtime_error(
                  "VertexFromTRBS: only defined for triangular cells.");
            trbs::Element el;
            for (int node = 0; node < 3; ++node)
              for (size_t d = 0; d < 3; ++d)
                el.pos[node][d] = vertexData[cell.vertices[node]][d];
            // Wall k joins vertex k and vertex k+1, which is exactly the
            // element's edge order (0-1, 1-2, then 0-2 for the closing wall).
            for (int k = 0; k < 3; ++k)
              el.rest[k] = wallData[cell.walls[k]][wallLengthIndex];
            el.cur[0] = trbs::nodeDistance(el, 0, 1);
            el.cur[1] = trbs::nodeDistance(el, 1, 2);
            el.cur[2] = trbs::nodeDistance(el, 0, 2);
            trbs::completeElement(el);
            double f[3][3];
            trbs::elementForces(el, trbs::stiffnessFrom(el, lambda + mio, mio),
                                f);
            for (int node = 0; node < 3; ++node)
              for (size_t d = 0; d < 3; ++d)
                vOut[cell.vertices[node]][d] += f[node][d];
          }
        });
  }
};
TISSUE_REGISTER_REACTION(VertexFromTRBS, "VertexFromTRBS")

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

// Pressure on a center-triangulated shell that *ramps in* over a set time
// rather than being applied at full strength from the start, which is how the
// published models inflate a tissue without kicking it. Per CT triangle the
// force is k_force * A * n_hat, distributed to all three nodes (not divided
// between them - each gets the whole thing, as legacy has it).
//
//   areaFlag 0: A = 1/3, i.e. no area weighting at all
//   areaFlag 1: A = Area/2
//   areaFlag 2: A = Area/2, and only the z component is applied, through a
//               second ramp
//
// With a fourth parameter the pressure ramps from p3 to p0 instead of 0 to p0,
// and the current value is written back into a cell variable named by a second
// index.
//
// Two legacy quirks kept. Under areaFlag 2 the *first* ramp is never advanced
// - `update` only steps timeFactor1 for flags 0 and 1 - so the main force term
// stays at zero for the whole run and only the z-only term does anything.
// And the ramp state is per-reaction rather than per-cell, so it is shared by
// every cell, which is the intent.
class Pressure3DCenterTriangulationLinear : public Reaction {
public:
  Pressure3DCenterTriangulationLinear(const ParameterList &p,
                                      const IndexLevels &i) {
    if (p.size() != 3 && p.size() != 4)
      throw std::runtime_error(
          "Pressure3D::CenterTriangulation::Linear: uses three or four "
          "parameters (k_force, areaFlag, deltaT, [k_force_start]).");
    if (p[1] != 0.0 && p[1] != 1.0 && p[1] != 2.0)
      throw std::runtime_error("Pressure3D::CenterTriangulation::Linear: "
                               "areaFlag must be 0, 1 or 2.");
    if (p[2] < 0.0)
      throw std::runtime_error("Pressure3D::CenterTriangulation::Linear: "
                               "deltaT must not be negative.");
    if (i.size() != 1 || i[0].empty())
      throw std::runtime_error(
          "Pressure3D::CenterTriangulation::Linear: one index level - the "
          "start of the center-triangulation cell variables, and with four "
          "parameters a second index to report the current pressure in.");
    configure("Pressure3D::CenterTriangulation::Linear", p, i, p.size(),
              {i[0].size()},
              p.size() == 4 ? std::vector<std::string>{"k_force", "areaFlag",
                                                       "deltaT", "k_force_0"}
                            : std::vector<std::string>{"k_force", "areaFlag",
                                                       "deltaT"});
  }
  void positionalCellVariables(std::vector<size_t> &out) const override {
    const size_t com = variableIndex(0, 0);
    out.push_back(com);
    out.push_back(com + 1);
    out.push_back(com + 2);
  }
  void update(Tissue &, Matrix &, Matrix &, Matrix &, double h) override {
    if (parameter(1) == 0.0 || parameter(1) == 1.0) {
      if (timeFactor1_ < 1.0)
        timeFactor1_ += h / parameter(2);
      if (timeFactor1_ > 1.0)
        timeFactor1_ = 1.0;
    }
    if (parameter(1) == 2.0) {
      if (timeFactor2_ < 1.0)
        timeFactor2_ += h / parameter(2);
      if (timeFactor2_ > 1.0)
        timeFactor2_ = 1.0;
    }
  }
  void derivs(Tissue &T, Matrix &cellData, Matrix &, Matrix &vertexData,
              Matrix &cellDerivs, Matrix &, Matrix &vertexDerivs) override {
    if (vertexData.cols() != 3)
      throw std::runtime_error(
          "Pressure3D::CenterTriangulation::Linear requires a 3D tissue.");
    const size_t comIndex = variableIndex(0, 0);
    const bool ramped = numParameter() == 4;
    const double pressure =
        ramped ? parameter(3) + timeFactor1_ * (parameter(0) - parameter(3))
               : timeFactor1_ * parameter(0);
    if (ramped && numVariableIndex(0) > 1)
      for (size_t c = 0; c < T.numCell(); ++c)
        cellData[c][variableIndex(0, 1)] = pressure;

    parallelScatter1(
        T.numCell(), vertexDerivs, [&](size_t b, size_t e, Matrix &vOut) {
          for (size_t c = b; c < e; ++c) {
            const CellTopo &cell = T.cell(c);
            const size_t n = cell.numWall();
            if (cell.numVertex() != n)
              throw std::runtime_error(
                  "Pressure3D::CenterTriangulation::Linear: needs the same "
                  "number of vertices and walls.");
            for (size_t k = 0; k < n; ++k) {
              const size_t v2 = cell.vertices[k];
              const size_t v3 = cell.vertices[(k + 1) % n];
              double pos[3][3];
              for (size_t d = 0; d < 3; ++d) {
                pos[0][d] = cellData[c][comIndex + d];
                pos[1][d] = vertexData[v2][d];
                pos[2][d] = vertexData[v3][d];
              }
              auto dist = [&](int a, int b2) {
                double s2 = 0;
                for (size_t d = 0; d < 3; ++d) {
                  const double diff = pos[a][d] - pos[b2][d];
                  s2 += diff * diff;
                }
                return std::sqrt(s2);
              };
              const double l0 = dist(0, 1), l1 = dist(1, 2), l2 = dist(0, 2);
              const double area =
                  0.25 * std::sqrt((l0 + l1 + l2) * (-l0 + l1 + l2) *
                                   (l0 - l1 + l2) * (l0 + l1 - l2));
              // Unit normal of the triangle, from its own frame.
              double X[3], B[3], nrm[3];
              double tA = 0, tB = 0;
              for (size_t d = 0; d < 3; ++d) {
                X[d] = pos[2][d] - pos[1][d];
                B[d] = pos[0][d] - pos[1][d];
                tA += X[d] * X[d];
                tB += B[d] * B[d];
              }
              tA = std::sqrt(tA);
              tB = std::sqrt(tB);
              for (size_t d = 0; d < 3; ++d) {
                X[d] /= tA;
                B[d] /= tB;
              }
              nrm[0] = X[1] * B[2] - X[2] * B[1];
              nrm[1] = X[2] * B[0] - X[0] * B[2];
              nrm[2] = X[0] * B[1] - X[1] * B[0];
              const double nn = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] +
                                          nrm[2] * nrm[2]);
              for (size_t d = 0; d < 3; ++d)
                nrm[d] /= nn;

              const double A = (parameter(1) == 1.0 || parameter(1) == 2.0)
                                   ? area / 2
                                   : 1.0 / 3;
              if (parameter(1) == 0.0 || parameter(1) == 1.0) {
                const double coeff = pressure * A;
                for (size_t d = 0; d < 3; ++d) {
                  cellDerivs[c][comIndex + d] += coeff * nrm[d];
                  vOut[v2][d] += coeff * nrm[d];
                  vOut[v3][d] += coeff * nrm[d];
                }
              }
              if (parameter(1) == 2.0) {
                const double coeff = timeFactor2_ * parameter(0) * A;
                cellDerivs[c][comIndex + 2] += coeff * nrm[2];
                vOut[v2][2] += coeff * nrm[2];
                vOut[v3][2] += coeff * nrm[2];
              }
            }
          }
        });
  }

private:
  double timeFactor1_ = 0.0, timeFactor2_ = 0.0;
};
TISSUE_REGISTER_REACTION(Pressure3DCenterTriangulationLinear,
                         "Pressure3D::CenterTriangulation::Linear",
                         "CenterTriangulation::Pressure3D::Linear",
                         "VertexFromCellPlaneLinearCenterTriangulation")

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
