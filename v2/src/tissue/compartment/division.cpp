//
// Cell division rules. Ported from legacy compartmentDivision.cc.
//
#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "tissue/compartment/compartment_change.h"
#include "tissue/core/random.h"
#include "tissue/core/tissue.h"

namespace tissue {
namespace {

constexpr double kPi = 3.14159265; // legacy myMath::pi()
inline int mathSign(double x) { return x >= 0 ? 1 : -1; }

// Divides cells above a volume threshold along a random direction through the
// center of mass (or a random interior point).
class DivisionVolumeRandomDirection : public CompartmentChange {
public:
  DivisionVolumeRandomDirection(const ParameterList &p, const IndexLevels &i) {
    configure("Division::VolumeRandomDirection", p, i, 4, {kAnyCount}, 1);
    // p: V_threshold, LWall_frac, Lwall_threshold, COM flag
  }

  int flag(Tissue &T, size_t i, Matrix &, Matrix &, Matrix &vertexData,
           Matrix &, Matrix &, Matrix &) override {
    if (T.cellVolume(i, vertexData) > parameter(0)) {
      std::cerr << "Cell " << i << " marked for division with volume "
                << T.cellVolume(i, vertexData) << std::endl;
      return 1;
    }
    return 0;
  }

  void update(Tissue &T, size_t i, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
              Matrix &vertexDerivs) override {
    const CellTopo &cell = T.cell(i);
    const size_t dimension = vertexData.cols();
    if (dimension != 2) {
      std::cerr << "Division::VolumeRandomDirection only supports 2D."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // Point on the division line.
    double com[2];
    if (parameter(3) == 1.0) {
      Vec3 c = T.cellPosition(i, vertexData);
      com[0] = c[0];
      com[1] = c[1];
    } else {
      auto r = T.randomPositionInCell(i, vertexData);
      if (!r)
        return; // legacy: silent no-division on sampling failure
      com[0] = (*r)[0];
      com[1] = (*r)[1];
    }
    // Random direction (legacy literal 2*3.14, preserved).
    double phi = 2.0 * 3.14 * random::Rnd();
    double n[2] = {std::sin(phi), std::cos(phi)};

    // Line/segment intersections with the cell's walls.
    size_t wI[2] = {cell.numWall(), cell.numWall()};
    double s[2] = {0.0, 0.0};
    size_t hits = 0;
    for (size_t k = 0; k < cell.numWall(); ++k) {
      const Wall &w = T.wall(cell.walls[k]);
      auto x1 = vertexData[w.vertex1];
      auto x2 = vertexData[w.vertex2];
      double w3[2] = {x2[0] - x1[0], x2[1] - x1[1]};
      double w0[2] = {com[0] - x1[0], com[1] - x1[1]};
      double a = n[0] * n[0] + n[1] * n[1];
      double b = n[0] * w3[0] + n[1] * w3[1];
      double c = w3[0] * w3[0] + w3[1] * w3[1];
      double d = n[0] * w0[0] + n[1] * w0[1];
      double e = w3[0] * w0[0] + w3[1] * w0[1];
      double fac = a * c - b * b;
      if (fac > 0.0) {
        double t = (a * e - b * d) / fac;
        if (t >= 0.0 && t < 1.0) {
          if (hits < 2) {
            wI[hits] = k;
            s[hits] = t;
          }
          ++hits;
        }
      }
    }
    if (hits != 2)
      return; // legacy: no division unless exactly two crossings

    Vec3 v1Pos{0, 0, 0}, v2Pos{0, 0, 0};
    {
      const Wall &w1 = T.wall(cell.walls[wI[0]]);
      const Wall &w2 = T.wall(cell.walls[wI[1]]);
      for (size_t d = 0; d < 2; ++d) {
        v1Pos[d] = vertexData[w1.vertex1][d] +
                   s[0] * (vertexData[w1.vertex2][d] - vertexData[w1.vertex1][d]);
        v2Pos[d] = vertexData[w2.vertex1][d] +
                   s[1] * (vertexData[w2.vertex2][d] - vertexData[w2.vertex1][d]);
      }
    }
    size_t numWallBefore = wallData.rows();
    T.divideCell(i, wI[0], wI[1], v1Pos, v2Pos, cellData, wallData, vertexData,
                 cellDerivs, wallDerivs, vertexDerivs, variableIndexLevel(0),
                 parameter(2));
    // New dividing wall's resting length scaled by LWall_frac.
    wallData[numWallBefore][0] *= parameter(1);
  }
};
TISSUE_REGISTER_COMPARTMENT_CHANGE(DivisionVolumeRandomDirection,
                                   "Division::VolumeRandomDirection",
                                   "DivisionVolumeRandomDirection")

// Divides cells above a volume threshold along the shortest straight cut
// through the center of mass (or a blend with a random interior point).
class DivisionShortestPath2D : public CompartmentChange {
public:
  DivisionShortestPath2D(const ParameterList &p, const IndexLevels &i) {
    // p: V_threshold, Lwall_fraction (unused, legacy quirk), Lwall_threshold,
    // CoM in [0,1]. Level 0: volume-dependent cell variable indices; optional
    // level 1: {age index} or {age, ownSize, sisterSize}.
    if (p.size() != 4)
      throw std::runtime_error(
          "Division::ShortestPath2D: uses four parameters (V_threshold, "
          "Lwall_fraction, Lwall_threshold, CoM).");
    if (p[3] < 0.0 || p[3] > 1.0)
      throw std::runtime_error(
          "Division::ShortestPath2D: CoM (parameter 3) must be in [0,1].");
    if ((i.size() == 2 && i[1].size() != 1 && i[1].size() != 3) ||
        (i.size() != 1 && i.size() != 2))
      throw std::runtime_error(
          "Division::ShortestPath2D: level 0 = volume-dependent indices; "
          "optional level 1 = {age index} or {age, size, sisterSize}.");
    configureLoose(p, i);
  }

  int flag(Tissue &T, size_t i, Matrix &, Matrix &, Matrix &vertexData,
           Matrix &, Matrix &, Matrix &) override {
    if (T.cellVolume(i, vertexData) > parameter(0)) {
      std::cerr << "Cell " << i << " marked for division with volume "
                << T.cellVolume(i, vertexData) << std::endl;
      return 1;
    }
    return 0;
  }

  struct Candidate {
    double distance;
    size_t wall1, wall2;
    double px, py, qx, qy;
  };

  void update(Tissue &T, size_t i, Matrix &cellData, Matrix &wallData,
              Matrix &vertexData, Matrix &cellDerivs, Matrix &wallDerivs,
              Matrix &vertexDerivs) override {
    if (vertexData.cols() != 2) {
      std::cerr << "Division::ShortestPath2D only supports two dimensions."
                << std::endl
                << "Consider using Division::ShortestPath." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::vector<Candidate> candidates = getCandidates(T, i, vertexData);
    if (candidates.empty()) {
      std::cerr << "Division::shortestPath2D.update() WARNING, cell " << i
                << " marked for division but no candidate shortest path found."
                << std::endl;
      return;
    }
    Candidate winner{DBL_MAX, 0, 0, 0, 0, 0, 0};
    for (const Candidate &c : candidates)
      if (c.distance < winner.distance)
        winner = c;

    Vec3 p{winner.px, winner.py, 0}, q{winner.qx, winner.qy, 0};
    if (numVariableIndexLevel() == 2) {
      size_t timeIndex = variableIndex(1, 0);
      std::cerr << "Cell age at division is " << cellData[i][timeIndex]
                << std::endl;
      cellData[i][timeIndex] = 0.0;
    }
    T.divideCell(i, winner.wall1, winner.wall2, p, q, cellData, wallData,
                 vertexData, cellDerivs, wallDerivs, vertexDerivs,
                 variableIndexLevel(0), parameter(2));

    if (numVariableIndexLevel() == 2 && numVariableIndex(1) == 3) {
      size_t sister = cellData.rows() - 1;
      double iSize = T.cellVolume(i, vertexData);
      double sisterSize = T.cellVolume(sister, vertexData);
      cellData[i][variableIndex(1, 1)] = iSize;
      cellData[i][variableIndex(1, 2)] = sisterSize;
      cellData[sister][variableIndex(1, 1)] = sisterSize;
      cellData[sister][variableIndex(1, 2)] = iSize;
    }
  }

private:
  // Validation already done in the constructor body; store members directly.
  void configureLoose(const ParameterList &p, const IndexLevels &i) {
    std::vector<size_t> counts;
    counts.push_back(kAnyCount);
    if (i.size() == 2)
      counts.push_back(kAnyCount);
    configure("Division::ShortestPath2D", p, i, 4, counts, 1);
  }

  static double f(double a, double sigma, double A, double B) {
    double sa = std::sin(a);
    double sb = std::sin(sigma - a);
    return -A * std::cos(a) / (sa * sa) +
           B * std::cos(kPi + sigma - a) / (sb * sb);
  }

  // Fixed 10-step bisection of f over (0, pi); returns 0 when no sign change.
  static double astar(double sigma, double A, double B) {
    double a = 0.0, b = kPi;
    double e = b - a;
    double u = f(a, sigma, A, B);
    double v = f(b, sigma, A, B);
    double c = 0.0;
    if (mathSign(u) == mathSign(v))
      return 0.0;
    for (int k = 0; k < 10; ++k) {
      e = 0.5 * e;
      c = a + e;
      double w = f(c, sigma, A, B);
      if (mathSign(w) != mathSign(u)) {
        b = c;
        v = w;
      } else {
        a = c;
        u = w;
      }
    }
    return c;
  }

  std::vector<Candidate> getCandidates(Tissue &T, size_t i,
                                       const Matrix &vertexData) {
    const CellTopo &cell = T.cell(i);
    Vec3 xc = T.cellPosition(i, vertexData);
    double ox = parameter(3) * xc[0];
    double oy = parameter(3) * xc[1];
    if (parameter(3) != 1.0) {
      auto r = T.randomPositionInCell(i, vertexData);
      if (!r)
        return {};
      ox += (1.0 - parameter(3)) * (*r)[0];
      oy += (1.0 - parameter(3)) * (*r)[1];
    }

    std::vector<Candidate> candidates;
    const size_t n = cell.numWall();
    for (size_t wi = 0; wi + 1 < n; ++wi) {
      for (size_t wj = wi + 1; wj < n; ++wj) {
        size_t wall1Index = wi, wall2Index = wj;
        // Orient the two edge vectors consistently relative to o; may swap
        // the two walls (legacy do/while loop). When o lies exactly on a
        // wall's supporting line the legacy loop cycles forever; v2 bounds it
        // and skips such degenerate pairs.
        double x1x, x1y, x2x, x2y, x1px, x1py, x2px, x2py;
        double vx, vy, ux, uy;
        bool flipped;
        int flipGuard = 0;
        do {
          if (++flipGuard > 8) {
            std::cerr << "Division::ShortestPath2D: skipping degenerate wall "
                         "pair (orientation does not converge)." << std::endl;
            break;
          }
          flipped = false;
          const Wall &w1 = T.wall(cell.walls[wall1Index]);
          const Wall &w2 = T.wall(cell.walls[wall2Index]);
          x1x = vertexData[w1.vertex1][0];
          x1y = vertexData[w1.vertex1][1];
          x2x = vertexData[w1.vertex2][0];
          x2y = vertexData[w1.vertex2][1];
          vx = x2x - x1x;
          vy = x2y - x1y;
          if (vx * (oy - x1y) - vy * (ox - x1x) > 0) {
            std::swap(x1x, x2x);
            std::swap(x1y, x2y);
            vx = -vx;
            vy = -vy;
          }
          x1px = vertexData[w2.vertex1][0];
          x1py = vertexData[w2.vertex1][1];
          x2px = vertexData[w2.vertex2][0];
          x2py = vertexData[w2.vertex2][1];
          ux = x2px - x1px;
          uy = x2py - x1py;
          if (ux * (oy - x1py) - uy * (ox - x1px) < 0) {
            std::swap(x1px, x2px);
            std::swap(x1py, x2py);
            ux = -ux;
            uy = -uy;
          }
          if (vx * uy - vy * ux > 0) {
            std::swap(wall1Index, wall2Index);
            flipped = true;
          }
        } while (flipped);
        if (flipGuard > 8)
          continue; // degenerate pair skipped

        double wx = ox - x1x, wy = oy - x1y;
        double wpx = ox - x1px, wpy = oy - x1py;
        double vv = vx * vx + vy * vy;
        double uu = ux * ux + uy * uy;
        double t = (vx * wx + vy * wy) / vv;
        double sPar = (ux * wpx + uy * wpy) / uu;
        double dvx = wx - t * vx, dvy = wy - t * vy;
        double dux = wpx - sPar * ux, duy = wpy - sPar * uy;
        double A = std::sqrt(dvx * dvx + dvy * dvy);
        double B = std::sqrt(dux * dux + duy * duy);
        double sigma = std::acos((vx * ux + vy * uy) /
                                 (std::sqrt(vv) * std::sqrt(uu)));
        double alpha = astar(sigma, A, B);
        double beta = kPi + sigma - alpha;

        double tp = t + (1.0 / std::sqrt(vv)) * A *
                            std::sin(alpha - 0.5 * kPi) / std::sin(alpha);
        double sp = sPar + (1.0 / std::sqrt(uu)) * B *
                               std::sin(beta - 0.5 * kPi) / std::sin(beta);
        if (tp <= 0.0 || tp >= 1.0 || sp <= 0.0 || sp >= 1.0)
          continue;
        double px = x1x + tp * vx, py = x1y + tp * vy;
        double qx = x1px + sp * ux, qy = x1py + sp * uy;
        double dx = qx - px, dy = qy - py;
        candidates.push_back(Candidate{std::sqrt(dx * dx + dy * dy),
                                       wall1Index, wall2Index, px, py, qx, qy});
      }
    }
    return candidates;
  }
};
TISSUE_REGISTER_COMPARTMENT_CHANGE(DivisionShortestPath2D,
                                   "Division::ShortestPath2D")

} // namespace
} // namespace tissue
