#ifndef TISSUE_REACTIONS_TRBS_CORE_H
#define TISSUE_REACTIONS_TRBS_CORE_H
//
// Shared kernel for the triangular biquadratic spring (TRBS) reactions
// (legacy mechanicalTRBS.cc, and the center-triangulated growth rules in
// growth.cc that read TRBS strain).
//
// Legacy carries this kernel as a verbatim copy inside each variant -
// mechanicalTRBS.cc is 14.5k lines, of which the same ~200-line triangle
// routine appears seven times with a handful of lines changed. Everything
// here is that routine, factored so a variant supplies only what it actually
// varies: the material coefficients, the rest lengths, and whether it wants
// forces, stress, or both.
//
// Edge convention throughout: edge 0 joins nodes 0-1, edge 1 joins nodes 1-2,
// edge 2 joins nodes 0-2; cot[i] is the cotangent of the rest angle at node i.
// For a center-triangulated cell node 0 is the cell centre and nodes 1, 2 are
// consecutive vertices of the sorted cycle, so edge 1 is the wall.
//
// The expressions are kept literally as legacy writes them, including the two
// places where the same quantity is summed in a different order for different
// nodes: rewriting those to share a subexpression changes the last bit, and
// this kernel is validated by bit-exact comparison against legacy.
//
#include <algorithm>
#include <cmath>

namespace tissue {
namespace trbs {

// One triangle: positions, rest and current edge lengths, and everything
// derived from them.
struct Element {
  double pos[3][3];  // node positions (3D)
  double rest[3];    // rest edge lengths
  double cur[3];     // current edge lengths
  double restArea;
  double cot[3];     // cotangent of the rest angle at each node
  double delta[3];   // biquadratic strains, cur^2 - rest^2
};

// Rest angles enter only as cotangents and as the sine and cosine of one
// angle, and the law of cosines hands us the cosine directly - so the
// acos/tan round trip legacy performs is avoidable. With c = cos(theta),
//     cot(acos(c)) = c / sqrt(1 - c^2),  sin(acos(c)) = sqrt(1 - c^2)
// both exact. Profiling put a quarter of this kernel in tan and acos; the
// algebraic forms cost one sqrt each. The clamp keeps 1 - c^2 above ~2e-9, so
// the divisor cannot vanish.
inline double cosFromLengths(double la, double lb, double opposite) {
  const double c =
      (la * la + lb * lb - opposite * opposite) / (2.0 * la * lb);
  return std::max(-1.0 + 1e-9, std::min(1.0 - 1e-9, c));
}

inline double cotFromLengths(double la, double lb, double opposite) {
  const double c = cosFromLengths(la, lb, opposite);
  return c / std::sqrt(1.0 - c * c);
}

// Numerically stable Heron (edges sorted), guarded against degenerate or
// inverted trial configurations so a bad adaptive-solver trial step is
// rejected by error control instead of poisoning the state with NaNs.
inline double stableHeronArea(double e0, double e1, double e2) {
  double l[3] = {e0, e1, e2};
  std::sort(l, l + 3);
  const double a = l[2], b = l[1], c = l[0];
  const double heron =
      ((b + c) + a) * (-(a - b) + c) * ((a - b) + c) * ((b - c) + a);
  return 0.25 * std::sqrt(std::max(heron, 1e-12));
}

// Fills restArea, cot[] and delta[] from pos, rest and cur.
inline void completeElement(Element &e) {
  e.restArea = stableHeronArea(e.rest[0], e.rest[1], e.rest[2]);
  e.cot[0] = cotFromLengths(e.rest[0], e.rest[2], e.rest[1]);
  e.cot[1] = cotFromLengths(e.rest[0], e.rest[1], e.rest[2]);
  e.cot[2] = cotFromLengths(e.rest[1], e.rest[2], e.rest[0]);
  for (int i = 0; i < 3; ++i)
    e.delta[i] = e.cur[i] * e.cur[i] - e.rest[i] * e.rest[i];
}

// Euclidean distance between two of the element's nodes.
inline double nodeDistance(const Element &e, int a, int b) {
  double s = 0.0;
  for (int d = 0; d < 3; ++d) {
    const double diff = e.pos[a][d] - e.pos[b][d];
    s += diff * diff;
  }
  return std::sqrt(s);
}

struct Stiffness {
  double tensile[3];
  double angular[3];
};

// coefA multiplies the cotangent products and coefB the isotropic term.
// Isotropic TRBS passes (lambda + mio, mio); the transversely isotropic
// variants pass (lambdaT + 2 mioT, 2 mioT), which is the only difference
// between legacy's two stiffness blocks.
inline Stiffness stiffnessFrom(const Element &e, double coefA, double coefB) {
  const double t = 1.0 / (e.restArea * 16.0);
  const double *c = e.cot;
  return Stiffness{{(2 * c[2] * c[2] * coefA + coefB) * t,
                    (2 * c[0] * c[0] * coefA + coefB) * t,
                    (2 * c[1] * c[1] * coefA + coefB) * t},
                   {(2 * c[1] * c[2] * coefA - coefB) * t,
                    (2 * c[0] * c[2] * coefA - coefB) * t,
                    (2 * c[0] * c[1] * coefA - coefB) * t}};
}

// Node forces from the biquadratic strains (legacy mechanicalTRBS.cc
// 1089-1169). out[node][dim].
inline void elementForces(const Element &e, const Stiffness &s,
                          double out[3][3]) {
  const double *T = s.tensile, *A = s.angular, *D = e.delta;
  const double f0c = T[0] * D[0] + A[1] * D[1] + A[0] * D[2];
  const double f0v = T[2] * D[2] + A[2] * D[1] + A[0] * D[0];
  const double f1c = T[0] * D[0] + A[0] * D[2] + A[1] * D[1];
  const double f1v = T[1] * D[1] + A[2] * D[2] + A[1] * D[0];
  const double f2c = T[2] * D[2] + A[0] * D[0] + A[2] * D[1];
  const double f2v = T[1] * D[1] + A[1] * D[0] + A[2] * D[2];
  for (int d = 0; d < 3; ++d) {
    out[0][d] = f0c * (e.pos[1][d] - e.pos[0][d]) +
                f0v * (e.pos[2][d] - e.pos[0][d]);
    out[1][d] = f1c * (e.pos[0][d] - e.pos[1][d]) +
                f1v * (e.pos[2][d] - e.pos[1][d]);
    out[2][d] = f2c * (e.pos[0][d] - e.pos[2][d]) +
                f2v * (e.pos[1][d] - e.pos[2][d]);
  }
}

// Per-triangle Cauchy stress in the element plane, rotated to the global
// frame and accumulated area-weighted onto `stressGlobal` (legacy
// mechanicalTRBS.cc 843-1055).
inline void addCauchyStress(const Element &e, double lambda, double mio,
                            double stressGlobal[3][3]) {
  const double trE =
      (e.delta[1] * e.cot[0] + e.delta[2] * e.cot[1] + e.delta[0] * e.cot[2]) /
      (4.0 * e.restArea);

  // Both configurations laid out in the element's own 2D frame.
  const double cCur = cosFromLengths(e.cur[0], e.cur[1], e.cur[2]);
  const double Qa = cCur * e.cur[0];
  const double Qc = std::sqrt(1.0 - cCur * cCur) * e.cur[0];
  const double Qb = e.cur[1];
  const double cRest = cosFromLengths(e.rest[0], e.rest[1], e.rest[2]);
  const double Pa = cRest * e.rest[0];
  const double Pc = std::sqrt(1.0 - cRest * cRest) * e.rest[0];
  const double Pb = e.rest[1];
  const double shapeResting[3][2] = {
      {0.0, 1.0 / Pc}, {-1.0 / Pb, (Pa - Pb) / (Pb * Pc)}, {1.0 / Pb, -Pa / (Pb * Pc)}};
  const double posLocal[3][2] = {{Qa, Qc}, {0, 0}, {Qb, 0}};

  double F[2][2] = {{0, 0}, {0, 0}};
  for (int i = 0; i < 3; ++i) {
    F[0][0] += posLocal[i][0] * shapeResting[i][0];
    F[1][0] += posLocal[i][1] * shapeResting[i][0];
    F[0][1] += posLocal[i][0] * shapeResting[i][1];
    F[1][1] += posLocal[i][1] * shapeResting[i][1];
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

  const double curArea =
      0.25 * std::sqrt(std::max((e.cur[0] + e.cur[1] + e.cur[2]) *
                                    (-e.cur[0] + e.cur[1] + e.cur[2]) *
                                    (e.cur[0] - e.cur[1] + e.cur[2]) *
                                    (e.cur[0] + e.cur[1] - e.cur[2]),
                                1e-12));
  const double fac = curArea / e.restArea;
  double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  S[0][0] = fac * ((lambda * trE - mio / 2) * Bc[0][0] + (mio / 2) * B2[0][0]);
  S[1][0] = fac * ((lambda * trE - mio / 2) * Bc[1][0] + (mio / 2) * B2[1][0]);
  S[0][1] = fac * ((lambda * trE - mio / 2) * Bc[0][1] + (mio / 2) * B2[0][1]);
  S[1][1] = fac * ((lambda * trE - mio / 2) * Bc[1][1] + (mio / 2) * B2[1][1]);

  // Rotation local->global from the triangle's own frame.
  double X[3], Bv[3], Z[3], Yv[3];
  double tA = 0, tB = 0;
  for (int d = 0; d < 3; ++d) {
    X[d] = e.pos[2][d] - e.pos[1][d];
    Bv[d] = e.pos[0][d] - e.pos[1][d];
    tA += X[d] * X[d];
    tB += Bv[d] * Bv[d];
  }
  tA = std::sqrt(tA);
  tB = std::sqrt(tB);
  for (int d = 0; d < 3; ++d) {
    X[d] /= tA;
    Bv[d] /= tB;
  }
  Z[0] = X[1] * Bv[2] - X[2] * Bv[1];
  Z[1] = X[2] * Bv[0] - X[0] * Bv[2];
  Z[2] = X[0] * Bv[1] - X[1] * Bv[0];
  const double zn = std::sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
  for (int d = 0; d < 3; ++d)
    Z[d] /= zn;
  Yv[0] = Z[1] * X[2] - Z[2] * X[1];
  Yv[1] = Z[2] * X[0] - Z[0] * X[2];
  Yv[2] = Z[0] * X[1] - Z[1] * X[0];
  double R[3][3];
  for (int d = 0; d < 3; ++d) {
    R[d][0] = X[d];
    R[d][1] = Yv[d];
    R[d][2] = Z[d];
  }
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t) {
      double v = 0.0;
      for (int u = 0; u < 3; ++u)
        for (int w = 0; w < 3; ++w)
          v += R[r][u] * S[u][w] * R[t][w];
      stressGlobal[r][t] += e.restArea * v;
    }
}

// Jacobi diagonalization of a symmetric 3x3 (legacy mechanicalTRBS.cc
// 1269-1334). A is overwritten with the diagonal form and eig's columns hold
// the eigenvectors. A symmetric 3x3 needs only a handful of sweeps; the cap
// bounds the cost when an eigenvalue pair is degenerate.
inline void jacobiEigen3(double A[3][3], double eig[3][3]) {
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      eig[r][t] = (r == t) ? 1.0 : 0.0;
  double pivot = 1.0;
  const double pi = 3.1415;
  int iterations = 0;
  while (pivot > 0.00001 && ++iterations < 20) {
    int I = 1, J = 0;
    pivot = std::fabs(A[1][0]);
    if (std::fabs(A[2][0]) > pivot) {
      pivot = std::fabs(A[2][0]);
      I = 2;
      J = 0;
    }
    if (std::fabs(A[2][1]) > pivot) {
      pivot = std::fabs(A[2][1]);
      I = 2;
      J = 1;
    }
    double rotAngle;
    if (std::fabs(A[I][I] - A[J][J]) < 0.00001)
      rotAngle = pi / 4;
    else
      rotAngle = 0.5 * std::atan((2 * A[I][J]) / (A[J][J] - A[I][I]));
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
          tmp[r][t] += A[r][w] * rot[w][t];
      }
    for (int r = 0; r < 3; ++r)
      for (int t = 0; t < 3; ++t) {
        A[r][t] = 0.0;
        for (int w = 0; w < 3; ++w)
          A[r][t] += rot[w][r] * tmp[w][t];
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
}

// Principal direction, anisotropy a = 1 - |s2|/|s1|, and the signed s1, from
// a symmetric 3x3.
//
// The two in-plane principal values are picked by |eigenvalue|, not by signed
// eigenvalue. A membrane carries no load through its normal, so the
// shell-normal eigenvalue is the one nearest zero either way - but ranking by
// signed value only finds the in-plane pair while the sheet is in *tension*.
// Under compression all three in-plane values are negative, the ~0 normal
// sorts to the top, and the rule silently reports the normal as the principal
// direction with zero anisotropy and zero stress rather than failing. Ranking
// by magnitude gives identical answers in tension (where the in-plane values
// are positive and the normal is ~0) and the correct pair in compression.
struct Principal {
  double dir[3];
  double anisotropy;
  double s1;
};

inline Principal principalOf(double A[3][3]) {
  double eig[3][3];
  jacobiEigen3(A, eig);
  int order[3] = {0, 1, 2};
  const double ev[3] = {A[0][0], A[1][1], A[2][2]};
  for (int r = 0; r < 3; ++r)
    for (int t = r + 1; t < 3; ++t)
      if (std::fabs(ev[order[t]]) > std::fabs(ev[order[r]]))
        std::swap(order[r], order[t]);
  const double s1 = ev[order[0]];
  const double s2 = ev[order[1]];
  double a = 0.0;
  if (std::fabs(s1) > 1e-12)
    a = 1.0 - std::fabs(s2) / std::fabs(s1);
  a = std::max(0.0, std::min(1.0, a));
  Principal p{{eig[0][order[0]], eig[1][order[0]], eig[2][order[0]]}, a, s1};
  return p;
}

} // namespace trbs
} // namespace tissue

#endif
