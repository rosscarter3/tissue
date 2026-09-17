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

#include "tissue/core/random.h"

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
  double cos[3];     // cosine of the same angles (the sliver test reads these)
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
  e.cos[0] = cosFromLengths(e.rest[0], e.rest[2], e.rest[1]);
  e.cos[1] = cosFromLengths(e.rest[0], e.rest[1], e.rest[2]);
  e.cos[2] = cosFromLengths(e.rest[1], e.rest[2], e.rest[0]);
  for (int i = 0; i < 3; ++i)
    e.cot[i] = e.cos[i] / std::sqrt(1.0 - e.cos[i] * e.cos[i]);
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


// The element in its own 2D frame: the rest shape vectors, the deformation
// gradient from rest to current, and the rotation that takes local vectors
// back to the global frame. Legacy rebuilds all of this inline in every
// variant - the isotropic stress path, each anisotropic force path, and again
// in the growth rules that read TRBS strain.
struct LocalFrame {
  double shapeResting[3][2];
  double F[2][2]; // deformation gradient, rest -> current
  double R[3][3]; // local -> global (columns are the element's x, y, z axes)
  double P[3];    // rest triangle in 2D: node 0 at (Pa, Pc), node 2 at (Pb, 0)
  double Q[3];    // the same for the current configuration
};

inline LocalFrame localFrameOf(const Element &e) {
  LocalFrame lf;
  // Both configurations laid out in 2D with node 1 at the origin and node 2
  // on the x axis.
  const double cCur = cosFromLengths(e.cur[0], e.cur[1], e.cur[2]);
  const double Qa = cCur * e.cur[0];
  const double Qc = std::sqrt(1.0 - cCur * cCur) * e.cur[0];
  const double Qb = e.cur[1];
  const double cRest = cosFromLengths(e.rest[0], e.rest[1], e.rest[2]);
  const double Pa = cRest * e.rest[0];
  const double Pc = std::sqrt(1.0 - cRest * cRest) * e.rest[0];
  const double Pb = e.rest[1];
  const double sv[3][2] = {{0.0, 1.0 / Pc},
                           {-1.0 / Pb, (Pa - Pb) / (Pb * Pc)},
                           {1.0 / Pb, -Pa / (Pb * Pc)}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 2; ++j)
      lf.shapeResting[i][j] = sv[i][j];

  const double posLocal[3][2] = {{Qa, Qc}, {0, 0}, {Qb, 0}};
  lf.F[0][0] = lf.F[1][0] = lf.F[0][1] = lf.F[1][1] = 0.0;
  for (int i = 0; i < 3; ++i) {
    lf.F[0][0] += posLocal[i][0] * sv[i][0];
    lf.F[1][0] += posLocal[i][1] * sv[i][0];
    lf.F[0][1] += posLocal[i][0] * sv[i][1];
    lf.F[1][1] += posLocal[i][1] * sv[i][1];
  }

  // Element frame from the current configuration: x along node1->node2, z
  // normal to the triangle, y completing the right-handed set.
  double X[3], B[3], Z[3], Y[3];
  double tA = 0, tB = 0;
  for (int d = 0; d < 3; ++d) {
    X[d] = e.pos[2][d] - e.pos[1][d];
    B[d] = e.pos[0][d] - e.pos[1][d];
    tA += X[d] * X[d];
    tB += B[d] * B[d];
  }
  tA = std::sqrt(tA);
  tB = std::sqrt(tB);
  for (int d = 0; d < 3; ++d) {
    X[d] /= tA;
    B[d] /= tB;
  }
  Z[0] = X[1] * B[2] - X[2] * B[1];
  Z[1] = X[2] * B[0] - X[0] * B[2];
  Z[2] = X[0] * B[1] - X[1] * B[0];
  const double zn = std::sqrt(Z[0] * Z[0] + Z[1] * Z[1] + Z[2] * Z[2]);
  for (int d = 0; d < 3; ++d)
    Z[d] /= zn;
  Y[0] = Z[1] * X[2] - Z[2] * X[1];
  Y[1] = Z[2] * X[0] - Z[0] * X[2];
  Y[2] = Z[0] * X[1] - Z[1] * X[0];
  for (int d = 0; d < 3; ++d) {
    lf.R[d][0] = X[d];
    lf.R[d][1] = Y[d];
    lf.R[d][2] = Z[d];
  }
  lf.P[0] = Pa;
  lf.P[1] = Pb;
  lf.P[2] = Pc;
  lf.Q[0] = Qa;
  lf.Q[1] = Qb;
  lf.Q[2] = Qc;
  return lf;
}

// Green-Lagrange strain E = (F^T F - I)/2 in the element frame.
inline void greenStrain(const LocalFrame &lf, double E[2][2]) {
  E[0][0] = 0.5 * (lf.F[0][0] * lf.F[0][0] + lf.F[1][0] * lf.F[1][0] - 1.0);
  E[0][1] = 0.5 * (lf.F[0][0] * lf.F[0][1] + lf.F[1][0] * lf.F[1][1]);
  E[1][0] = 0.5 * (lf.F[0][1] * lf.F[0][0] + lf.F[1][1] * lf.F[1][0]);
  E[1][1] = 0.5 * (lf.F[0][1] * lf.F[0][1] + lf.F[1][1] * lf.F[1][1] - 1.0);
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
  const LocalFrame lf = localFrameOf(e);

  double Bc[2][2]; // left Cauchy-Green B = F F^T
  Bc[0][0] = lf.F[0][0] * lf.F[0][0] + lf.F[0][1] * lf.F[0][1];
  Bc[1][0] = lf.F[1][0] * lf.F[0][0] + lf.F[1][1] * lf.F[0][1];
  Bc[0][1] = lf.F[0][0] * lf.F[1][0] + lf.F[0][1] * lf.F[1][1];
  Bc[1][1] = lf.F[1][0] * lf.F[1][0] + lf.F[1][1] * lf.F[1][1];
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

  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t) {
      double v = 0.0;
      for (int u = 0; u < 3; ++u)
        for (int w = 0; w < 3; ++w)
          v += lf.R[r][u] * S[u][w] * lf.R[t][w];
      stressGlobal[r][t] += e.restArea * v;
    }
}

// The transversely isotropic correction to the isotropic TRBS force.
//
// The isotropic part of the material is evaluated with the *transverse*
// moduli; this adds what the fibre direction contributes, through the extra
// Lame pair (dLambda, dMio) = (lambdaL - lambdaT, mioL - mioT). The fibre
// direction arrives in the global frame, is rotated into the element frame
// and pulled back to the rest configuration with the cofactor of F (so it
// follows the material rather than space), then normalized.
//
// A direction that pulls back to nothing - a fibre lying along the element
// normal - is replaced by a random one. Legacy draws that from libc `rand()`;
// this uses the seeded generator the rest of the rewrite shares, so a run
// stays reproducible. It is a degenerate fallback, and the validation runs do
// not reach it.
struct AnisoTerms {
  double trE;      // I1
  double I2, I4, I5;
  double aRest[2]; // fibre direction in the rest frame, normalized
  double aa[2][2]; // the dyad a (x) a
  double E[2][2];  // Green-Lagrange strain
  double Eaa[2][2], aaE[2][2];
  double deltaS[2][2]; // filled by whichever form the caller asks for
};

// Everything the fibre contributes that is common to the variants: the
// direction pulled back into the rest frame and the strain invariants built
// from it. The anisotropic stress deltaS is *not* set here - legacy writes it
// two different ways (see below) and each reaction is validated against its
// own.
inline AnisoTerms anisotropyInvariants(const LocalFrame &lf,
                                       const double dirGlobal[3]) {
  AnisoTerms out;
  double aLocal[3];
  for (int j = 0; j < 3; ++j)
    aLocal[j] = lf.R[0][j] * dirGlobal[0] + lf.R[1][j] * dirGlobal[1] +
                lf.R[2][j] * dirGlobal[2];
  double a[2] = {lf.F[1][1] * aLocal[0] - lf.F[0][1] * aLocal[1],
                 -lf.F[1][0] * aLocal[0] + lf.F[0][0] * aLocal[1]};
  const double measure = std::sqrt(a[0] * a[0] + a[1] * a[1]);
  if (measure < 0.0001) {
    const double angle = random::Rnd() * 2.0 * 3.14159265;
    a[0] = std::cos(angle);
    a[1] = std::sin(angle);
  } else {
    a[0] /= measure;
    a[1] /= measure;
  }
  out.aRest[0] = a[0];
  out.aRest[1] = a[1];
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t)
      out.aa[r][t] = a[r] * a[t];

  greenStrain(lf, out.E);
  const double(&E)[2][2] = out.E;
  out.trE = E[0][0] + E[1][1];
  const double E2[2][2] = {
      {E[0][0] * E[0][0] + E[0][1] * E[1][0],
       E[0][0] * E[0][1] + E[0][1] * E[1][1]},
      {E[1][0] * E[0][0] + E[1][1] * E[1][0],
       E[1][0] * E[0][1] + E[1][1] * E[1][1]}};
  out.I2 = E2[0][0] + E2[1][1];
  out.I5 = a[0] * a[0] * E2[0][0] + a[0] * a[1] * (E2[0][1] + E2[1][0]) +
           a[1] * a[1] * E2[1][1];
  out.I4 = a[0] * a[0] * E[0][0] + a[0] * a[1] * (E[0][1] + E[1][0]) +
           a[1] * a[1] * E[1][1];
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t) {
      out.Eaa[r][t] = E[r][0] * out.aa[0][t] + E[r][1] * out.aa[1][t];
      out.aaE[r][t] = out.aa[r][0] * E[0][t] + out.aa[r][1] * E[1][t];
    }
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t)
      out.deltaS[r][t] = 0.0;
  return out;
}


// The *other* way legacy pulls the fibre direction back to the rest frame,
// used by VertexFromTRBScenterTriangulationConcentrationHillMT. Instead of
// applying the cofactor of F to the direction, it takes the point one unit
// along the fibre from the element's centroid, expresses it in barycentric
// coordinates of the current triangle, maps those onto the rest triangle, and
// subtracts the rest centroid.
//
// The two agree for a rigid motion and differ under shear, but the more
// consequential difference is that this one keeps the *length* of the result:
// the caller scales the anisotropic Lame pair by it, so an element whose
// deformation shortens the fibre direction also weakens its anisotropy.
// `measure` is that length, returned before normalization.
struct BarycentricFibre {
  double aRest[2];
  double measure;
};

inline BarycentricFibre barycentricFibre(const LocalFrame &lf,
                                         const double dirGlobal[3]) {
  const double Pa = lf.P[0], Pb = lf.P[1], Pc = lf.P[2];
  const double Qa = lf.Q[0], Qb = lf.Q[1], Qc = lf.Q[2];
  double aLocal[3];
  for (int j = 0; j < 3; ++j)
    aLocal[j] = lf.R[0][j] * dirGlobal[0] + lf.R[1][j] * dirGlobal[1] +
                lf.R[2][j] * dirGlobal[2];
  const double cmCur[2] = {(Qa + Qb) / 3, Qc / 3};
  const double aCur[2] = {cmCur[0] + aLocal[0], cmCur[1] + aLocal[1]};
  const double svCur[3][3] = {{0, 1 / Qc, 0},
                              {-1 / Qb, (Qa - Qb) / (Qb * Qc), 1},
                              {1 / Qb, -Qa / (Qb * Qc), 0}};
  double bari[3];
  for (int i = 0; i < 3; ++i)
    bari[i] = svCur[i][0] * aCur[0] + svCur[i][1] * aCur[1] + svCur[i][2];
  const double aRestPoint[2] = {Pa * bari[0] + Pb * bari[2], Pc * bari[0]};
  BarycentricFibre out;
  out.aRest[0] = aRestPoint[0] - (Pa + Pb) / 3;
  out.aRest[1] = aRestPoint[1] - Pc / 3;
  out.measure = std::sqrt(out.aRest[0] * out.aRest[0] +
                          out.aRest[1] * out.aRest[1]);
  if (out.measure < 0.001) { // degenerate: pick a direction at random
    const double angle = random::Rnd() * 2.0 * 3.14159265;
    out.aRest[0] = std::cos(angle);
    out.aRest[1] = std::sin(angle);
  } else {
    out.aRest[0] /= out.measure;
    out.aRest[1] /= out.measure;
  }
  return out;
}


// A third way legacy pulls the fibre back to the rest frame, in
// VertexFromTRBSMT: F^T a, rather than the cofactor of F
// (VertexFromTRBScenterTriangulationMT) or the barycentric map
// (VertexFromTRBScenterTriangulationConcentrationHillMT). Legacy normalizes
// the result twice - the second pass is a no-op, and it makes that variant's
// degenerate-direction fallback unreachable, so there is none here.
inline void fibrePullbackTransposeF(const LocalFrame &lf,
                                    const double dirGlobal[3],
                                    double aRest[2]) {
  double aLocal[3];
  for (int j = 0; j < 3; ++j)
    aLocal[j] = lf.R[0][j] * dirGlobal[0] + lf.R[1][j] * dirGlobal[1] +
                lf.R[2][j] * dirGlobal[2];
  aRest[0] = lf.F[0][0] * aLocal[0] + lf.F[1][0] * aLocal[1];
  aRest[1] = lf.F[0][1] * aLocal[0] + lf.F[1][1] * aLocal[1];
  const double n = std::sqrt(aRest[0] * aRest[0] + aRest[1] * aRest[1]);
  aRest[0] /= n;
  aRest[1] /= n;
}

// Build the strain invariants from a fibre direction already in the rest
// frame (the barycentric pullback above supplies one).
inline AnisoTerms anisotropyInvariantsFrom(const LocalFrame &lf,
                                           const double aRest[2]) {
  AnisoTerms out;
  out.aRest[0] = aRest[0];
  out.aRest[1] = aRest[1];
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t)
      out.aa[r][t] = aRest[r] * aRest[t];
  greenStrain(lf, out.E);
  const double(&E)[2][2] = out.E;
  out.trE = E[0][0] + E[1][1];
  const double E2[2][2] = {
      {E[0][0] * E[0][0] + E[0][1] * E[1][0],
       E[0][0] * E[0][1] + E[0][1] * E[1][1]},
      {E[1][0] * E[0][0] + E[1][1] * E[1][0],
       E[1][0] * E[0][1] + E[1][1] * E[1][1]}};
  out.I2 = E2[0][0] + E2[1][1];
  out.I5 = aRest[0] * aRest[0] * E2[0][0] +
           aRest[0] * aRest[1] * (E2[0][1] + E2[1][0]) +
           aRest[1] * aRest[1] * E2[1][1];
  out.I4 = aRest[0] * aRest[0] * E[0][0] +
           aRest[0] * aRest[1] * (E[0][1] + E[1][0]) +
           aRest[1] * aRest[1] * E[1][1];
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t) {
      out.Eaa[r][t] = E[r][0] * out.aa[0][t] + E[r][1] * out.aa[1][t];
      out.aaE[r][t] = out.aa[r][0] * E[0][t] + out.aa[r][1] * E[1][t];
    }
  for (int r = 0; r < 2; ++r)
    for (int t = 0; t < 2; ++t)
      out.deltaS[r][t] = 0.0;
  return out;
}

// The equipartitioned anisotropic stress, as VertexFromTRBScenterTriangulationMT
// writes it. The atEa term sits on the diagonal only.
inline void equipartitionedDeltaS(AnisoTerms &t, double dLambda, double dMio) {
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 2; ++c)
      t.deltaS[r][c] = (dLambda / 2) * (t.trE * t.aa[r][c]) +
                       dMio * (t.Eaa[r][c] + t.aaE[r][c]);
  t.deltaS[0][0] += (dLambda / 2) * t.I4;
  t.deltaS[1][1] += (dLambda / 2) * t.I4;
}

// The form VertexFromTRBScenterTriangulationConcentrationHillMT uses instead:
// the shear term is halved and a -(dLambda + dMio) atEa a(x)a term subtracted.
// The two are not the same material law; each is kept with the reaction it
// belongs to.
inline void hillVariantDeltaS(AnisoTerms &t, double dLambda, double dMio) {
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 2; ++c)
      t.deltaS[r][c] = dLambda * (t.trE * t.aa[r][c]) +
                       (dMio / 2) * (t.Eaa[r][c] + t.aaE[r][c]) -
                       (dLambda + dMio) * t.I4 * t.aa[r][c];
  t.deltaS[0][0] += dLambda * t.I4;
  t.deltaS[1][1] += dLambda * t.I4;
}

// Push an anisotropic stress through to the three node forces: the 2nd
// Piola-Kirchhoff tensor, contracted with the rest shape vectors and rotated
// back to the global frame, negated.
inline void pushDeltaS(const Element &e, const LocalFrame &lf,
                       const double deltaS[2][2], double deltaF[3][3]) {
  double TPK[2][2];
  for (int r = 0; r < 2; ++r)
    for (int c = 0; c < 2; ++c)
      TPK[r][c] =
          e.restArea * (lf.F[r][0] * deltaS[0][c] + lf.F[r][1] * deltaS[1][c]);
  for (int i = 0; i < 3; ++i) {
    const double lx =
        TPK[0][0] * lf.shapeResting[i][0] + TPK[0][1] * lf.shapeResting[i][1];
    const double ly =
        TPK[1][0] * lf.shapeResting[i][0] + TPK[1][1] * lf.shapeResting[i][1];
    for (int d = 0; d < 3; ++d)
      deltaF[i][d] = -(lf.R[d][0] * lx + lf.R[d][1] * ly);
  }
}

// The equipartitioned correction end to end, for callers that want it in one
// step.
inline AnisoTerms anisotropicDeltaForce(const Element &e, const LocalFrame &lf,
                                        double dLambda, double dMio,
                                        const double dirGlobal[3],
                                        double deltaF[3][3]) {
  AnisoTerms t = anisotropyInvariants(lf, dirGlobal);
  equipartitionedDeltaS(t, dLambda, dMio);
  pushDeltaS(e, lf, t.deltaS, deltaF);
  return t;
}


// The second anisotropic force legacy uses, in
// VertexFromTRBScenterTriangulationConcentrationHillMT. Rather than pushing a
// stress tensor through the shape vectors, it differentiates the invariants
// I1, I4 and I5 directly with respect to the node positions and assembles
//   dF = -[ dLambda (I4 dI1 + I1 dI4) + dMio dI5
//           - (dMio + dLambda) I4 dI4 ] A_rest
//
// Transcribed as legacy writes it, including two things that look like slips
// but change the answer, so are not "corrected" here: the rest and current
// edge arrays are cyclically rotated just before this block while `cotan` and
// `Delta` are not, and the derIprim1 accumulation multiplies by
// position[m] inside a loop over i (so the i-sum only scales one position
// vector).
inline void invariantAnisotropicForce(const Element &e, const LocalFrame &lf,
                                      double area, const double aRest[2],
                                      double dLambda, double dMio,
                                      double deltaF[3][3]) {
  // Angle between the fibre and each rest shape vector.
  double teta[3];
  for (int i = 0; i < 3; ++i) {
    const double sx = lf.shapeResting[i][0], sy = lf.shapeResting[i][1];
    teta[i] = std::acos((sx * aRest[0] + sy * aRest[1]) /
                        std::sqrt(sx * sx + sy * sy + 0.0000001));
  }
  // Legacy's cyclic rotation of the edge arrays.
  const double rest[3] = {e.rest[1], e.rest[2], e.rest[0]};
  const double len[3] = {e.cur[1], e.cur[2], e.cur[0]};
  const double A = e.restArea;
  const double Rcirc2 = (0.25 * len[0] * len[1] * len[2] / area) *
                        (0.25 * len[0] * len[1] * len[2] / area);

  // The index not in {a, b}; only read when a != b.
  auto other = [](int a, int b) { return 3 - a - b; };
  auto DD = [&](int a, int b) {
    return a == b ? 0.25 * rest[a] * rest[a] / (A * A)
                  : -0.5 * e.cot[other(a, b)] / A;
  };
  auto aD = [&](int i) { return 0.5 * std::cos(teta[i]) * rest[i] / A; };

  double derIprim1[3][3], derIprim4[3][3], derIprim5[3][3];
  for (int m = 0; m < 3; ++m) {
    for (int coor = 0; coor < 3; ++coor)
      derIprim1[m][coor] = 0;
    for (int i = 0; i < 3; ++i) {
      const double DiDm = DD(i, m);
      for (int coor = 0; coor < 3; ++coor)
        derIprim1[m][coor] += 2 * DiDm * e.pos[m][coor];
    }
  }
  double Iprim4 = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double QiQj =
          (i == j) ? Rcirc2 : Rcirc2 - len[other(i, j)] * len[other(i, j)] * 0.5;
      Iprim4 += QiQj * aD(i) * aD(j);
    }
  for (int p = 0; p < 3; ++p) {
    for (int coor = 0; coor < 3; ++coor)
      derIprim4[p][coor] = 0;
    for (int m = 0; m < 3; ++m)
      for (int coor = 0; coor < 3; ++coor)
        derIprim4[p][coor] += aD(m) * e.pos[m][coor];
    for (int coor = 0; coor < 3; ++coor)
      derIprim4[p][coor] *= 2 * aD(p);
  }
  for (int p = 0; p < 3; ++p) {
    for (int coor = 0; coor < 3; ++coor)
      derIprim5[p][coor] = 0;
    for (int n = 0; n < 3; ++n)
      for (int r = 0; r < 3; ++r)
        for (int sIdx = 0; sIdx < 3; ++sIdx) {
          const double QrQs = e.pos[r][0] * e.pos[sIdx][0] +
                              e.pos[r][1] * e.pos[sIdx][1] +
                              e.pos[r][2] * e.pos[sIdx][2];
          const double DnDr = DD(n, r);
          const double DsDp = DD(sIdx, p);
          const double w =
              2 * (DnDr * aD(sIdx) * aD(p) + DsDp * aD(r) * aD(n)) * QrQs;
          for (int coor = 0; coor < 3; ++coor)
            derIprim5[p][coor] += w * e.pos[n][coor];
        }
  }

  const double I1 =
      (e.delta[1] * e.cot[0] + e.delta[2] * e.cot[1] + e.delta[0] * e.cot[2]) /
      (4 * A);
  const double I4 = 0.5 * Iprim4 - 0.5;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const double derI1 = 0.5 * derIprim1[i][j];
      const double derI4 = 0.5 * derIprim4[i][j];
      const double derI5 = 0.25 * derIprim5[i][j] - 0.5 * derIprim4[i][j];
      deltaF[i][j] = (-dLambda * (I4 * derI1 + I1 * derI4) - dMio * derI5 +
                      (dMio + dLambda) * I4 * derI4) *
                     A;
    }
}

// Jacobi diagonalization of a symmetric 3x3. A is overwritten with the
// diagonal form and eig's columns hold the eigenvectors.
//
// Legacy writes this loop four times, and no two are the same. The
// differences are small, none of them look deliberate, and all of them change
// the last digits - so rather than keep four near-identical copies they are
// one function and four named styles, which at least makes the differences
// legible:
//
//   - the rotation angle comes either from 0.5 atan(2 A_IJ/(A_JJ - A_II)) or
//     from half-angle formulas on the same quantity;
//   - the angle used when the two diagonal entries are within eps is written
//     as pi/4 with pi spelled 3.1415 in one place and 3.14159265 in another;
//   - the pivot is found either at the top of the loop (so one extra rotation
//     always runs, the last one below threshold) or before it and refreshed
//     at the bottom (so it stops without that rotation);
//   - the eigenvector columns are renormalized afterwards, or not.
//
// Legacy bounds none of these loops, so a degenerate tensor hangs the run.
// All four are capped here.
struct JacobiStyle {
  double eps;
  bool halfAngle;         // half-angle formulas instead of 0.5*atan
  double degenerateAngle; // used when |A_II - A_JJ| < eps (atan style only)
  bool pivotFirst;        // find the pivot at the top of the loop
  bool normalize;         // renormalize the eigenvector columns at the end
};

// The isotropic TRBS stress path (mechanicalTRBS.cc:1269).
inline constexpr JacobiStyle kJacobiIsotropic{1e-5, false, 3.1415 / 4, true,
                                              false};
// VertexFromTRBScenterTriangulationMT's update, stress pass.
inline constexpr JacobiStyle kJacobiMtStress{1e-6, true, 0.0, true, true};
// ... and its strain pass.
inline constexpr JacobiStyle kJacobiMtStrain{1e-6, true, 0.0, false, true};
// VertexFromTRBSMT's strain pass, and ...
inline constexpr JacobiStyle kJacobiMtPlain{1e-6, false, 3.14159265 / 4, false,
                                            false};
// ... its stress pass, which takes the pivot at the top instead.
inline constexpr JacobiStyle kJacobiMtPlainTop{1e-6, false, 3.14159265 / 4,
                                               true, false};

namespace detail {

inline bool jacobiPivot(const double A[3][3], double eps, int &I, int &J,
                        double &pivot) {
  pivot = std::fabs(A[1][0]);
  I = 1;
  J = 0;
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
  return pivot > eps;
}

inline void jacobiRotate(double A[3][3], double eig[3][3], int I, int J,
                         const JacobiStyle &st) {
  double Si, Co;
  const bool degenerate = std::fabs(A[I][I] - A[J][J]) < st.eps;
  if (st.halfAngle) {
    const double root2 = 0.70710678118;
    if (degenerate) {
      Si = root2;
      Co = root2;
    } else {
      double t = (2 * A[I][J]) / (A[J][J] - A[I][I]);
      t = 1 / std::sqrt(1 + t * t);
      Si = t < 1 ? root2 * std::sqrt(1 - t) : root2 * std::sqrt(t - 1);
      Co = root2 * std::sqrt(1 + t);
    }
  } else {
    const double angle =
        degenerate ? st.degenerateAngle
                   : 0.5 * std::atan((2 * A[I][J]) / (A[J][J] - A[I][I]));
    Si = std::sin(angle);
    Co = std::cos(angle);
  }
  double rot[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  rot[I][I] = Co;
  rot[J][J] = Co;
  rot[I][J] = Si;
  rot[J][I] = -Si;
  double tmp[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      for (int w = 0; w < 3; ++w)
        tmp[r][t] += A[r][w] * rot[w][t];
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      A[r][t] = 0;
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      for (int w = 0; w < 3; ++w)
        A[r][t] += rot[w][r] * tmp[w][t];
  double old[3][3];
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      old[r][t] = eig[r][t];
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t) {
      eig[r][t] = 0;
      for (int w = 0; w < 3; ++w)
        eig[r][t] += old[r][w] * rot[w][t];
    }
}

} // namespace detail

inline void jacobiEigen3(double A[3][3], double eig[3][3],
                         const JacobiStyle &st) {
  for (int r = 0; r < 3; ++r)
    for (int t = 0; t < 3; ++t)
      eig[r][t] = (r == t) ? 1.0 : 0.0;
  int I = 1, J = 0, sweeps = 0;
  double pivot = 1.0;
  if (st.pivotFirst) {
    while (pivot > st.eps && ++sweeps < 50) {
      detail::jacobiPivot(A, st.eps, I, J, pivot);
      detail::jacobiRotate(A, eig, I, J, st);
    }
  } else {
    bool go = detail::jacobiPivot(A, st.eps, I, J, pivot);
    while (go && ++sweeps < 50) {
      detail::jacobiRotate(A, eig, I, J, st);
      go = detail::jacobiPivot(A, st.eps, I, J, pivot);
    }
  }
  if (st.normalize)
    for (int col = 0; col < 3; ++col) {
      const double nrm = std::sqrt(eig[0][col] * eig[0][col] +
                                   eig[1][col] * eig[1][col] +
                                   eig[2][col] * eig[2][col]);
      if (nrm > 0)
        for (int r = 0; r < 3; ++r)
          eig[r][col] /= nrm;
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
  jacobiEigen3(A, eig, kJacobiIsotropic);
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
