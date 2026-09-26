#include "tissue/core/meshing.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

namespace tissue {
namespace {

using Pt = std::array<double, 2>;
using Tri = std::array<std::size_t, 3>;

double cross2(const Pt &o, const Pt &a, const Pt &b) {
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

double signedArea(const std::vector<Pt> &P) {
  double s = 0.0;
  const std::size_t n = P.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Pt &a = P[i], &b = P[(i + 1) % n];
    s += a[0] * b[1] - b[0] * a[1];
  }
  return 0.5 * s;
}

bool pointInTriangle(const Pt &p, const Pt &a, const Pt &b, const Pt &c) {
  const double d1 = cross2(a, b, p), d2 = cross2(b, c, p), d3 = cross2(c, a, p);
  const bool neg = d1 < 0 || d2 < 0 || d3 < 0;
  const bool pos = d1 > 0 || d2 > 0 || d3 > 0;
  return !(neg && pos);
}

bool pointInPolygon(const Pt &p, const std::vector<Pt> &P) {
  bool in = false;
  const std::size_t n = P.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Pt &a = P[i], &b = P[(i + 1) % n];
    if ((a[1] > p[1]) != (b[1] > p[1])) {
      const double t = (p[1] - a[1]) / (b[1] - a[1]);
      if (p[0] < a[0] + t * (b[0] - a[0]))
        in = !in;
    }
  }
  return in;
}

double segDistance(const Pt &p, const Pt &a, const Pt &b) {
  const double dx = b[0] - a[0], dy = b[1] - a[1];
  const double l2 = dx * dx + dy * dy;
  double t = 0.0;
  if (l2 > 0.0)
    t = std::min(1.0, std::max(0.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2));
  const double qx = a[0] + t * dx - p[0], qy = a[1] + t * dy - p[1];
  return std::sqrt(qx * qx + qy * qy);
}

// Is d strictly inside the circumcircle of triangle a,b,c?
bool inCircumcircle(Pt a, Pt b, const Pt &c, const Pt &d) {
  if (cross2(a, b, c) <= 0.0)
    std::swap(a, b);
  const double ax = a[0] - d[0], ay = a[1] - d[1];
  const double bx = b[0] - d[0], by = b[1] - d[1];
  const double cx = c[0] - d[0], cy = c[1] - d[1];
  return (ax * ax + ay * ay) * (bx * cy - by * cx) -
             (bx * bx + by * by) * (ax * cy - ay * cx) +
             (cx * cx + cy * cy) * (ax * by - ay * bx) >
         1e-12;
}

// Radius ratio: circumradius over twice the inradius. It is 1 for an
// equilateral triangle and cannot be less, so anything below 1 is a
// degenerate element and is reported as infinite rather than as a very good
// triangle.
//
// That floor is not pedantry. A triangle with two coincident corners has a
// cross product of order 1e-18 rather than exactly zero, so an `area <= 0`
// guard does not catch it; the computation then continues with one edge of
// length zero, the circumradius comes out as zero, and the triangle scores a
// quality of 0 -- better than equilateral. Every check built on this measure
// therefore passed it, including this file's own unit test, an offline
// validation against a Python reference, and CellMesh::Initiate's report of
// its worst element. 1064 degenerate triangles across 30 real cells went
// unnoticed that way.
double radiusRatio(const Pt &a, const Pt &b, const Pt &c) {
  const double a2 = cross2(a, b, c);
  auto len = [](const Pt &p, const Pt &q) {
    const double dx = q[0] - p[0], dy = q[1] - p[1];
    return std::sqrt(dx * dx + dy * dy);
  };
  const double la = len(b, c), lb = len(c, a), lc = len(a, b);
  const double inf = std::numeric_limits<double>::infinity();
  if (la <= 0.0 || lb <= 0.0 || lc <= 0.0)
    return inf;
  // Area has to be judged against the element's own size, not against zero.
  const double scale = std::max(la, std::max(lb, lc));
  if (a2 <= 1e-12 * scale * scale)
    return inf;
  const double A = 0.5 * a2;
  const double s = 0.5 * (la + lb + lc);
  if (s <= 0.0)
    return inf;
  const double R = la * lb * lc / (4.0 * A);
  const double r = A / s;
  if (!(r > 0.0))
    return inf;
  const double ratio = R / (2.0 * r);
  return ratio >= 1.0 ? ratio : inf;
}

// Ear clipping. Always succeeds for a simple polygon, convex or not, and
// uses no vertices but the polygon's own -- which is what makes it the safe
// starting point before any interior point is considered.
std::vector<Tri> earClip(const std::vector<Pt> &P) {
  const std::size_t n = P.size();
  std::vector<std::size_t> idx(n);
  for (std::size_t i = 0; i < n; ++i)
    idx[i] = i;
  std::vector<Tri> tris;
  std::size_t guard = 0;
  while (idx.size() > 3 && guard < 12 * n) {
    ++guard;
    bool clipped = false;
    for (std::size_t t = 0; t < idx.size(); ++t) {
      const std::size_t i0 = idx[(t + idx.size() - 1) % idx.size()];
      const std::size_t i1 = idx[t];
      const std::size_t i2 = idx[(t + 1) % idx.size()];
      if (cross2(P[i0], P[i1], P[i2]) <= 0.0)
        continue; // reflex, not an ear
      bool blocked = false;
      for (std::size_t j : idx) {
        if (j == i0 || j == i1 || j == i2)
          continue;
        if (pointInTriangle(P[j], P[i0], P[i1], P[i2])) {
          blocked = true;
          break;
        }
      }
      if (blocked)
        continue;
      tris.push_back({i0, i1, i2});
      idx.erase(idx.begin() + static_cast<long>(t));
      clipped = true;
      break;
    }
    if (!clipped)
      break; // numerically stuck; return what we have
  }
  if (idx.size() == 3)
    tris.push_back({idx[0], idx[1], idx[2]});
  return tris;
}

// Lawson flips: make the triangulation Delaunay wherever it is allowed to be.
//
// Boundary edges are never flipped -- they are the cell outline and are
// shared with a neighbour -- so what comes out is the constrained Delaunay
// triangulation, which for a fixed vertex set is the quality-optimal one.
//
// This keeps a queue of edges to test and pushes the neighbours of each edge
// it flips, which is the textbook algorithm and terminates. The first version
// here rebuilt the whole edge map and performed one flip per pass under a
// fixed budget of 200 passes. On a mesh with a thousand interior points that
// runs out long before the triangulation is Delaunay, and the symptom was
// quality stuck at 11 however finely the outline was resolved: refining
// bought more points and no more flips with which to organise them.
void delaunayFlips(const std::vector<Pt> &P, std::vector<Tri> &tris,
                   const std::set<std::pair<std::size_t, std::size_t>> &fixed) {
  auto key = [](std::size_t a, std::size_t b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
  };
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> adj;
  auto addTri = [&](std::size_t t) {
    const Tri &T = tris[t];
    adj[key(T[0], T[1])].push_back(t);
    adj[key(T[1], T[2])].push_back(t);
    adj[key(T[2], T[0])].push_back(t);
  };
  auto dropTri = [&](std::size_t t) {
    const Tri &T = tris[t];
    const std::pair<std::size_t, std::size_t> es[3] = {
        key(T[0], T[1]), key(T[1], T[2]), key(T[2], T[0])};
    for (const auto &e : es) {
      auto &v = adj[e];
      v.erase(std::remove(v.begin(), v.end(), t), v.end());
    }
  };
  for (std::size_t t = 0; t < tris.size(); ++t)
    addTri(t);

  std::vector<std::pair<std::size_t, std::size_t>> queue;
  for (const auto &kv : adj)
    if (!fixed.count(kv.first))
      queue.push_back(kv.first);

  std::size_t guard = 0;
  const std::size_t budget = 200 * tris.size() + 1000;
  while (!queue.empty() && guard++ < budget) {
    const auto e = queue.back();
    queue.pop_back();
    auto it = adj.find(e);
    if (it == adj.end() || it->second.size() != 2 || fixed.count(e))
      continue;
    const std::size_t t1 = it->second[0], t2 = it->second[1];
    const std::size_t u = e.first, v = e.second;
    std::size_t p = SIZE_MAX, q = SIZE_MAX;
    for (std::size_t x : tris[t1])
      if (x != u && x != v)
        p = x;
    for (std::size_t x : tris[t2])
      if (x != u && x != v)
        q = x;
    if (p == SIZE_MAX || q == SIZE_MAX || p == q)
      continue;
    if (!inCircumcircle(P[p], P[u], P[v], P[q]))
      continue;
    if (cross2(P[p], P[u], P[q]) <= 0.0 || cross2(P[p], P[q], P[v]) <= 0.0)
      continue; // the quadrilateral is not convex; the flip would invert
    dropTri(t1);
    dropTri(t2);
    tris[t1] = {p, u, q};
    tris[t2] = {p, q, v};
    addTri(t1);
    addTri(t2);
    const std::pair<std::size_t, std::size_t> ns[4] = {
        key(p, u), key(u, q), key(q, v), key(v, p)};
    for (const auto &ne : ns)
      if (!fixed.count(ne))
        queue.push_back(ne);
  }
}

// Delaunay triangulation of a point set by incremental insertion.
//
// Used in place of ear clipping the outline and then inserting the interior
// points one at a time into it. That order leaves slivers the insertion
// cannot reach: on a segmented outline, which has long runs of nearly
// collinear vertices, ear clipping bridges across such a run and the
// resulting sliver has all three corners on the boundary, so no interior
// point is adjacent to it and no flip can remove it -- its quadrilateral is
// not convex. Measured on a real lobed cell, that left a worst element of
// 957 where triangulating the whole point set at once gives 2.8.
std::vector<Tri> bowyerWatson(const std::vector<Pt> &P, std::size_t nReal) {
  double lox = P[0][0], hix = P[0][0], loy = P[0][1], hiy = P[0][1];
  for (std::size_t i = 0; i < nReal; ++i) {
    lox = std::min(lox, P[i][0]); hix = std::max(hix, P[i][0]);
    loy = std::min(loy, P[i][1]); hiy = std::max(hiy, P[i][1]);
  }
  const double cx = 0.5 * (lox + hix), cy = 0.5 * (loy + hiy);
  const double r = std::hypot(hix - lox, hiy - loy) + 1.0;

  std::vector<Pt> pts(P.begin(), P.begin() + static_cast<long>(nReal));
  const std::size_t s0 = pts.size();
  pts.push_back({cx, cy + 3 * r});
  pts.push_back({cx - 3 * r, cy - 2 * r});
  pts.push_back({cx + 3 * r, cy - 2 * r});

  std::vector<Tri> tris{{s0, s0 + 1, s0 + 2}};
  auto key = [](std::size_t a, std::size_t b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
  };
  for (std::size_t i = 0; i < nReal; ++i) {
    std::vector<std::size_t> bad;
    for (std::size_t t = 0; t < tris.size(); ++t)
      if (inCircumcircle(pts[tris[t][0]], pts[tris[t][1]], pts[tris[t][2]],
                         pts[i]))
        bad.push_back(t);
    if (bad.empty())
      continue;
    std::map<std::pair<std::size_t, std::size_t>, int> count;
    for (std::size_t t : bad) {
      ++count[key(tris[t][0], tris[t][1])];
      ++count[key(tris[t][1], tris[t][2])];
      ++count[key(tris[t][2], tris[t][0])];
    }
    std::set<std::size_t> badset(bad.begin(), bad.end());
    std::vector<Tri> kept;
    kept.reserve(tris.size());
    for (std::size_t t = 0; t < tris.size(); ++t)
      if (!badset.count(t))
        kept.push_back(tris[t]);
    for (const auto &kv : count) {
      if (kv.second != 1)
        continue;
      const std::size_t u = kv.first.first, v = kv.first.second;
      if (cross2(pts[u], pts[v], pts[i]) > 0.0)
        kept.push_back({u, v, i});
      else if (cross2(pts[v], pts[u], pts[i]) > 0.0)
        kept.push_back({v, u, i});
    }
    tris.swap(kept);
  }
  std::vector<Tri> out;
  for (const Tri &t : tris)
    if (t[0] < s0 && t[1] < s0 && t[2] < s0)
      out.push_back(t);
  return out;
}

double bestSpacingHint(const std::vector<Pt> &outline) {
  const std::size_t n = outline.size();
  std::vector<double> seg(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Pt &a = outline[i], &b = outline[(i + 1) % n];
    seg[i] = std::hypot(b[0] - a[0], b[1] - a[1]);
  }
  std::nth_element(seg.begin(), seg.begin() + static_cast<long>(n / 2),
                   seg.end());
  return 1.6 * seg[n / 2];
}

} // namespace

std::vector<double> triangleQuality(const PolygonMesh &m) {
  std::vector<double> q;
  q.reserve(m.tris.size());
  for (const Tri &t : m.tris)
    q.push_back(radiusRatio(m.points[t[0]], m.points[t[1]], m.points[t[2]]));
  return q;
}

namespace {
PolygonMesh triangulateAt(const std::vector<Pt> &outline, double spacing);
} // namespace

// Sweep the interior spacing and keep the best mesh.
//
// One spacing is not enough, for two reasons found by testing. Too coarse and
// no interior point fits clear of the walls at all, so the result falls back
// to a bare ear clipping whose slivers reach a radius ratio of 99. Too fine
// and points crowd the boundary and make slivers of their own. Quality is not
// monotone in the spacing either -- on a narrow-toothed test polygon it runs
// 22, 10.9, 7.7, 17.3 as the spacing falls -- so this searches rather than
// shrinks until satisfied.
//
// Refining near a reflex corner is what a Ruppert-style mesher would do by
// splitting the boundary segment, which is not available here: those vertices
// are shared walls, and adding one would change the tissue's topology. The
// interior spacing is the only lever, and where the outline is too coarse
// around a narrow neck for any spacing to work, the bound cannot be met --
// that is a property of the outline, not of this code, and the caller is told
// by the quality it gets back.
PolygonMesh triangulatePolygon(const std::vector<Pt> &outline,
                               double qualityBound, double spacing) {
  const double h0 = spacing > 0.0 ? spacing : bestSpacingHint(outline);
  std::vector<PolygonMesh> tried;
  std::vector<double> worsts;
  double h = h0;
  for (int attempt = 0; attempt < 8; ++attempt, h *= 0.72) {
    PolygonMesh m = triangulateAt(outline, h);
    if (m.tris.empty())
      continue;
    double worst = 0.0;
    for (double q : triangleQuality(m))
      worst = std::max(worst, q);
    tried.push_back(std::move(m));
    worsts.push_back(worst);
    if (worst <= qualityBound)
      break; // the coarsest spacing that meets the bound; no reason to refine
  }
  if (tried.empty())
    return PolygonMesh{};

  // Take the cheapest mesh of acceptable quality, not the best mesh at any
  // price. Where the bound is reachable the loop above has already stopped at
  // the coarsest spacing that reaches it. Where it is not -- a cell whose
  // outline is too coarse around its own narrowest neck -- the difference
  // between the best attempt and one slightly worse is a handful of
  // hundredths in radius ratio and can be several times the vertex count,
  // and every one of those vertices is a degree of freedom the relaxation
  // then has to carry. Keeping the best regardless took a 40-cell patch to
  // 4818 interior vertices for a worst element of 11.66, where a 15% wider
  // tolerance buys nearly all of it for a fraction of the cost.
  // A ceiling on interior vertices, because every one of them is a degree of
  // freedom the relaxation carries for the rest of the run, and because rows
  // are padded to the widest cell, so a single greedy cell sets the width for
  // the whole tissue. Where the bound is unreachable -- an outline too coarse
  // around its own neck -- the sweep will otherwise refine to the end of its
  // range chasing it: one cell of a 40-cell patch reached 2338 interior
  // vertices, padding 40 rows to 93520 slots to hold 4983 real ones.
  //
  // Stopping the sweep early instead was tried and is wrong: quality is not
  // monotone in the spacing, so a step that fails to improve is not evidence
  // that later ones will not, and breaking on one cost the 10-point star its
  // bound entirely.
  const std::size_t n = outline.size();
  const std::size_t ceiling = 4 * n;

  double bestWorst = worsts[0];
  for (double w : worsts)
    bestWorst = std::min(bestWorst, w);
  const double allow = std::max(qualityBound, 1.15 * bestWorst);
  // Order of preference: inside both limits, then inside the vertex ceiling
  // at whatever quality that allows, and only then over the ceiling. The
  // ceiling outranks the quality allowance deliberately. A degree of freedom
  // is paid for on every force evaluation for the rest of the run, and paid
  // for by every other cell too through the row padding, whereas a worse
  // element costs conditioning once. Ordering it the other way round left a
  // single cell holding 2338 interior vertices because no mesh satisfied both.
  std::size_t pick = 0;
  bool have = false;
  for (int pass = 0; pass < 3 && !have; ++pass) {
    for (std::size_t i = 0; i < tried.size(); ++i) {
      const std::size_t inner = tried[i].points.size() - tried[i].numBoundary;
      if (pass < 2 && inner > ceiling)
        continue;
      if (pass != 1 && worsts[i] > allow)
        continue;
      const bool better =
          !have || (pass == 1 ? worsts[i] < worsts[pick]
                              : tried[i].points.size() < tried[pick].points.size());
      if (better) {
        pick = i;
        have = true;
      }
    }
  }
  if (!have)
    for (std::size_t i = 0; i < tried.size(); ++i)
      if (worsts[i] == bestWorst) { pick = i; break; }
  return std::move(tried[pick]);
}

namespace {
PolygonMesh triangulateAt(const std::vector<Pt> &outline, double spacing) {
  PolygonMesh out;
  if (outline.size() < 3)
    return out;

  std::vector<Pt> P = outline;
  if (signedArea(P) < 0.0)
    std::reverse(P.begin(), P.end());
  const std::size_t n = P.size();

  // Edge lengths set the interior spacing, so the mesh is as fine as the
  // outline it has to match and no finer.
  std::vector<double> seg(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Pt &a = P[i], &b = P[(i + 1) % n];
    seg[i] = std::hypot(b[0] - a[0], b[1] - a[1]);
  }
  std::vector<double> sorted = seg;
  std::nth_element(sorted.begin(), sorted.begin() + static_cast<long>(n / 2),
                   sorted.end());
  const double median = sorted[n / 2];
  const double h = spacing > 0.0 ? spacing : 1.6 * median;

  out.points = P;
  out.numBoundary = n;

  // Interior points on a lattice, kept clear of the outline: a Steiner point
  // close to the boundary makes exactly the sliver this is here to avoid.
  if (h > 0.0) {
    double lox = P[0][0], hix = P[0][0], loy = P[0][1], hiy = P[0][1];
    for (const Pt &p : P) {
      lox = std::min(lox, p[0]); hix = std::max(hix, p[0]);
      loy = std::min(loy, p[1]); hiy = std::max(hiy, p[1]);
    }
    const double dy = h * std::sqrt(3.0) / 2.0; // triangular lattice
    std::size_t row = 0;
    for (double y = loy + dy * 0.5; y < hiy; y += dy, ++row) {
      const double off = (row % 2) ? 0.5 * h : 0.0;
      for (double x = lox + h * 0.5 + off; x < hix; x += h) {
        const Pt p{x, y};
        if (!pointInPolygon(p, P))
          continue;
        double dmin = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; ++i)
          dmin = std::min(dmin, segDistance(p, P[i], P[(i + 1) % n]));
        if (dmin > 0.55 * h)
          out.points.push_back(p);
      }
    }
  }

  // Triangulate the whole point set at once and keep what lies inside the
  // polygon. Ear clipping is the fallback, not the starting point: it is
  // valid for any simple polygon but its choices cannot be undone by later
  // insertion, and on a segmented outline it bridges nearly collinear runs
  // into slivers nothing can remove.
  out.tris = bowyerWatson(out.points, out.points.size());
  {
    std::vector<Tri> inside;
    for (const Tri &t : out.tris) {
      const Pt c{(out.points[t[0]][0] + out.points[t[1]][0] + out.points[t[2]][0]) / 3.0,
                 (out.points[t[0]][1] + out.points[t[1]][1] + out.points[t[2]][1]) / 3.0};
      if (pointInPolygon(c, P) && cross2(out.points[t[0]], out.points[t[1]],
                                         out.points[t[2]]) > 0.0)
        inside.push_back(t);
    }
    out.tris.swap(inside);
  }

  std::set<std::pair<std::size_t, std::size_t>> fixed;
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t a = i, b = (i + 1) % n;
    fixed.insert({std::min(a, b), std::max(a, b)});
  }

  // Every outline edge has to be an edge of the mesh, or the face does not
  // meet its own walls. The interior clearance makes the outline edges
  // locally Delaunay so they survive; if one does not, take the ear clipping,
  // which conforms by construction.
  {
    std::set<std::pair<std::size_t, std::size_t>> have;
    for (const Tri &t : out.tris) {
      have.insert({std::min(t[0], t[1]), std::max(t[0], t[1])});
      have.insert({std::min(t[1], t[2]), std::max(t[1], t[2])});
      have.insert({std::min(t[2], t[0]), std::max(t[2], t[0])});
    }
    bool conforms = true;
    for (const auto &e : fixed)
      if (!have.count(e)) {
        conforms = false;
        break;
      }
    if (!conforms) {
      out.points.resize(n);
      out.tris = earClip(P);
    }
  }

  // Insert each interior point by a constrained Bowyer-Watson cavity.
  //
  // The cavity has to be grown by walking outward from the triangle that
  // contains the point, never stepping across an outline edge. Taking every
  // triangle whose circumcircle contains the point -- the unconstrained test
  // -- is wrong for a non-convex polygon: a circumcircle happily reaches
  // across a reflex corner and picks up triangles on the far side of the
  // cell, and the cavity is then not star-shaped about the point, so
  // retriangulating it leaves holes and overlaps. That showed up as a jigsaw
  // outline whose triangles did not add up to its own area.
  auto key = [](std::size_t a, std::size_t b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
  };
  // Adjacency is built once and maintained. Rebuilding it per inserted point
  // made this quadratic in the interior count and cost twenty seconds a cell,
  // against the milliseconds the whole scheme is predicated on.
  std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> adj;
  auto adjAdd = [&](std::size_t t) {
    const Tri &T = out.tris[t];
    adj[key(T[0], T[1])].push_back(t);
    adj[key(T[1], T[2])].push_back(t);
    adj[key(T[2], T[0])].push_back(t);
  };
  auto adjDrop = [&](std::size_t t) {
    const Tri &T = out.tris[t];
    const std::pair<std::size_t, std::size_t> es[3] = {
        key(T[0], T[1]), key(T[1], T[2]), key(T[2], T[0])};
    for (const auto &e : es) {
      auto &v = adj[e];
      v.erase(std::remove(v.begin(), v.end(), t), v.end());
    }
  };
  for (std::size_t t = 0; t < out.tris.size(); ++t)
    adjAdd(t);

  for (std::size_t k = n; k < out.points.size(); ++k) {
    const Pt &p = out.points[k];

    std::size_t seed = SIZE_MAX;
    for (std::size_t t = 0; t < out.tris.size(); ++t) {
      const Tri &T = out.tris[t];
      if (pointInTriangle(p, out.points[T[0]], out.points[T[1]],
                          out.points[T[2]])) {
        seed = t;
        break;
      }
    }
    if (seed == SIZE_MAX)
      continue; // outside the current mesh; skip rather than corrupt it

    std::set<std::size_t> cavity{seed};
    std::vector<std::size_t> stack{seed};
    while (!stack.empty()) {
      const std::size_t t = stack.back();
      stack.pop_back();
      const Tri &T = out.tris[t];
      const std::pair<std::size_t, std::size_t> es[3] = {
          key(T[0], T[1]), key(T[1], T[2]), key(T[2], T[0])};
      for (const auto &e : es) {
        if (fixed.count(e))
          continue; // an outline edge bounds the cavity, always
        for (std::size_t nb : adj[e]) {
          if (cavity.count(nb))
            continue;
          const Tri &N = out.tris[nb];
          if (!inCircumcircle(out.points[N[0]], out.points[N[1]],
                              out.points[N[2]], p))
            continue;
          cavity.insert(nb);
          stack.push_back(nb);
        }
      }
    }

    std::map<std::pair<std::size_t, std::size_t>, int> count;
    for (std::size_t t : cavity) {
      const Tri &T = out.tris[t];
      ++count[key(T[0], T[1])];
      ++count[key(T[1], T[2])];
      ++count[key(T[2], T[0])];
    }
    std::vector<Tri> added;
    bool ok = true;
    for (const auto &kv : count) {
      if (kv.second != 1)
        continue; // interior to the cavity, disappears
      const std::size_t u = kv.first.first, v = kv.first.second;
      const double or1 = cross2(out.points[u], out.points[v], p);
      if (or1 > 0.0)
        added.push_back({u, v, k});
      else if (cross2(out.points[v], out.points[u], p) > 0.0)
        added.push_back({v, u, k});
      else {
        ok = false; // point on the cavity boundary; leave it out
        break;
      }
    }
    if (!ok || added.empty())
      continue;

    // Overwrite the cavity's slots in place and append the rest, so triangle
    // indices outside the cavity keep their meaning and the adjacency stays
    // valid without being rebuilt.
    std::vector<std::size_t> slots(cavity.begin(), cavity.end());
    for (std::size_t t : slots)
      adjDrop(t);
    std::size_t w = 0;
    for (; w < added.size() && w < slots.size(); ++w) {
      out.tris[slots[w]] = added[w];
      adjAdd(slots[w]);
    }
    for (std::size_t i = w; i < added.size(); ++i) {
      out.tris.push_back(added[i]);
      adjAdd(out.tris.size() - 1);
    }
    // Any cavity slot left over becomes a degenerate placeholder, removed
    // below; leaving a stale triangle in place would corrupt the mesh.
    for (std::size_t i = w; i < slots.size(); ++i)
      out.tris[slots[i]] = {SIZE_MAX, SIZE_MAX, SIZE_MAX};
  }
  out.tris.erase(std::remove_if(out.tris.begin(), out.tris.end(),
                                [](const Tri &t) { return t[0] == SIZE_MAX; }),
                 out.tris.end());

  delaunayFlips(out.points, out.tris, fixed);

  // Drop degenerate triangles. They should not be produced at all, and the
  // assertion in tests/meshing_test.cpp now says so, but a mesh that leaves
  // with one is worse than a mesh that is one triangle smaller: the element
  // it becomes has a rest state no triangle can adopt, which is the failure
  // this whole file exists to prevent.
  {
    std::vector<Tri> ok;
    ok.reserve(out.tris.size());
    for (const Tri &t : out.tris) {
      if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2])
        continue;
      // Only genuine degeneracy: a zero-length edge or an area that is zero
      // against the element's own size. Filtering on the quality measure
      // instead also throws away thin but valid slivers, and the holes that
      // leaves are worse than the slivers were -- it cost this file's own
      // area check on a 7-lobed test polygon.
      const Pt &p0 = out.points[t[0]], &p1 = out.points[t[1]],
               &p2 = out.points[t[2]];
      const double e0 = std::hypot(p1[0]-p0[0], p1[1]-p0[1]);
      const double e1 = std::hypot(p2[0]-p1[0], p2[1]-p1[1]);
      const double e2 = std::hypot(p0[0]-p2[0], p0[1]-p2[1]);
      const double sc = std::max(e0, std::max(e1, e2));
      if (e0 <= 0.0 || e1 <= 0.0 || e2 <= 0.0)
        continue;
      if (std::fabs(cross2(p0, p1, p2)) <= 1e-12 * sc * sc)
        continue;
      ok.push_back(t);
    }
    out.tris.swap(ok);
  }

  // Drop any point the insertion refused, so the mesh has no orphans.
  std::vector<bool> used(out.points.size(), false);
  for (const Tri &t : out.tris)
    for (std::size_t v : t)
      used[v] = true;
  if (std::find(used.begin() + static_cast<long>(n), used.end(), false) !=
      used.end()) {
    std::vector<std::size_t> remap(out.points.size(), SIZE_MAX);
    std::vector<Pt> pts;
    for (std::size_t i = 0; i < out.points.size(); ++i)
      if (i < n || used[i]) {
        remap[i] = pts.size();
        pts.push_back(out.points[i]);
      }
    for (Tri &t : out.tris)
      for (std::size_t &v : t)
        v = remap[v];
    out.points.swap(pts);
  }
  return out;
}
} // namespace

} // namespace tissue
