//
// Deterministic parallel scatter-accumulation.
//
// Wall/edge loops accumulate forces into shared vertex/cell derivative rows,
// so a naive parallel-for races. parallelScatter gives each partition a
// private zeroed copy of the target tables and reduces the partials in
// partition order afterwards — deterministic for a fixed thread count.
// Below the grain threshold it runs serially on the shared targets (zero
// overhead, bit-identical to a serial loop).
//
// Two costs decide that threshold, and only one of them scales with the loop:
//   - the parallel region itself (mutex + broadcast + join), and
//   - the private copies: zeroing parts*|target| doubles and reducing them
//     back in partition order.
// The second is independent of n, so a wall loop whose target is the vertex
// table only pays off once n is comfortably larger than parts*|target|. On the
// hook shell (1040 walls scattering into 528x3 vertex derivatives) the
// reduction alone moves more doubles than the walls do, which is why the
// serial path is taken there. The private copies are kept in a reusable
// scratch cache: reshapeLike() re-zeroes without reallocating when the shape
// is unchanged, so the steady state performs no allocation at all.
//
#ifndef TISSUE2_PARALLEL_SCATTER_H
#define TISSUE2_PARALLEL_SCATTER_H

#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/parallel/thread_pool.h"

namespace tissue {

namespace detail {

// Scratch for the private per-partition copies. parallelScatter is always
// called from the driving thread (the pool workers run inside it), so a
// function-local cache per target arity is safe and keeps the steady state
// allocation-free.
inline std::vector<Matrix> &scatterScratch(size_t which, size_t parts) {
  static thread_local std::vector<std::vector<Matrix>> cache(3);
  std::vector<Matrix> &v = cache[which];
  if (v.size() < parts)
    v.resize(parts);
  return v;
}

// Parallelizing a scatter only pays when the element work dominates the fixed
// cost of zeroing and reducing the private copies. Balancing
//   n*w*(1 - 1/t)  against  t*S*r + F
// with per-element work w ~ 30 ns, reduce cost r ~ 3 ns and region overhead
// F ~ 20 us puts the crossover near n ~ S, so require the loop to be at least
// as long as the target table is wide, on top of the plain grain.
//
// Note this rarely fires: the private-copy strategy costs O(threads * |target|)
// no matter how the loop is split, which is the real reason wall-force
// threading shows little gain even at 40k cells. A two-pass gather
// (per-wall coefficients, then per-vertex summation over the incidence lists)
// has no such term and is the scalable alternative.
inline bool scatterWorthIt(size_t n, size_t grain, size_t targetElems) {
  return n >= grain && n >= targetElems;
}

} // namespace detail

// fn(begin, end, targets...) accumulates the contribution of items
// [begin, end) into the given target tables (+= only).
template <typename Fn>
void parallelScatter2(size_t n, Matrix &a, Matrix &b, Fn fn,
                      size_t grain = ThreadPool::defaultGrain) {
  ThreadPool &pool = ThreadPool::instance();
  size_t threads = pool.numThreads();
  size_t parts = std::min(threads, n);
  if (threads == 1 ||
      !detail::scatterWorthIt(n, grain,
                              a.flat().size() + b.flat().size())) {
    fn(size_t{0}, n, a, b);
    return;
  }
  std::vector<Matrix> &aParts = detail::scatterScratch(0, parts);
  std::vector<Matrix> &bParts = detail::scatterScratch(1, parts);
  for (size_t p = 0; p < parts; ++p) {
    aParts[p].reshapeLike(a);
    bParts[p].reshapeLike(b);
  }
  pool.parallelFor(n, 1, [&](size_t begin, size_t end, size_t p) {
    fn(begin, end, aParts[p], bParts[p]);
  });
  // Deterministic reduction in partition order.
  auto aFlat = a.flat();
  auto bFlat = b.flat();
  for (size_t p = 0; p < parts; ++p) {
    auto ap = aParts[p].flat();
    auto bp = bParts[p].flat();
    for (size_t k = 0; k < aFlat.size(); ++k)
      aFlat[k] += ap[k];
    for (size_t k = 0; k < bFlat.size(); ++k)
      bFlat[k] += bp[k];
  }
}

template <typename Fn>
void parallelScatter1(size_t n, Matrix &a, Fn fn,
                      size_t grain = ThreadPool::defaultGrain) {
  ThreadPool &pool = ThreadPool::instance();
  size_t threads = pool.numThreads();
  size_t parts = std::min(threads, n);
  if (threads == 1 || !detail::scatterWorthIt(n, grain, a.flat().size())) {
    fn(size_t{0}, n, a);
    return;
  }
  std::vector<Matrix> &aParts = detail::scatterScratch(2, parts);
  for (size_t p = 0; p < parts; ++p)
    aParts[p].reshapeLike(a);
  pool.parallelFor(n, 1, [&](size_t begin, size_t end, size_t p) {
    fn(begin, end, aParts[p]);
  });
  auto aFlat = a.flat();
  for (size_t p = 0; p < parts; ++p) {
    auto ap = aParts[p].flat();
    for (size_t k = 0; k < aFlat.size(); ++k)
      aFlat[k] += ap[k];
  }
}

} // namespace tissue

#endif
