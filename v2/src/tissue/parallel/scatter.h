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
#ifndef TISSUE2_PARALLEL_SCATTER_H
#define TISSUE2_PARALLEL_SCATTER_H

#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/parallel/thread_pool.h"

namespace tissue {

// fn(begin, end, targets...) accumulates the contribution of items
// [begin, end) into the given target tables (+= only).
template <typename Fn>
void parallelScatter2(size_t n, Matrix &a, Matrix &b, Fn fn,
                      size_t grain = ThreadPool::defaultGrain) {
  ThreadPool &pool = ThreadPool::instance();
  size_t threads = pool.numThreads();
  if (threads == 1 || n < grain) {
    fn(size_t{0}, n, a, b);
    return;
  }
  size_t parts = std::min(threads, n);
  std::vector<Matrix> aParts(parts), bParts(parts);
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
  if (threads == 1 || n < grain) {
    fn(size_t{0}, n, a);
    return;
  }
  size_t parts = std::min(threads, n);
  std::vector<Matrix> aParts(parts);
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
