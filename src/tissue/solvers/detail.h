//
// Flat element-wise kernels shared by the solvers. All state tables in one
// group (cell/wall/vertex) share shapes with their derivative/scratch
// counterparts, so the kernels run over the contiguous storage — the hot
// loops the compiler auto-vectorizes, parallelized above a grain threshold.
//
#ifndef TISSUE2_SOLVERS_DETAIL_H
#define TISSUE2_SOLVERS_DETAIL_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "tissue/core/matrix.h"
#include "tissue/parallel/thread_pool.h"

namespace tissue::detail {

// Elements below which loops run serial; shares the pool's measured
// crossover so one knob (TISSUE_GRAIN) tunes every parallel region.
inline const size_t kGrain = ThreadPool::defaultGrain;

// y := x + a*d
inline void assignAxpy(Matrix &y, const Matrix &x, const Matrix &d, double a) {
  auto yf = y.flat();
  auto xf = x.flat();
  auto df = d.flat();
  parallelFor(
      yf.size(),
      [&](size_t b, size_t e) {
        for (size_t k = b; k < e; ++k)
          yf[k] = xf[k] + a * df[k];
      },
      kGrain);
}

// y += a*d
inline void addScaled(Matrix &y, const Matrix &d, double a) {
  auto yf = y.flat();
  auto df = d.flat();
  parallelFor(
      yf.size(),
      [&](size_t b, size_t e) {
        for (size_t k = b; k < e; ++k)
          yf[k] += a * df[k];
      },
      kGrain);
}

// max over elements of |err/scale|; max is order-independent, so the parallel
// reduction is deterministic for any thread count.
inline double maxErrRatio(const Matrix &err, const Matrix &scale) {
  auto ef = err.flat();
  auto sf = scale.flat();
  ThreadPool &pool = ThreadPool::instance();
  size_t parts = pool.numThreads();
  if (parts == 1 || ef.size() < kGrain) {
    double errMax = 0.0;
    for (size_t k = 0; k < ef.size(); ++k) {
      double aux = std::fabs(ef[k] / sf[k]);
      if (aux > errMax)
        errMax = aux;
    }
    return errMax;
  }
  std::vector<double> partMax(parts, 0.0);
  pool.parallelFor(ef.size(), 1, [&](size_t b, size_t e, size_t p) {
    double m = 0.0;
    for (size_t k = b; k < e; ++k) {
      double aux = std::fabs(ef[k] / sf[k]);
      if (aux > m)
        m = aux;
    }
    partMax[p] = m;
  });
  double errMax = 0.0;
  for (double m : partMax)
    errMax = std::max(errMax, m);
  return errMax;
}

} // namespace tissue::detail

#endif
