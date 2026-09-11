//
// Persistent thread pool with deterministic static partitioning.
//
// Design goals:
//  - Zero external dependencies (std::thread only).
//  - Deterministic results for a fixed thread count: work is split into one
//    contiguous block per worker, and reductions sum per-worker buffers in
//    worker order, so floating point addition order never depends on timing.
//  - Small problems run serially (grain threshold) so tiny models pay no
//    synchronization cost.
//
#ifndef TISSUE2_PARALLEL_THREAD_POOL_H
#define TISSUE2_PARALLEL_THREAD_POOL_H

#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace tissue {

class ThreadPool {
public:
  // Global pool. Thread count: TISSUE_NUM_THREADS env var if set (>0),
  // otherwise hardware concurrency. TISSUE_NUM_THREADS=1 disables threading.
  static ThreadPool &instance();

  size_t numThreads() const { return workers_.size() + 1; } // workers + caller

  // Runs fn(begin, end, part) over [0, n) split into numThreads() contiguous
  // blocks; part is the block index in [0, numThreads()). Falls back to a
  // single serial call fn(0, n, 0) when n < grain or only one thread exists.
  // Blocks until all parts are done. Not reentrant.
  void parallelFor(size_t n, size_t grain,
                   const std::function<void(size_t, size_t, size_t)> &fn);

  // Default work-item threshold below which parallelFor runs serially.
  static constexpr size_t defaultGrain = 1024;

  ~ThreadPool();

private:
  explicit ThreadPool(size_t numWorkers);
  void workerLoop(size_t workerIndex);

  struct Job {
    const std::function<void(size_t, size_t, size_t)> *fn = nullptr;
    size_t n = 0;
    size_t parts = 0;
    uint64_t epoch = 0;
  };

  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable wake_;
  std::condition_variable done_;
  Job job_; // job_.epoch increments per parallelFor call
  size_t pending_ = 0; // workers still running current job
  bool stop_ = false;
};

// Convenience: parallel loop over [0, n) ignoring the partition index.
inline void parallelFor(size_t n, const std::function<void(size_t, size_t)> &fn,
                        size_t grain = ThreadPool::defaultGrain) {
  ThreadPool::instance().parallelFor(
      n, grain, [&fn](size_t b, size_t e, size_t) { fn(b, e); });
}

} // namespace tissue

#endif
