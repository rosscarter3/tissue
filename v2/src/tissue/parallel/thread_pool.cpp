#include "tissue/parallel/thread_pool.h"

#include <algorithm>
#include <string>

namespace tissue {

ThreadPool &ThreadPool::instance() {
  static ThreadPool pool([] {
    if (const char *env = std::getenv("TISSUE_NUM_THREADS")) {
      long v = std::strtol(env, nullptr, 10);
      if (v > 0)
        return static_cast<size_t>(v - 1); // v total threads incl. caller
    }
    unsigned hw = std::thread::hardware_concurrency();
    return static_cast<size_t>(hw > 1 ? hw - 1 : 0);
  }());
  return pool;
}

ThreadPool::ThreadPool(size_t numWorkers) {
  workers_.reserve(numWorkers);
  for (size_t i = 0; i < numWorkers; ++i)
    workers_.emplace_back([this, i] { workerLoop(i); });
}

ThreadPool::~ThreadPool() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stop_ = true;
  }
  wake_.notify_all();
  for (auto &w : workers_)
    w.join();
}

void ThreadPool::workerLoop(size_t workerIndex) {
  uint64_t seenEpoch = 0;
  for (;;) {
    Job job;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      wake_.wait(lock, [&] { return stop_ || job_.epoch > seenEpoch; });
      if (stop_)
        return;
      job = job_;
      seenEpoch = job.epoch;
    }
    // Partition p = workerIndex + 1 (caller runs partition 0).
    size_t part = workerIndex + 1;
    if (part < job.parts) {
      size_t chunk = (job.n + job.parts - 1) / job.parts;
      size_t begin = std::min(part * chunk, job.n);
      size_t end = std::min(begin + chunk, job.n);
      if (begin < end)
        (*job.fn)(begin, end, part);
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (--pending_ == 0)
        done_.notify_one();
    }
  }
}

void ThreadPool::parallelFor(
    size_t n, size_t grain,
    const std::function<void(size_t, size_t, size_t)> &fn) {
  size_t threads = numThreads();
  if (n == 0)
    return;
  if (threads == 1 || n < grain) {
    fn(0, n, 0);
    return;
  }
  size_t parts = std::min(threads, n);
  {
    std::lock_guard<std::mutex> lock(mutex_);
    job_.fn = &fn;
    job_.n = n;
    job_.parts = parts;
    ++job_.epoch;
    pending_ = workers_.size();
  }
  wake_.notify_all();
  // Caller executes partition 0.
  size_t chunk = (n + parts - 1) / parts;
  fn(0, std::min(chunk, n), 0);
  std::unique_lock<std::mutex> lock(mutex_);
  done_.wait(lock, [&] { return pending_ == 0; });
}

} // namespace tissue
