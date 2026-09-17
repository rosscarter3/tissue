#ifndef TISSUE2_CORE_TYPES_H
#define TISSUE2_CORE_TYPES_H

#include <array>
#include <cmath>
#include <cstddef>
#include <span>

namespace tissue {

// Fixed-capacity position/direction; only the first `dimension` (2 or 3)
// components are meaningful. Avoids per-call heap allocation in hot loops.
using Vec3 = std::array<double, 3>;

inline double squaredDistance(std::span<const double> a,
                              std::span<const double> b) {
  double s = 0.0;
  for (size_t d = 0; d < a.size(); ++d) {
    double diff = a[d] - b[d];
    s += diff * diff;
  }
  return s;
}

inline double distance(std::span<const double> a, std::span<const double> b) {
  return std::sqrt(squaredDistance(a, b));
}

} // namespace tissue

#endif
