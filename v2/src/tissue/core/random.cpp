#include "tissue/core/random.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace tissue::random {

// Numerical Recipes ran3 (Knuth subtractive generator), ported verbatim from
// the legacy myRandom.cc so seeded runs reproduce identical sequences —
// including its quirks: default seed 1 with lazy initialization, and 0 being
// mapped to FAC.
namespace {
constexpr long MBIG = 1000000000;
constexpr long MSEED = 161803398;
constexpr long MZ = 0;
constexpr double FAC = 1.0 / MBIG;
constexpr double kPi = 3.14159265; // legacy myMath::pi() literal

long idum = 1;
int inext, inextp;
long ma[56];
int iff = 0;
} // namespace

double ran3() {
  long mj, mk;
  if (idum < 0 || iff == 0) {
    iff = 1;
    mj = MSEED - (idum < 0 ? -idum : idum);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (int i = 1; i <= 54; ++i) {
      int ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (int k = 1; k <= 4; ++k)
      for (int i = 1; i <= 55; ++i) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    inext = 0;
    inextp = 31;
    idum = 1;
  }
  if (++inext == 56)
    inext = 1;
  if (++inextp == 56)
    inextp = 1;
  mj = ma[inext] - ma[inextp];
  if (mj < MZ)
    mj += MBIG;
  ma[inext] = mj;
  double ret = mj * FAC;
  if (mj == 0)
    ret = FAC;
  return ret;
}

void sran3(long seed) {
  std::cerr << "sran3()\n";
  if (seed < 1 || seed > 10000000) {
    std::cerr << "tissue::random::sran3() Seed (idum) provided in dangerous "
              << "region. Use values in [1:10000000]." << std::endl;
    std::exit(-1);
  }
  idum = -seed;
}

double Rnd() { return ran3(); }

double Grand() {
  // Explicit sequencing (legacy relies on unspecified evaluation order; this
  // matches left-to-right, which is what clang produced for the old binary).
  double u1 = Rnd();
  double u2 = Rnd();
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(kPi * u2);
}

} // namespace tissue::random
