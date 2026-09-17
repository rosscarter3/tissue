//
// Legacy-compatible global random number generator (Numerical Recipes ran3,
// the subtractive Knuth generator used by the original myRandom.cc). Kept
// bit-compatible so reruns of old models with seeded reactions reproduce the
// same sequences.
//
#ifndef TISSUE2_CORE_RANDOM_H
#define TISSUE2_CORE_RANDOM_H

namespace tissue::random {

// Seeds the global ran3 state; the seed must lie in [1, 10000000] (exits
// otherwise), matching legacy myRandom::sran3 semantics.
void sran3(long seed);

// Uniform deviate in (0, 1) (a raw 0 is mapped to 1/MBIG, as in legacy).
double ran3();

// Legacy aliases used across reactions and division rules.
double Rnd();   // = ran3()
double Grand(); // legacy quasi Box-Muller: sqrt(-2 ln u1) * cos(pi u2)

} // namespace tissue::random

#endif
