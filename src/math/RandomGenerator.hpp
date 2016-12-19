#ifndef RANDOM_GENERATOR_DEFINED
#define RANDOM_GENERATOR_DEFINED

#include <cstdint>

/**
 * Provides a random number generator based on the mzran algorithm
 * marsaglia, g. and a. zaman. 1994.
 * some portable very-long period random number generators.
 * computers in physics.  8(1): 117-121.
 **/
namespace cc4s {
  class RandomGenerator {
  public:
    RandomGenerator(
      unsigned int seed = 1, int shaking = 32
    ):
      x(seed ^  521288629),
      y(seed ^  362436069),
      z(seed ^   16163801),
      n(seed ^ 1131199209),
      c(y > x)
    {
      for (int i(0); i < shaking; ++i) nextInteger();
    }

    double nextUniform() {
      uint64_t i(nextCoarseInteger() >> 20);
      for (int j(1); j < 6; ++j) {
        i <<= 8;
        i |= nextCoarseInteger() >> 24;
      }
      return static_cast<double>(i) / static_cast<double>(1l << 52);
    }

    /**
     * Returns the next integer with 32 random bits where all bits show good
     * randomness.
     **/
    unsigned int nextInteger() {
      unsigned int i(nextCoarseInteger() >> 24);
      // draw 4 random numbers using 8 of their most significant bits
      for (int j(1); j < 4; ++j) {
        i <<= 8;
        i |= nextCoarseInteger() >> 24;
      }
      return i;
    }

    /**
     * Returns the next integer with 32 random bits where only the most
     * significant bits show good randomness.
     **/
    unsigned int nextCoarseInteger() {
      unsigned int s;
      if (y - c > x) {
        s = y - c - x;
        c = 0;
      } else {
        s = y - c - x - 18;
        c = 1;
      }
      x = y;
      y = z;
      z = s;
      n = 69069 * n + 1013904243;
      return s + n;
    }

  protected:
    unsigned int x, y, z, n, c;
  };
}

#endif

