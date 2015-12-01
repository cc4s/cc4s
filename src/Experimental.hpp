/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef EXPERIMENTAL_DEFINED
#define EXPERIMENTAL_DEFINED

namespace cc4s {
  class Experimental {
    public:
      // TODO: check if the cft framework offers such a function natively
      /**
       * \brief Converts a global index of a tensor entry into the positions
       * along each dimension.
       */
      // FIXME: cannot handle out of bound indicies
      static void from_index(CTF::Tensor<> const &t, int64_t index, int *pos) {
        for (int d(0); d < t.order; ++d) {
          pos[d] = index % int64_t(t.lens[d]);
          index /= t.lens[d];
        }
      }

      // TODO: check if the cft framework offers such a function natively
      /**
       * \brief Converts given positions along each dimension into a global index.
       */
      static int64_t to_index(int dim, int const *lens, int const *pos) {
        int64_t index(0);
        for (int d(dim-1); d >= 0; --d) {
          index *= lens[d];
          index += pos[d];
        }
        return index;
      }
      static void readRandom(CTF::Tensor<> *tensor, int seed);
      static void testSymmetries();
  };
}

#endif
