/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_DEFINED
#define COULOMB_INTEGRALS_DEFINED

#include "PerturbationTensor.hpp"
#include "Chi.hpp"
#include <ctf.hpp>

class CoulombIntegrals: public PerturbationTensor {
  public:
    CoulombIntegrals(Chi *chiReal, Chi *chiImag);
    virtual ~CoulombIntegrals();

    virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);

    CTF::Tensor<> getSlice(int a, int b);

    void fetch();

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts a global index of a tensor entry into the positions
     * along each dimension.
     */
    // FIXME: cannot handle out of bound indicies
    void from_index(CTF::Tensor<> const &t, int64_t index, int *pos) {
      for (int d(0); d < t.order; ++d) {
        pos[d] = index % (int64_t)t.lens[d];
        index /= t.lens[d];
      }
    }

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts given positions along each dimension into a global index.
     */
    int64_t to_index(int dim, int const *lens, int const *pos) {
      int64_t index(0);
      for (int d(dim-1); d >= 0; --d) {
        index *= lens[d];
        index += pos[d];
      }
      return index;
    }

    Chi *chiReal, *chiImag;
    CTF::Tensor<> *a, *i, *ai, *abij;
// NOTE: only for testing
    CTF::Tensor<> *abcd;
};

#endif

