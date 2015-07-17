/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_DEFINED
#define COULOMB_INTEGRALS_DEFINED

#include "Chi.hpp"
#include <ctf.hpp>

/**
 * \brief Parts of the full tensors V_p, V_pq, V_pqrs
 */
enum Part {
  A, I,
  AB, AI, IJ,
  ABCD, ABIJ
};

class CoulombIntegrals {
  public:
    CoulombIntegrals(Chi *chi);
    CTF::Tensor<> &get(Part part);
    CTF::Tensor<> getSlice(Part part, int a, int b);
    /**
     * \deprecated  
     */
    void calculate_xycd(CTF::Tensor<> &xycd, int a, int b);

  protected:
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

    Chi *chi;
    CTF::Tensor<> *a, *i, *ai, *abij;
// NOTE: only for testing
    CTF::Tensor<> *abcd;
};

#endif

