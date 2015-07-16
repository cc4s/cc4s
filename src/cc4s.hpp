/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CC4S_DEFINED
#define CC4S_DEFINED

#include <ctf.hpp>
#include "CoulombIntegrals.hpp"
#include "Chi.hpp"

class cc4s {
  public:
    static CTF::World *world;
    static int rank, np, no, nv, nG, niter;
    static CoulombIntegrals *V;
    static CTF::Tensor<> *Tabij, *Tai;
    static Chi *chi;
    static bool profile;

    static void startup();
    static void cleanup();
    static void run();

  protected:
    static void iterateAmplitudes();
    static void add_Vxyef_T21efij(CTF::Tensor<> &Zabij, CTF::Tensor<> &T21);

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts a total index of a tensor entry into the positions
     * along each dimension.
     */
    // FIXME: cannot handle out of bound indicies
    static void from_index(CTF::Tensor<> const &t, int64_t index, int *pos) {
      for (int d(0); d < t.order; ++d) {
        pos[d] = index % (int64_t)t.lens[d];
        index /= t.lens[d];
      }
    }

    // TODO: check if the cft framework offers such a function natively
    /**
     * \brief Converts given positions along each dimension into a total index.
     */
    static int64_t to_index(int dim, int const *lens, int const *pos) {
      int64_t index(0);
      for (int d(dim-1); d >= 0; --d) {
        index *= lens[d];
        index += pos[d];
      }
      return index;
    }
  friend class Chi;
  friend class CoulombIntegrals;
};

#endif

