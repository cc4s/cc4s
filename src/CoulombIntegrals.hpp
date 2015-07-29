/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_DEFINED
#define COULOMB_INTEGRALS_DEFINED

#include "PerturbationTensor.hpp"
#include "Chi.hpp"
#include "Options.hpp"
#include <ctf.hpp>

class CoulombIntegrals: public PerturbationTensor {
  public:
    CoulombIntegrals(
      Chi *chiReal, Chi *chiImag, CTF::World *world, Options const *options
    );
    virtual ~CoulombIntegrals();

    virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);

    CTF::Tensor<> getSlice(int a, int b);

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

    Chi *chiR;
    Chi *chiI;
    CTF::Tensor<> *i;
    CTF::Tensor<> *a;
    CTF::Tensor<> *ij;
    CTF::Tensor<> *ia;
    CTF::Tensor<> *ai;
    CTF::Tensor<> *ab;
    CTF::Tensor<> *ijkl;
    CTF::Tensor<> *ijak;
    CTF::Tensor<> *aijk;
    CTF::Tensor<> *ijab;
    CTF::Tensor<> *abij;
    CTF::Tensor<> *aibj;
    // NOTE: only allocated if storeV is enabled
    CTF::Tensor<> *aibc;
    CTF::Tensor<> *abci;
    CTF::Tensor<> *abcd;

  private:
    /**
     * \brief Fetch all integrals in memory from chi.
     */
    void fetch();

    /**
     * \brief Fetch the given tensor from the chi tensors. The index map
     * contained in the tensor's name determine the indices taken from chi.
     */
    void fetch(CTF::Tensor<> *t, char const *indexMap);

    /**
     * \brief Add the given tensor to the tensor map v. The standard index
     * map is taken from the tensor's name.
     */
    void add(CTF::Tensor<> *t);

    // map of tensors: standard index map -> tensor
    std::map<std::string, CTF::Tensor<> *> v;
};

#endif

