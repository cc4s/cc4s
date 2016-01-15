/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_DEFINED
#define COULOMB_INTEGRALS_DEFINED

#include <PerturbationTensor.hpp>
#include <Chi.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \deprecated This will no longer be used in the Data-Algorithm
   * design.
   */
  class CoulombIntegrals: public PerturbationTensor {
  public:
    CoulombIntegrals(Chi *chiReal, Chi *chiImag);
    virtual ~CoulombIntegrals();
    virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);

    /**
     * \brief Fetch all integrals in memory from chi.
     */
    void fetch();

    /**
     * \brief Fetch and return a slice of the Vxycd starting at x=a, y=b
     */
    CTF::Tensor<> getSlice(int a, int b);

    Chi *chiR, *chiI;
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
    CTF::Tensor<> *aijb;
    // NOTE: only allocated if storeV is enabled
    CTF::Tensor<> *aibc;
    CTF::Tensor<> *abci;
    CTF::Tensor<> *abcd;

  private:
    /**
     * \brief Fetch the given tensor from the chi tensors. The index map
     * contained in the tensor's name determine the indices taken from chi.
     */
    void fetch(CTF::Tensor<> &t, char const *indexMap);

    /**
     * \brief Add the given tensor to the tensor map v. The standard index
     * map is taken from the tensor's name.
     */
    void add(CTF::Tensor<> *t);

    // map of tensors: standard index map -> tensor
    std::map<std::string, CTF::Tensor<> *> v;
  };
}

#endif

