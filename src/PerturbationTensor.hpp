/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PERTURBATION_TENSOR_DEFINED
#define PERTURBATION_TENSOR_DEFINED

#include <ctf.hpp>
#include "Options.hpp"

class PerturbationTensor {
  public:
    PerturbationTensor(CTF::World *world, Options const *options);
    /**
     * \brief Returns the indexed tensor where stdIndexMap refers to
     * the standard names of the indices, irrespective of their use,
     * and indexMap refers to the index map used for the subsequent tensor
     * operation. E.g. T["klab"] results in a call to get("ijab", "klab").
     * This method is to be overriden by any deriving class.
     */
    virtual CTF::Idx_Tensor get(
      char const *stdIndexMap, char const *indexMap
    ) = 0;

    CTF::Idx_Tensor operator[](char const *indexMap);

    CTF::World *world;
    Options const *options;
};

#endif

