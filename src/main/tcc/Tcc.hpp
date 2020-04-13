/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/Indexing.hpp>
#include <tcc/Move.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Map.hpp>
#include <tcc/Sequence.hpp>
#include <tcc/Slice.hpp>
#include <tcc/Tensor.hpp>

#include <util/SharedPointer.hpp>
#include <math/Real.hpp>

// TODO: support slicing
// TODO: expression definitions with local index renaming
// TODO: binary function application
// TODO: support hard memory limit for costs
// TODO: looping over indices for memory reduction
// TODO: common subexpression optimization
// TODO: heuristics: limit number of simultaneously considered intermediates
// TODO: fix max memory assessment
// TODO: permutation and anti-permutation operator

namespace cc4s {
  template <typename TensorEngine>
  class Tcc {
  protected:
    class ProtectedToken {
    };

  public:
    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(
      const std::vector<size_t> &lens, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(lens, name);
    }

    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(
      const Ptr<Tensor<F,TensorEngine>> &source, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(source->getLens(), name);
    }

    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(const std::string &name) {
      return Tensor<F,TensorEngine>::create(name);
    }

    static Ptr<Sequence<TensorEngine>> nothing() {
      return New<Sequence<TensorEngine>>();
    }
  };
}

#endif

