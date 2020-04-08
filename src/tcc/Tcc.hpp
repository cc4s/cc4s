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

// TODO: support slicing
// TODO: expression definitions with local index renaming
// TODO: binary function application
// TODO: support hard memory limit for costs
// TODO: looping over indices for memory reduction
// TODO: common subexpression optimization
// TODO: heuristics: limit number of simultaneously considered intermediates
// TODO: fix max memory assessment
// TODO: permutation and anti-permutation operator

namespace tcc {
  template <typename TensorEngine>
  class Tcc {
  protected:
    class ProtectedToken {
    };

  public:
    template <typename F=cc4s::real>
    static PTR(ESC(Tensor<F,TensorEngine>)) tensor(
      const std::vector<size_t> &lens, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(lens, name);
    }

    template <typename F=cc4s::real>
    static PTR(ESC(Tensor<F,TensorEngine>)) tensor(
      const PTR(ESC(Tensor<F,TensorEngine>)) &source, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(source->getLens(), name);
    }

    template <typename F=cc4s::real>
    static PTR(ESC(Tensor<F,TensorEngine>)) tensor(const std::string &name) {
      return Tensor<F,TensorEngine>::create(name);
    }

    static PTR(Sequence<TensorEngine>) nothing() {
      return NEW(Sequence<TensorEngine>);
    }
  };
}

#endif

