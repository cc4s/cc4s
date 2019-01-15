/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/IndexedTensor.hpp>
#include <tcc/Move.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Sequence.hpp>
#include <tcc/Tensor.hpp>

#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

#include <vector>
#include <string>
#include <memory>

// FIXME: tensor <<= sclar * tensor throws static exception
// TODO: *= operator with scalars
// TODO: heuristics: limit number of simultaneously considered intermediates
// TODO: fix max memory assessment
// TODO: mixed type tensor operations
// TODO: unary and binary function application
// TODO: expression definitions with local index renaming
// TODO: permutation and anti-permutation operator
// TODO: common subexpression optimization
// TODO: support hard memory limit for costs
// TODO: support slicing and looping over indices for memory reduction

namespace tcc {
  template <typename F>
  class Tcc: public THISABLE(Tcc<F>) {
  protected:
    class ProtectedToken {
    };

  public:
    Tcc(
      const PTR(MachineTensorFactory<F>) machineTensorFactory_,
      const ProtectedToken &
    ): machineTensorFactory(machineTensorFactory_) {
    }

    static PTR(Tcc<F>) create(
      const PTR(MachineTensorFactory<F>) machineTensorFactory_
    ) {
      return NEW(Tcc<F>, machineTensorFactory_, ProtectedToken());
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Symmetryies are not supported at the moment.
     * Note that tensor objects should only be created by the Tcc object
     * which specifies the environment the tensor lives in.
     **/
    PTR(Tensor<F>) createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return NEW(Tensor<F>,
        lens, name, THIS, typename Tensor<F>::ProtectedToken()
      );
    }

    /**
     * \brief Create an empty tensor of identical shape as the given tensor.
     * The name, however, should differ.
     **/
    PTR(Tensor<F>) createTensor(
      const PTR(Tensor<F>) &tensor,
      const std::string &name
    ) {
      return createTensor(tensor->lens, name);
    }

    /**
     * \brief Create a tcc tensor that uses a given machine tensor as its
     * machine tensor for representing the tensor data.
     **/
    PTR(Tensor<F>) createTensor(
      const PTR(MachineTensor<F>) &machineTensor
    ) {
      return NEW(Tensor<F>,
        machineTensor, THIS, typename Tensor<F>::ProtectedToken()
      );
    }

    PTR(MachineTensor<F>) createMachineTensor(
      const PTR(Tensor<F>) &tensor
    ) {
      return machineTensorFactory->createTensor(tensor->lens, tensor->name);
    }

    PTR(Sequence<F>) emptySequence() {
      return NEW(Sequence<F>);
    }

    PTR(Operation<F>) compile(const PTR(Expression<F>) &expression) {
      IndexCounts indexCounts;
      return expression->compile(indexCounts);
    }

  protected:
    PTR(MachineTensorFactory<F>) machineTensorFactory;
  };

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename Rhs>
  inline PTR(Move<typename Rhs::FieldType>) operator <<=(
    const PTR(IndexedTensor<typename Rhs::FieldType>) &lhs,
    const PTR(Rhs) &rhs
  ) {
    return Move<typename Rhs::FieldType>::create(lhs, rhs, 0);
  }

  template <typename Lhs, typename Rhs>
  inline PTR(Move<typename Lhs::FieldType>) operator <<=(
    const PTR(Lhs) &, const PTR(Rhs) &rhs
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename Lhs::FieldType, typename Rhs::FieldType
      >::EQUALS,
      "Move operations requires tensors of same type."
    );
    static_assert(
      cc4s::StaticAssert<Lhs>::FALSE,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return PTR(Move<typename Lhs::FieldType>)();
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->getName() << "[" << t.indices << "]";
  }
}

#endif

