/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_OPERATION_DEFINED
#define TCC_CONTRACTION_OPERATION_DEFINED

#include <tcc/TensorResultOperation.hpp>
#include <tcc/Costs.hpp>
#include <tcc/Tensor.hpp>

#include <util/SharedPointer.hpp>

#include <string>

namespace tcc {
  template <typename F>
  class ContractionOperation: public TensorResultOperation<F> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    ContractionOperation(
      const PTR(Operation<F>) &left_,
      const PTR(Operation<F>) &right_,
      const PTR(Tensor<F>) &result_,
      const char *resultIndices_,
      Costs contractionCosts,
      const typename Operation<F>::ProtectedToken &
    ):
      TensorResultOperation<F>(
        result_, resultIndices_, left_->costs + right_->costs,
        typename Operation<F>::ProtectedToken()
      ),
      left(left_), right(right_)
    {
      // so far, costs contains costs involved to get left and right factors
      // during contraction all elements of left,right and result are present
      contractionCosts.maxElementsCount =
        contractionCosts.elementsCount + this->costs.elementsCount;
      // the intermediate results are, however, no longer needed afterwards
      this->costs.elementsCount = 0;
      this->costs += contractionCosts;
    }

    virtual ~ContractionOperation() {
    }

    virtual void execute() {
      left->execute();
      right->execute();
      this->getResult()->getMachineTensor()->contract(
        this->alpha,
        left->getResult()->getMachineTensor(), left->getResultIndices(),
        right->getResult()->getMachineTensor(), right->getResultIndices(),
        this->beta,
        this->getResultIndices()
      );
    }

  protected:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     **/
    static PTR(ContractionOperation<F>) create(
      const PTR(Operation<F>) &left_,
      const PTR(Operation<F>) &right_,
      const PTR(Tensor<F>) &result_,
      const char *resultIndices_,
      const Costs &contractionCosts
    ) {
      return NEW(ContractionOperation<F>,
        left_, right_,
        result_, resultIndices_,
        contractionCosts,
        typename Operation<F>::ProtectedToken()
      );
    }

    PTR(Operation<F>) left;
    PTR(Operation<F>) right;

    friend class Tcc<F>;
  };
}

#endif

