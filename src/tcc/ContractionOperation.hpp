/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_OPERATION_DEFINED
#define TCC_CONTRACTION_OPERATION_DEFINED

#include <tcc/TensorResultOperation.hpp>
#include <tcc/Costs.hpp>
#include <tcc/Tensor.hpp>

#include <string>
#include <memory>

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
      const std::shared_ptr<Operation<F>> &left_,
      const std::shared_ptr<Operation<F>> &right_,
      const std::shared_ptr<Tensor<F>> &result_,
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
    static std::shared_ptr<ContractionOperation<F>> create(
      const std::shared_ptr<Operation<F>> &left_,
      const std::shared_ptr<Operation<F>> &right_,
      const std::shared_ptr<Tensor<F>> &result_,
      const char *resultIndices_,
      const Costs &contractionCosts
    ) {
      return std::make_shared<ContractionOperation<F>>(
        left_, right_,
        result_, resultIndices_,
        contractionCosts,
        typename Operation<F>::ProtectedToken()
      );
    }

    std::shared_ptr<Operation<F>> left;
    std::shared_ptr<Operation<F>> right;

    friend class Tcc<F>;
  };
}

#endif

