/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_OPERATION_DEFINED
#define TENSOR_CONTRACTION_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>
#include <tcc/Costs.hpp>
#include <tcc/DryTensor.hpp>

#include <string>

namespace cc4s {
  template <typename F>
  class TensorContractionOperation: public TensorOperation<F> {
  public:
    TensorContractionOperation(
      TensorOperation<F> *left_,
      TensorOperation<F> *right_,
      DryTensor<F> *result_, std::string const &resultIndices_,
      Costs &contractionCosts
    ):
      TensorOperation<F>(left_.costs + right_.costs),
      left(left_), right(right_),
      result(result_),
      resultIndices(resultIndices_)
    {
      // so far, costs contains costs involved to get left and right factors
      // during contraction all elements of left,right and result are present
      contractionCosts.maxElementsCount =
        contractionCosts.elementsCount + this->costs.elementsCount;
      // the intermediate results are, however, no longer needed afterwards
      this->costs.elementsCount = 0;
      this->costs += contractionCosts;
    }

    virtual ~TensorContractionOperation() {
      // the suboperations and the intermediate result are dependent entities
      // unless they have been taken over by another entity
      if (left) delete left;
      if (right) delete right;
      if (result) delete result;
    }

    virtual void execute() {
      // TODO: call CTF functions or think of other binding mechanism
    }

    virtual DryTensor<F> *getResult() {
      return result;
    }

    virtual std::string const &getResultIndices() {
      return resultIndices;
    }

  protected:
    virtual void clearLeavingFetches() {
      left.clearLeavingFetches();
      delete left; left = nullptr;
      right.clearLeavingFetches();
      delete right; right = nullptr;
    }

    TensorOperation<F> *left;
    TensorOperation<F> *right;

    DryTensor<F> *result;
    std::string resultIndices;
  };
}

#endif

