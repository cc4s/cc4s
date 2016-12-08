/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_OPERATION_DEFINED
#define TENSOR_CONTRACTION_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>
#include <tcc/Costs.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>

#include <string>
#include <memory>
using std::shared_ptr;

namespace cc4s {
  template <typename F>
  class TensorContractionOperation: public TensorOperation<F> {
  public:
    TensorContractionOperation(
      const shared_ptr<TensorOperation<F>> &left_,
      const shared_ptr<TensorOperation<F>> &right_,
      DryTensor<F> *result_, const char *resultIndices_,
      Costs &contractionCosts
    ):
      TensorOperation<F>(left_->costs + right_->costs),
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
      if (result) delete result;
    }

    virtual void execute() {
      left->execute();
      right->execute();
      LOG(1, "TCC") << "executing " <<
        result->get_name() << "[" << resultIndices << "] = " <<
        left->getResult()->get_name() << "[" << left->getResultIndices() <<
        "] * " <<
        right->getResult()->get_name() << "[" << right->getResultIndices() <<
        "]" << std::endl;
    }

    virtual DryTensor<F> *getResult() {
      return result;
    }

    virtual std::string const &getResultIndices() {
      return resultIndices;
    }

  protected:
    shared_ptr<TensorOperation<F>> left;
    shared_ptr<TensorOperation<F>> right;

    DryTensor<F> *result;
    std::string resultIndices;

    friend class TensorContraction<F>;
  };
}

#endif

