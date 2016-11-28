/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_OPERATION_DEFINED
#define TENSOR_OPERATION_DEFINED

#include <tcc/TensorExpression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>
#include <algorithm>

namespace cc4s {
  template <typename F>
  class TensorOperation {
  public:
    virtual ~TensorOperation() {
    }

    virtual void execute() = 0;

    /**
     * \brief Number of tensor elements of storage required by the result.
     **/
    int64_t elementsCount;
    /**
     * \brief Maximum number of tensor elements of storage required
     * during the evaluation of this operation.
     **/
    int64_t maxMemoryCount;
    /**
     * \brief Number of tensor element multiplication required for
     * the evaluation of this operation.
     **/
    int64_t multiplicationsCount;
    /**
     * \brief Number of tensor elements additions required for
     * the evaluation of this operation.
     **/
    int64_t additionsCount;
  };

  template <typename F>
  TensorOperation<F> *compile(TensorExpression<F> &expression) {
    return expression.compile("");
  }
}

// include all known operation types
#include <tcc/TensorContractionOperation.hpp>
#include <tcc/TensorSumOperation.hpp>

#endif

