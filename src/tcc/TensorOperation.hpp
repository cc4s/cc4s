/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_OPERATION_DEFINED
#define TENSOR_OPERATION_DEFINED

#include <tcc/TensorExpression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>
#include <algorithm>

namespace cc4s {
  template <typename F=double>
  class TensorOperation {
  public:
    TensorOperation(TensorExpression<F> &expression) {
      static_assert(
        StaticAssert<F>::False,
        "A operation requires a single tensor assignment."
      );
    }
    TensorOperation(
      TensorAssignment<F> &assignment
    ) {
      // compile time rules guarantee exprected structure
      TensorContraction<F> *contraction(
        dynamic_cast<TensorContraction<F> *>(assignment.rhs)
      );
      if (contraction) buildContraction(assignment.lhs->indices, contraction);
      assignment.log();
    }

    virtual void execute() {
    }

  protected:
    void buildContraction(
      std::string const &lhsIndices,
      TensorContraction<F> *contraction
    ) {
      LOG(0, "TCC") << "building contraction..." << std::endl;
      buildLocalContraction(
        lhsIndices,
        contraction->factors[0], contraction->factors[1]
      );
    }
    void buildLocalContraction(
      std::string const &lhsIndices,
      IndexedTensor<F> *a, IndexedTensor<F> *b
    ) {
      char contractedIndices[
        std::min(a->indices.length(), b->indices.length()) + 1
      ];
      char outerIndices[
        a->indices.length() + b->indices.length() + 1
      ];
      int c(0), o(0);
      for (int i(0); i < a->indices.length(); ++i) {
        const char index(a->indices[i]);
        // go through unique indices of a
        if (a->indices.find(index, i+1) == std::string::npos) {
          if (
            b->indices.find(index) != std::string::npos &&  // in both
            lhsIndices.find(index) == std::string::npos     // but not on lhs
          ) {
            contractedIndices[c++] = index;
          } else {
            outerIndices[o++] = index;
          }
        }
      }
      contractedIndices[c] = 0;
      for (int i(0); i < b->indices.length(); ++i) {
        const char index(b->indices[i]);
        // go through unique indices of b
        if (b->indices.find(index, i+1) == std::string::npos) {
          if (
            std::find(contractedIndices, contractedIndices+c, index) ==
              contractedIndices+c
          ) {
            outerIndices[o++] = index;
          }
        }
      }
      outerIndices[o] = 0;
      LOG(0, "TCC") << *a << "*" << *b << ": " <<
        "contrated indices=\"" << contractedIndices << "\"" << ", " <<
        "outer indices=\"" << outerIndices << "\"" << std::endl;
    }


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
}

// include all known operation types
#include <tcc/TensorContractionOperation.hpp>
#include <tcc/TensorSumOperation.hpp>

#endif

