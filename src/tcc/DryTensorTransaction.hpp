/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_TRANSACTION_DEFINED
#define DRY_TENSOR_TRANSACTION_DEFINED

#include <tcc/DryTensorExpression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F=double>
  class DryTensorTransaction {
  public:
    DryTensorTransaction(DryTensorExpression<F> &expression) {
      static_assert(
        StaticAssert<F>::False,
        "A transaction requires a single tensor assignment."
      );
    }
    DryTensorTransaction(
      DryTensorAssignment<F> &assignment
    ) {
      // compile time rules guarantee structure
      DryTensorContraction<F> *contraction(
        dynamic_cast<DryTensorContraction<F> *>(assignment.rhs)
      );
      if (contraction) buildContraction(assignment.lhs->indices, contraction);
      assignment.log();
    }
  protected:
    void buildContraction(
      std::string const &lhsIndices,
      DryTensorContraction<F> *contraction
    ) {
      LOG(0, "TCC") << "building contraction..." << std::endl;
      buildLocalContraction(
        lhsIndices,
        contraction->factors[0], contraction->factors[1]
      );
    }
    void buildLocalContraction(
      std::string const &lhsIndices,
      IndexedDryTensor<F> *a, IndexedDryTensor<F> *b
    ) {
      char contractedIndices[
        std::min(a->indices.length(), b->indices.length()) + 1
      ];
      char outerIndices[
        a->indices.length() + b->indices.length() + 1
      ];
      int j(0);
      for (int i(0); i < a->indices.length(); ++i) {
        const char index(a->indices[i]);
        if (
          b->indices.find(index) &&     // contracted index must occurr in both
          !lhsIndices.find(index) &&    // but not on the lhs
          !a->indices.find(index, i+1)  // and it must be unique
        ) {
          contractedIndices[j++] = index;
        }
      }
      contractedIndices[j] = 0;
      LOG(0, "TCC") << a->indices << "*" << b->indices << std::endl;
    }
  };
}

#endif

