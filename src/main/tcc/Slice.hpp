/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SLICE_DEFINED
#define TCC_SLICE_DEFINED

#include <tcc/ClosedTensorExpression.hpp>

#include <tcc/SliceOperation.hpp>
#include <tcc/SliceIntoOperation.hpp>
#include <util/SharedPointer.hpp>
#include <vector>

namespace cc4s {
  /**
   * \brief 
   **/
  template <typename F, typename TE>
  class Slice: public ClosedTensorExpression<F,TE> {
  public:
    Slice(
      const PTR(ESC(ClosedTensorExpression<F,TE>)) &source_,
      const std::vector<size_t> &begins_,
      const std::vector<size_t> &ends_,
      const typename Expression<TE>::ProtectedToken &
    ): source(source_), begins(begins_), ends(ends_) {
    }

    static PTR(ESC(Slice<F,TE>)) create(
      const PTR(ESC(ClosedTensorExpression<F,TE>)) &source,
      const std::vector<size_t> &begins,
      const std::vector<size_t> &ends
    ) {
      return NEW(ESC(Slice<F,TE>),
        source, begins, ends, typename Expression<TE>::ProtectedToken()
      );
    }

    PTR(Operation<TE>) compile(Scope &scope) override {
      Scope sourceScope(scope.file, scope.line);
      auto sourceOperation(
        dynamicPtrCast<TensorOperation<F,TE>>(source->compile(sourceScope))
      );
      return SliceOperation<F,TE>::create(
        sourceOperation,
        Tensor<F,TE>::create(
          getLens(), sourceOperation->getResult()->getName()+"$"
        ),
        begins, ends,
        scope
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    ) override {
      auto lhsTensor(DYNAMIC_PTR_CAST(ESC(Tensor<F,TE>), source));
      if (!lhsTensor) {
        throw New<Exception>(
          "Expecting tensor for slice operation on left-hand-side.",
          SourceLocation(rhsOperation->file, rhsOperation->line)
        );
      }
      if (!rhsOperation->getResult()->assumedShape) {
        // create intermediate tensor as result tensor for rhsOperation
        auto resultLens(ends);
        for (size_t i(0); i < resultLens.size(); ++i) resultLens[i] -= begins[i];
        auto intermediateTensor(
          Tensor<F,TE>::create(resultLens, lhsTensor->getName() + "'")
        );
        rhsOperation->result = intermediateTensor;
        WARNING_LOCATION(
          SourceLocation(rhsOperation->file, rhsOperation->line)
        ) <<
          "updating parts of a slice on the left-hand-side "
          "currently requires the entire slice as intermediate tensor. "
          "This is less efficient than using write or read "
          "for the intended indices." << std::endl;
        if (rhsOperation->beta != F(1)) {
          WARNING_LOCATION(
            SourceLocation(rhsOperation->file, rhsOperation->line)
          ) <<
            "updating slice parts with operations other than += or -= "
            "may currently give wrong results as the entire slice is not "
            "read from the left-hand-side tensor before doing the update "
            "with the right-hand-side." << std::endl;
        }
      }
      auto sliceIntoOperation(
        SliceIntoOperation<F,TE>::create(
          rhsOperation,
          lhsTensor,
          begins, ends,
          Scope(rhsOperation->file, rhsOperation->line)
        )
      );
      // transfer beta from inner to outer operation
      sliceIntoOperation->beta = rhsOperation->beta;
      rhsOperation->beta = F(0);
      // TODO: transfer alpha in case of moves or contractions
      return sliceIntoOperation;
    }

    operator std::string () const override {
      return std::string(*source) + "( " +
        SliceOperation<F,TE>::coordinateString(begins) + "-" +
        SliceOperation<F,TE>::coordinateString(ends) + " )";
    }

  protected:
    std::vector<size_t> getLens() {
      std::vector<size_t> lens(ends);
      for (size_t i(0); i < lens.size(); ++i) {
        lens[i] -= begins[i];
      }
      return lens;
    }

    Ptr<ClosedTensorExpression<F,TE>> source;
    std::vector<size_t> begins, ends;
  };
}

#endif

