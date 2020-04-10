/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SLICE_DEFINED
#define TCC_SLICE_DEFINED

#include <tcc/ClosedTensorExpression.hpp>

#include <tcc/SliceOperation.hpp>
#include <tcc/SliceIntoOperation.hpp>
#include <util/SharedPointer.hpp>
#include <vector>

namespace tcc {
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

    virtual PTR(Operation<TE>) compile(Scope &) {
      auto sourceOperation(
        DYNAMIC_PTR_CAST(ESC(TensorOperation<F,TE>), source->compile())
      );
      return SliceOperation<F,TE>::create(
        sourceOperation,
        Tensor<F,TE>::create(
          getLens(), sourceOperation->getResult()->getName()+"$"
        ),
        begins, ends
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    virtual PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    ) {
      auto targetTensor(DYNAMIC_PTR_CAST(ESC(Tensor<F,TE>), source));
      if (!targetTensor) {
        throw new EXCEPTION("Expecting tensor for slice operation on left-hand-side.");
      }
      auto sliceIntoOperation(
          SliceIntoOperation<F,TE>::create(
          rhsOperation,
          targetTensor,
          begins, ends
        )
      );
      // transfer beta from inner to outer operation
      sliceIntoOperation->beta = rhsOperation->beta;
      rhsOperation->beta = 0.0;
      // TODO: transfer alpha in case of moves or contractions
      return sliceIntoOperation;
    }

    virtual operator std::string () const {
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

    PTR(ESC(ClosedTensorExpression<F,TE>)) source;
    std::vector<size_t> begins, ends;
  };
}

#endif

