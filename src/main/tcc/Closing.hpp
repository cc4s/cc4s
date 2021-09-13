#ifndef TCC_CLOSING_DEFINED
#define TCC_CLOSING_DEFINED

#include <tcc/ClosedTensorExpression.hpp>

#include <tcc/IndexedTensorExpression.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE>
  class Closing: public ClosedTensorExpression<F,TE> {
  public:
    /**
     * \brief Creates an expression with ordered indices from the indexed tensor
     * source expression for further operations on the whole tensor such as
     * slicing.
     * Not for direct invocation. Use close on IndexedTensorExpressions instead.
     **/
    Closing(
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &source_,
      const std::string &indices_,
      const typename Expression<TE>::ProtectedToken &
    ): source(source_), indices(indices_) {
    }

    virtual ~Closing() {
    }

    /**
     * \brief Creates an expression with ordered indices from the indexed tensor
     * source expression for further operations on the whole tensor such as
     * slicing.
     * \param[in] source The indexed tensor whose indices should be closed
     * \param[in] indices The index character string where each character
     * specifies the index name in the source expression.
     **/
    static PTR(ESC(Closing<F,TE>)) create(
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &source,
      const std::string &indices
    ) {
      return NEW(ESC(Closing<F,TE>),
        source, indices, typename Expression<TE>::ProtectedToken()
      );
    }

    static PTR(ESC(Closing<F,TE>)) create(
      const PTR(ESC(Expression<F,TE>)) &source,
      const std::string &indices
    ) {
      static_assert(
        StaticAssert<RHS>::FALSE,
        "Only indexed tensor expressions can be closed."
      );
      return NEW(ESC(Closing<F,TE>),
        nullptr, "", typename Expression<TE>::ProtectedToken()
      );
    }

    virtual PTR(Operation<TE>) compile(Scope &scope) {
//      return ClosingOperation<F,TE>::create(result, resultIndices);
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    virtual void countIndices(Scope &scope) {
      // the closing indices are free indices and therefore must be counted
      scope.add(indices);
    }

    virtual operator std::string () const {
      return std::string(*source) + "^" + indices;
    }

  protected:
// TODO: to be done in operation
/*
    PTR(ESC(Tensor<F,TE>)) closedTensor(
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &source,
      const std::string &indices
    ) {
      if (indices.length() != source->getResultIndices().length()) {
        throw new EXCEPTION(
          "Number of indices for closing the indexed tensor " +
          source->getResult()->getName() +
          " must match the number of indices of the result tensor " +
          ->getResult()->getName()
        );
      }

      // determine shape of source
      std::vector<size_t> lens(source->getResult()->getLens().size());
      size_t lenOfIndex[std::numeric_limits<char>::max()+1];
      // enter which indices correspond to which length on the source
      for (unsigned int i(0); i < lhs->indices.length(); ++i) {
        lenOfIndex[static_cast<unsigned int>(operation->getResultIndices()[i])]=
          operation->result->getLens()[i];
      }
      // use for lhs
      for (unsigned int i(0); i < lhs->indices.length(); ++i) {
        lens[i] = lenOfIndex[static_cast<unsigned int>(lhs->indices[i])];
      }
      if (!lhs->tensor->assumedShape) {
        // assume
        lhs->tensor->lens = lens;
        lhs->tensor->assumedShape = true;
        // or check shape
      } else if (lhs->tensor->getLens() != lens) {
    }
*/
    PTR(ESC(IndexedTensorExpression<F,TE>)) source;
    std::string indices;
  };


  template <typename LHS>
  inline
  PTR(ESC(Closing<typename LHS::FieldType, typename LHS::TensorEngine>))
  operator ^(
    const PTR(LHS) &source, const std::string &indices
  ) {
    return Closing<typename LHS::FieldType, typename LHS::TensorEngine>::create(
      source, indices
    );
  }
}

#endif

