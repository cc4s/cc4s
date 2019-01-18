/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <util/SharedPointer.hpp>
#include <tcc/Operation.hpp>
#include <tcc/IndexCounts.hpp>

namespace tcc {
  template <typename TE>
  class Expression: public THISABLE(Expression<TE>) {
  public:
    /**
     * \brief Allow inferring the tensor engine TE given any expression types.
     **/
    typedef TE TensorEngine;

    virtual ~Expression() {
    }

    /**
     * \brief Compiles this expression and its subexpressions and returns
     * the resulting operation.
     **/
    virtual PTR(Operation<TE>) compile(
      IndexCounts &indexCounts = IndexCounts()
    ) {
      throw new EXCEPTION("Sequence (,) of move operation (<<=, +=, -=) expected.");
    }

    // TODO: should be protected
    /**
     * \brief Count all indices occurring in this expression and its
     * subexpressions.
     **/
    virtual void countIndices(IndexCounts &indexCounts) = 0;

  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };
  };
}

#endif

