/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <util/SharedPointer.hpp>
#include <tcc/Operation.hpp>
#include <tcc/IndexCounts.hpp>

namespace tcc {
  template <typename F>
  class Tcc;

  template <typename F>
  class Expression: public THISABLE(Expression<F>) {
  public:
    /**
     * \brief Allow inferring the field type F given any of its expression
     * types.
     **/
    typedef F FieldType;

    virtual ~Expression() {
    }

    /**
     * \brief Compiles this expression and its subexpressions and returns
     * the resulting operation.
     **/
    virtual PTR(Operation<F>) compile(IndexCounts &indexCounts) {
      throw new EXCEPTION("Sequence (,) of move operation (<<=, +=, -=) expected.");
    }

    /**
     * \brief Count all indices occurring in this expression and its
     * subexpressions.
     **/
    virtual void countIndices(IndexCounts &indexCounts) = 0;

  protected:
    // FIXME: protect or discard
    /**
     * \brief The enclosing expression using the result of this expression or
     * nullptr if it is the outermost expression node.
     * Note that the object is not owned by the parent.
     **/
/*
    WEAK_PTR(Expression<F>) parent;
*/
    friend class Tcc<F>;

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

