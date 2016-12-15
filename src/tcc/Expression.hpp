/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <util/StaticAssert.hpp>
#include <memory>

namespace tcc {
  // forward class declaration of interdependent expression types and methods
  template <typename F=double>
  class Expression;

  template <typename F=double>
  class IndexedTensor;

  template <typename F=double>
  class Contraction;

  template <typename F=double>
  class Move;

  template <typename F=double>
  class Operation;

  template <typename F>
  class Expression {
  public:
    /**
     * \brief Allow inferring the field type F given any of its expression
     * types.
     **/
    typedef F FieldType;

    virtual ~Expression() {
    }

    // FIXME: protect
    /**
     * \brief The enclosing expression using the result of this expression or
     * nullptr if it is the outermost expression node.
     * Note that the object is not owned by the parent.
     **/
    std::weak_ptr<Expression<F>> parent;

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

