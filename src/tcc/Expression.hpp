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
  class Assignment;

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

    /**
     * \brief Compile this tensor expression into a Operation.
     **/
    virtual std::shared_ptr<Operation<F>> compile(
      std::string const &lhsIndices
    ) = 0;

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

