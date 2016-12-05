/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <util/StaticAssert.hpp>
#include <memory>

namespace tcc {
  // forward class declaration of interdependent expression types
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
    typedef F FieldType;
    virtual ~Expression() {
    }

    template <typename Rhs>
    Assignment<F> &operator =(Rhs &rhs) {
      static_assert(
        cc4s::StaticAssert<F>::False,
        "Only indexed tensors may be used as the left hand side of an assignment."
      );
    }

    /**
     * \brief Compile this tensor expression into a Operation.
     **/
    virtual std::shared_ptr<Operation<F>> compile(
      std::string const &lhsIndices
    ) = 0;
  };
}

#endif

