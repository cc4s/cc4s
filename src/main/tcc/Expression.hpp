#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <tcc/Operation.hpp>
#include <tcc/Scope.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class TensorRecipe;
  template <typename F, typename TE>
  Ptr<TensorRecipe<F,TE>> createTensorRecipe(
    const Ptr<Tensor<F,TE>> &result, const Ptr<Operation<TE>> &recipe
  );

  template <typename TE>
  class Expression: public Thisable<Expression<TE>> {
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
    virtual PTR(Operation<TE>) compile(Scope &scope) {
      throw New<Exception>(
        "Sequence (,) of move operation (<<=, +=, -=) expected.",
        SourceLocation(scope.file, scope.line)
      );
    }

    virtual PTR(Operation<TE>) compile(
      const std::string &file, const size_t line
    ) {
      Scope scope(file, line);
      return this->compile(scope);
    }

    template <typename F>
    Ptr<TensorRecipe<F,TE>> compileRecipe(
      const Ptr<Tensor<F,TE>> &result,
      const std::string &file,
      const size_t line
    ) {
      return createTensorRecipe(result, compile(file, line));
    }

    // TODO: should be protected
    /**
     * \brief Count all indices occurring in this expression and its
     * subexpressions.
     **/
    virtual void countIndices(Scope &) = 0;

    virtual operator std::string () const = 0;

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

