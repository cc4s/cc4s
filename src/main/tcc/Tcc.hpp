#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/Indexing.hpp>
#include <tcc/Move.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Map.hpp>
#include <tcc/Sequence.hpp>
#include <tcc/Slice.hpp>
#include <tcc/TensorRecipe.hpp>
#include <tcc/Tensor.hpp>
#include <util/SharedPointer.hpp>

// TODO: binary function application
// TODO: support hard memory limit for costs
// TODO: looping over indices for memory reduction
// TODO: automatic common subexpression optimization
// TODO: heuristics: limit number of simultaneously considered intermediates
// TODO: fix max memory assessment

/**
 * \breif Compiles one or more expressions for later execution.
 * \return Returns an Operation object whose execute() method will
 * execute the compiled expressions.
 **/
#define COMPILE(...) ((__VA_ARGS__)->compile(__FILE__, __LINE__))
/**
 * \breif Executes one or more expressions.
 **/
#define EXECUTE(...) ((__VA_ARGS__)->compile(__FILE__, __LINE__)->execute())
/**
 * \deprecated use PLAN() instead
 **/
#define COMPILE_RECIPE(RESULT, ...) \
  ((__VA_ARGS__)->compileRecipe(RESULT, __FILE__, __LINE__))
/**
 * \brief Plans the evaluation of the given expressions for later evaluation.
 * It is assumed that the execution of the given expressions updates the
 * tensor given as RESULT.
 * \return Returns a TensorRecipe object whose evaluate() method
 * ensures that the RESULT tensor is up-to-date with all source
 * tensors occurring in the given expressions. If any of the source
 * tensors is newer than the RESULT tensor the given expressions are executed
 * for updating the RESULT tensor.
 **/
#define PLAN(RESULT, ...) \
  ((__VA_ARGS__)->compileRecipe(RESULT, __FILE__, __LINE__))
/**
 * \brief Evaluates the given expression, assuming that they given
 * expressions update the tensor given as RESULT.
 * \return Returns the updated RESULT tensor.
 **/
#define EVALUATE(RESULT, ...) \
  ((__VA_ARGS__)->compileRecipe(RESULT, __FILE__, __LINE__)->evaluate())

namespace cc4s {
  template <typename TensorEngine>
  class Tcc {
  protected:
    class ProtectedToken {
    };

  public:
    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(
      const std::vector<size_t> &lens, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(lens, name);
    }

    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(
      const Ptr<Tensor<F,TensorEngine>> &source, const std::string &name
    ) {
      return Tensor<F,TensorEngine>::create(source->getLens(), name);
    }

    template <typename F=Real<>>
    static Ptr<Tensor<F,TensorEngine>> tensor(const std::string &name) {
      return Tensor<F,TensorEngine>::create(name);
    }

    static Ptr<Sequence<TensorEngine>> sequence() {
      return New<Sequence<TensorEngine>>();
    }
  };
}

#endif

