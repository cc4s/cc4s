/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_RECIPE_DEFINED
#define TENSOR_RECIPE_DEFINED

#include <tcc/ClosedTensorExpression.hpp>
#include <tcc/TensorOperation.hpp>

namespace cc4s {
  /**
   * \brief
   **/
  template <typename F, typename TE>
  class TensorRecipe:
    public ClosedTensorExpression<F,TE>,
    public TensorOperation<F,TE>
  {
  protected:
    /**
     * \brief Compiled operation to compute the result of this expression
     **/
    Ptr<Operation<TE>> recipe;
    bool evaluated;

  public:
    TensorRecipe(
      const Ptr<Tensor<F,TE>> &result_,
      const Ptr<Operation<TE>> &recipe_,
      const typename Expression<TE>::ProtectedToken &
    ): ClosedTensorExpression<F,TE>(
    ),
    TensorOperation<F,TE>(
      result_, Costs(0), typename Operation<TE>::ProtectedToken()
    ),
    recipe(
      recipe_
    ),
    evaluated(
      false
    ) {
    }

    static Ptr<TensorRecipe<F,TE>> create(
      const Ptr<Tensor<F,TE>> &result,
      const Ptr<Operation<TE>> &recipe
    ) {
      return New<TensorRecipe<F,TE>>(
        result, recipe, typename Expression<TE>::ProtectedToken()
      );
    }

    Ptr<Operation<TE>> compile(Scope &) override {
      // when using this object in another expression, use this object
      // for execution.
      // in the execute method it is decided whether to execute the recipe
      return this->template toPtr<Operation<TE>>();
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    operator std::string () const override {
      return this->getResult()->getName();
    }

    void execute(const size_t targetVersion) override {
      if (!evaluated) {
        recipe->execute();
        evaluated = true;
      }
    }
  };

  template <typename F, typename TE>
  inline Ptr<TensorRecipe<F,TE>> createTensorRecipe(
    const Ptr<Tensor<F,TE>> &result,
    const Ptr<Operation<TE>> &recipe
  ) {
    return TensorRecipe<F,TE>::create(result, recipe);
  }
}

#endif

