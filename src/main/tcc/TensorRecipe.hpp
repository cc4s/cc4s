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

  public:
    TensorRecipe(
      const Ptr<Tensor<F,TE>> &result_,
      const Ptr<Operation<TE>> &recipe_,
      const typename Expression<TE>::ProtectedToken &
    ): ClosedTensorExpression<F,TE>(
    ),
    TensorOperation<F,TE>(
      result_, recipe_->costs,
      recipe_->file, recipe_->line, typename Operation<TE>::ProtectedToken()
    ),
    recipe(
      recipe_
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

    void execute() override {
      if (getLatestSourceVersion() >= this->result->getVersion()) {
        recipe->execute();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() << " up-to-date with all sources." << std::endl;
      }
    }

    // return latest version of any tensor of the recipe
    size_t getLatestSourceVersion() override {
      return recipe->getLatestSourceVersion();
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

