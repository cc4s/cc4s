/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TENSOR_RECIPE_DEFINED
#define TENSOR_RECIPE_DEFINED

#include <tcc/TensorExpression.hpp>
#include <tcc/TensorOperation.hpp>

namespace cc4s {
  /**
   * \brief
   **/
  template <typename F, typename TE>
  class TensorRecipe:
    public TensorExpression<F,TE>,
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
    ): TensorExpression<F,TE>(
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

    Ptr<Tensor<F,TE>> inspect() override {
      return this->getResult();
    }
    Ptr<Tensor<F,TE>> evaluate() override {
      execute();
      return this->getResult();
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
          this->getName() << " up-to-date with all sources. " <<
          "Operations: " << this->getFloatingPointOperations() << std::endl;
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

