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

#ifndef TCC_TENSOR_EXPRESSION_DEFINED
#define TCC_TENSOR_EXPRESSION_DEFINED

#include <tcc/TensorExpression.hpp>

#include <tcc/Indexing.hpp>
#include <tcc/Slice.hpp>
#include <util/SharedPointer.hpp>
#include <initializer_list>

namespace cc4s {
  /**
   * \brief Tensor expressions are expression
   * whose dimensions are referred to by their slot number rather than
   * by a specific index label. Tensors, tensor recipes and slices of tensors
   * are examples.
   * All tensor expressions can be used as source tensors within
   * other expressions. Tensors and slices can also be used as destination
   * tensors wihtin other expression, e.g. on the left-hand-side of a move
   * expression.
   * Access to the tensor data of tensor expressions is given via
   * the possibly modifying methods inspect() and evaluate(). Each returning
   * a Tensor object containing the requested data. inspect() only determines
   * size, type and unit, while evaluate() also determines all its data.
   **/
  template <typename F, typename TE>
  class TensorExpression: public AbstractTensorExpression<F,TE> {
  public:
    void countIndices(Scope &scope) override {
      // a tensor expression has no indices
    }

    /**
     * \brief Determines the size, the dimension types and the unit
     * of this tensor expression.
     * \return Returns a pointer to a tensor object containing the
     * requested information. The tensor may, however, contain no
     * data yet, indicated by isEvaluated() being false for the returned
     * tensor. In this case, subsequent calls to read() or write() fail.
     * Use evaluate() if you also want to access the tensor data.
     **/
    virtual Ptr<Tensor<F,TE>> inspect() {
      throw New<Exception>(
        "Only Tensor and TensorRecipe are supported for inspect().",
        SOURCE_LOCATION
      );
    }
    /**
     * \brief Evaluates all data of this tensor expression
     * including its size, dimension types and unit.
     * \return Returns a pointer to a tensor object containing the
     * requested information. The tensor is guaranteed to be
     * ready for read(), write(), getLengths() and more, indicated
     * by isShaped() and isEvaluated() being true for the returned tensor.
     **/
    virtual Ptr<Tensor<F,TE>> evaluate() {
      throw New<Exception>(
        "Only Tensor and TensorRecipe are supported for evaluate().",
        SOURCE_LOCATION
      );
    }

    /**
     * \brief Specify named indices of this tensor to be used in a
     * tensor expression. Indexed tensors are atomic types of tensor
     * expressions.
     **/
    Ptr<Indexing<F,TE>> operator [](
      const std::string &indices
    ) {
      return Indexing<F,TE>::create(
        this->template toPtr<TensorExpression<F,TE>>(), indices
      );
    }

    Ptr<Slice<F,TE>> operator ()(
      std::initializer_list<size_t> min,
      std::initializer_list<size_t> max
    ) {
      return Slice<F,TE>::create(
        this->template toPtr<TensorExpression<F,TE>>(),
        min, max
      );
    }
    Ptr<Slice<F,TE>> operator ()(
      const std::vector<size_t> &min,
      const std::vector<size_t> &max
    ) {
      return Slice<F,TE>::create(
        this->template toPtr<TensorExpression<F,TE>>(),
        min, max
      );
    }
  };
}

#endif

