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

#ifndef TCC_INDEXED_TENSOR_EXPRESSION_DEFINED
#define TCC_INDEXED_TENSOR_EXPRESSION_DEFINED

#include <tcc/AbstractTensorExpression.hpp>

namespace cc4s {
  /**
   * \brief Indexed tensor expressions are tensor valued expressions
   * whose dimensions are referred to by a specific index labels
   * for indentifying related indices in other indexed tensors.
   * Results of contractions are examples.
   * Indexing and closing converts between closed and indexed tensor
   * expressions. The transposition of a matrix tensor shared pointer A is
   * for instance done by: (*A)["ij"]^"ji"
   **/
  template <typename F, typename TE>
  class IndexedTensorExpression: public AbstractTensorExpression<F,TE> {
  public:
  };
}

#endif

