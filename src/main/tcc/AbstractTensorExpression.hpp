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

#ifndef TCC_ABSTRACT_TENSOR_EXPRESSION_DEFINED
#define TCC_ABSTRACT_TENSOR_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>

#include <tcc/TensorOperation.hpp>
#include <Exception.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class AbstractTensorExpression: public Expression<TE> {
  public:
    typedef F FieldType;

    virtual Ptr<TensorOperation<F,TE>> lhsCompile(
      const Ptr<TensorOperation<F,TE>> &rhsOperation
    ) {
      throw New<Exception>(
        "Writable expression expected on the "
        "left-hand-side of " + rhsOperation->getResult()->getName(),
        SourceLocation(rhsOperation->file, rhsOperation->line)
      );
    }
  };
}

#endif

