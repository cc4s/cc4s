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

#ifndef TCC_INDEXING_OPERATION_DEFINED
#define TCC_INDEXING_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <SharedPointer.hpp>

#include <string>

namespace cc4s {
  template <typename F, typename TE> class Indexing;

  template <typename F, typename TE>
  class IndexingOperation: public IndexedTensorOperation<F,TE> {
  public:
    typedef F FieldType;

    IndexingOperation(
      const Ptr<TensorOperation<F,TE>> &source_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_,
        Costs(0), costs_,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      source(source_)
    {
    }

    void execute() override {
      source->execute();
    }

    operator std::string () const override {
      return std::string(*source) + "[" + this->resultIndices + "]";
    }

  protected:
    Ptr<TensorOperation<F,TE>> source;

    static Ptr<IndexingOperation<F,TE>> create(
      const Ptr<TensorOperation<F,TE>> &source,
      const Ptr<Tensor<F,TE>> &result,
      const char *resultIndices,
      const Costs &costs,
      const Scope &scope
    ) {
      return New<IndexingOperation<F,TE>>(
        source, result, resultIndices, costs,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Indexing<F,TE>;
  };
}

#endif

