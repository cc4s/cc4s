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

#ifndef TCC_INDEXED_TENSOR_OPERATION_DEFINED
#define TCC_INDEXED_TENSOR_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <SharedPointer.hpp>

#include <string>

namespace cc4s {
  template <typename F, typename TE> class Indexing;
  template <typename F, typename TE> class Move;
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class IndexedTensorOperation: public TensorOperation<F,TE> {
  public:
    typedef F FieldType;

    IndexedTensorOperation(
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &operationCosts_,
      const Costs &operandsCosts_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, operandsCosts_, file_, line_,
        typename Operation<TE>::ProtectedToken()
      ),
      resultIndices(resultIndices_)
    {
      // costs field constructed with costs of evaluating all operands
      // account for additional costs of this operation:
      this->costs.maxStorageCount = std::max(
        operationCosts_.maxStorageCount,  // more storage needed by result
        this->costs.maxStorageCount       // more storage needed by operands
      );
      this->costs.accessCount += operationCosts_.accessCount;
      this->costs.multiplicationsCount += operationCosts_.multiplicationsCount;
      this->costs.additionsCount += operationCosts_.additionsCount;
    }

    const std::string &getResultIndices() {
      return resultIndices;
    }

    std::string getName() const {
      return TensorOperation<F,TE>::getName() + "[" + resultIndices + "]";
    }

  protected:
    std::string resultIndices;

    friend class Indexing<F,TE>;
    friend class Move<F,TE>;
    friend class Contraction<F,TE>;
  };
}

#endif

