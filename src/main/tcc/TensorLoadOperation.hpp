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

#ifndef TCC_TENSOR_LOAD_OPERATION_DEFINED
#define TCC_TENSOR_LOAD_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <SharedPointer.hpp>

#include <string>

namespace cc4s {
  template <typename F, typename TE> class Tensor;

  template <typename F, typename TE>
  class TensorLoadOperation: public TensorOperation<F,TE> {
  public:
    TensorLoadOperation(
      const Ptr<Tensor<F,TE>> &source_,
      const Costs &costs_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      // the source tensor is also the result tensor, unless altered
      TensorOperation<F,TE>(
        source_,
        costs_,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      source(source_)
    {
    }

    void execute() override {
      if (source == this->getResult()) return;
      if (source->getVersion() > this->getResult()->getVersion()) {
        // move the data only if source and result tensors are different
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: move " <<
          this->getName() << " <<= " << source->getName() << std::endl;

        *this->getResult() = *source;
        this->updated();
        this->accountFlops();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() <<
          " up-to-date with " << source->getName() << std::endl;
      }
    }

    operator std::string () const override {
      if (source == this->getResult()) return std::string(*source);
      std::stringstream stream;
      stream << "load( " <<
        std::string(*source) <<
        " from " << source.get() <<
        " into " << this->result.get() << " )";
      return stream.str();
    }

  protected:
    Ptr<Tensor<F,TE>> source;

    static Ptr<TensorOperation<F,TE>> create(
      const Ptr<Tensor<F,TE>> &source,
      const Costs &costs,
      const Scope &scope
    ) {
      return New<TensorLoadOperation<F,TE>>(
        source, costs,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Tensor<F,TE>;
  };
}

#endif

