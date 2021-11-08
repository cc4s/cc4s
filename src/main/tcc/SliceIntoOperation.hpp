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

#ifndef TCC_SLICE_INTO_OPERATION_DEFINED
#define TCC_SLICE_INTO_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/SliceOperation.hpp>
#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Slice;

  template <typename F, typename TE>
  class SliceIntoOperation: public TensorOperation<F,TE> {
  public:
    SliceIntoOperation(
      const Ptr<TensorOperation<F,TE>> &source_,
      const Ptr<Tensor<F,TE>> &result_,
      const std::vector<size_t> begins_,
      const std::vector<size_t> ends_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, source_->costs,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      source(source_), begins(begins_), ends(ends_)
    {
    }

    void execute() override {
      source->execute();
      if (this->template isOlderThan<F>(source)) {
        // read from the entire result tensor of the source
        auto aEnds(source->getResult()->getLens());
        auto aBegins(std::vector<size_t>(aEnds.size()));

        std::stringstream beginsStream, endsStream, aBeginsStream, aEndsStream;
        for (auto d: begins) { beginsStream << " " << d; }
        for (auto d: ends) { endsStream << " " << d; }
        for (auto d: aBegins) { aBeginsStream << " " << d; }
        for (auto d: aEnds) { aEndsStream << " " << d; }

        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: slice " <<
          this->getName() << "(" <<
            beginsStream.str() << "," << endsStream.str() <<
          ") <<= " << this->alpha << " * " << source->getName() << "(" <<
            aBeginsStream.str() << "," << aEndsStream.str() << ") + " <<
          this->beta << " * " << this->getName() << "(" <<
            beginsStream.str() << "," << endsStream.str() <<
          ")" << std::endl;

        this->getResult()->getMachineTensor()->slice(
          F(1), source->getResult()->getMachineTensor(), aBegins, aEnds,
          F(this->beta), begins, ends
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() <<
          "up-to-date with " << source->getName() << std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return source->getLatestSourceVersion();
    }

    operator std::string () const override {
      return "SliceInto( " + std::string(*source) + ", " +
        SliceOperation<F,TE>::coordinateString(begins) + "-" +
        SliceOperation<F,TE>::coordinateString(ends) + " )";
    }

  protected:
    Ptr<TensorOperation<F,TE>> source;
    std::vector<size_t> begins, ends;

    static Ptr<SliceIntoOperation<F,TE>> create(
      const Ptr<TensorOperation<F,TE>> &source,
      const Ptr<Tensor<F,TE>> &result,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends,
      const Scope &scope
    ) {
      return New<SliceIntoOperation<F,TE>>(
        source, result, begins, ends,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

