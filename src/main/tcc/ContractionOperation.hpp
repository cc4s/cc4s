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

#ifndef TCC_CONTRACTION_OPERATION_DEFINED
#define TCC_CONTRACTION_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <tcc/Tensor.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class ContractionOperation: public IndexedTensorOperation<F,TE> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    ContractionOperation(
      const Ptr<IndexedTensorOperation<F,TE>> &left_,
      const Ptr<IndexedTensorOperation<F,TE>> &right_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      Costs contractionCosts,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_,
        contractionCosts, left_->costs + right_->costs,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      left(left_), right(right_)
    {
    }

    void execute() override {
      left->execute();
      right->execute();
      if (
        this->template isOlderThan<F>(left) ||
        this->template isOlderThan<F>(right)
      ) {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: contraction " <<
          this->getName() << " <<= " <<
          this->alpha << " * " <<
          left->getName() << " * " << right->getName() << " + " <<
          this->beta << " * " << this->getName() << std::endl;

        this->getResult()->getMachineTensor()->contract(
          this->alpha,
          left->getResult()->getMachineTensor(), left->getResultIndices(),
          right->getResult()->getMachineTensor(), right->getResultIndices(),
          this->beta,
          this->getResultIndices()
        );
        this->updated();
        this->accountFlops();
      } else {
          LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getResult()->getName() << " up-to-date with " <<
          "(" << left->getName() << ", " << right->getName() << ")" <<
          std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return std::max(
        left->getLatestSourceVersion(),
        right->getLatestSourceVersion()
      );
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "contraction( " << this->alpha << ", " <<
        std::string(*left) << ", " << std::string(*right) << ", " <<
        this->beta << " )";
      return stream.str();
    }

  protected:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     **/
    static Ptr<ContractionOperation<F,TE>> create(
      const Ptr<IndexedTensorOperation<F,TE>> &left_,
      const Ptr<IndexedTensorOperation<F,TE>> &right_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &contractionCosts,
      const Scope &scope
    ) {
      return New<ContractionOperation<F,TE>>(
        left_, right_,
        result_, resultIndices_,
        contractionCosts,
        scope.file, scope.line,
        typename Operation<TE>::ProtectedToken()
      );
    }

    Ptr<IndexedTensorOperation<F,TE>> left;
    Ptr<IndexedTensorOperation<F,TE>> right;

    friend class Contraction<F,TE>;
  };
}

#endif

