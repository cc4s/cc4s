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

#ifndef TCC_MOVE_OPERATION_DEFINED
#define TCC_MOVE_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class MoveOperation: public IndexedTensorOperation<F,TE> {
  public:
    /**
     * \brief Creates a move operation moving the results of
     * the right hand side operation into given tensor of the
     * left hand side after applying the function f.
     * The function f defaults to the identity operation.
     * Not intended for direct invocation. Use expression->compile()
     * to generate operations.
     **/
    MoveOperation(
      const Ptr<IndexedTensorOperation<F,TE>> &rhs_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      Costs moveCosts_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_,
        moveCosts_, rhs_->costs,
        file_, line_,
        typename Operation<TE>::ProtectedToken()
      ),
      rhs(rhs_)
    {
    }

    void execute() override {
      rhs->execute();
      if (this->template isOlderThan<F>(rhs)) {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: sum " <<
          this->getName() << " <<= " <<
          this->alpha << " * " << rhs->getName() << " + " <<
          this->beta << " * " << this->getName() << std::endl;

        this->getResult()->getMachineTensor()->sum(
          this->alpha,
          rhs->getResult()->getMachineTensor(), rhs->getResultIndices(),
          this->beta,
          this->resultIndices
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() << " up-to-date with " << rhs->getName() << std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return rhs->getLatestSourceVersion();
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "move( " << this->alpha << ", " <<
        std::string(*this->result) << ", " << std::string(*rhs) << ", " <<
        this->beta << " )";
      return stream.str();
    }

  protected:
    static Ptr<MoveOperation<F,TE>> create(
      const Ptr<IndexedTensorOperation<F,TE>> &rhs_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &moveCosts,
      const Scope &scope
    ) {
      return New<MoveOperation<F,TE>>(
        rhs_,
        result_, resultIndices_, moveCosts,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    Ptr<IndexedTensorOperation<F,TE>> rhs;
//    F alpha, beta;

    friend class Contraction<F,TE>;
    friend class Indexing<F,TE>;
  };
}

#endif

