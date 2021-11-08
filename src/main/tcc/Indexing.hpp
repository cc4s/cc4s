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

#ifndef TCC_INDEXING_DEFINED
#define TCC_INDEXING_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/IndexingOperation.hpp>
#include <tcc/MoveOperation.hpp>
#include <util/SharedPointer.hpp>
#include <algorithm>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class TensorExpression;

  template <typename F, typename TE>
  class Indexing: public IndexedTensorExpression<F,TE> {
  public:
    /**
     * \brief Creates an expression with named indices from the closed tensor
     * source expression for further operations such as sums or contractions.
     * Not for direct invocation. Use the operator [] on
     * TensorExpression<F,TE> objects instead.
     **/
    Indexing(
      const Ptr<TensorExpression<F,TE>> &source_,
      const std::string &indices_,
      const typename Expression<TE>::ProtectedToken &
    ): source(source_), indices(indices_), canonicalIndices(indices_) {
      std::sort(canonicalIndices.begin(), canonicalIndices.end());
    }

    virtual ~Indexing() {
    }

    /**
     * \brief Creates an expression with named indices from a closed tensor
     * source expression for further operations such as sums or contractions.
     * \param[in] source The closed tensor expression to operate on.
     * \param[in] indices The index character string where each character
     * specifies the index name of the respective dimension index in the
     * source expression.
     **/
    static Ptr<Indexing<F,TE>> create(
      const Ptr<TensorExpression<F,TE>> &source,
      const std::string &indices
    ) {
      return New<Indexing<F,TE>>(
        source, indices, typename Expression<TE>::ProtectedToken()
      );
    }

    Ptr<Operation<TE>> compile(Scope &scope) override {
      Scope sourceScope(scope.file, scope.line);
      auto sourceOperation(
        dynamicPtrCast<TensorOperation<F,TE>>(source->compile(sourceScope))
      );
      return IndexingOperation<F,TE>::create(
        sourceOperation,
        sourceOperation->getResult(), // indexing operates on the closed tensor
        indices.c_str(),
        sourceOperation->costs,
        scope
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    Ptr<TensorOperation<F,TE>> lhsCompile(
      const Ptr<TensorOperation<F,TE>> &rhsOperation
    ) override {
      auto indexedRhs(
        dynamicPtrCast<IndexedTensorOperation<F,TE>>(rhsOperation)
      );
      ASSERT_LOCATION(
        indexedRhs,
        "Expecting indexed expression ([\"...\"]) on the right-hand-side of "
        "an assignment, got " + static_cast<std::string>(*rhsOperation),
        SourceLocation(rhsOperation->file, rhsOperation->line)
      );

      auto canonicalRhsIndices(indexedRhs->getResultIndices());
      std::sort(canonicalRhsIndices.begin(),canonicalRhsIndices.end());
      if (canonicalIndices == canonicalRhsIndices) {
        // rhs can be transposed to match lhs:
        // transposes the result of the indexed rhs tensor to match the index
        // order of this lhs indexed expression
        // determine shape of lhs
        size_t lenOfIndex[std::numeric_limits<char>::max()+1];
        // enter which indices correspond to which length on the rhs
        for (unsigned int i(0); i < indices.length(); ++i) {
          lenOfIndex[static_cast<unsigned int>(
            indexedRhs->getResultIndices()[i])
          ]= indexedRhs->getResult()->lens[i];
        }
        // transpose rhs result tensor
        for (unsigned int i(0); i < indices.length(); ++i) {
          indexedRhs->getResult()->lens[i] = lenOfIndex[
            static_cast<unsigned int>(indices[i])
          ];
        }
      } else {
        // otherwise: create new intermediate result tensor of unknown shape
        // it is expected to be overwritten by lhsCompile of lhs result tensor
        auto indexingResult(
          Tensor<F,TE>::create(
            indexedRhs->getResult()->getName() + "`"
          )
        );
        // use this result instead of the intermediate rhs result
        // in the last rhs operation together with the indices of the lhs
        indexedRhs->result = indexingResult;
        // TODO: assess costs of copy to larger lhs
      }
      // in any case, use the indices of the lhs in the last operation
      indexedRhs->resultIndices = indices;
      // continue lhs compiling of contained closed tensor expression
      return source->lhsCompile(indexedRhs);
    }

    void countIndices(Scope &scope) override {
      scope.add(indices);
    }

    operator std::string () const override {
      return std::string(*source) + "[" + indices + "]";
    }

    Ptr<TensorExpression<F,TE>> source;
    std::string indices, canonicalIndices;

  protected:
  };
}


#endif

