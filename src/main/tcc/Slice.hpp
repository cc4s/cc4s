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

#ifndef TCC_SLICE_DEFINED
#define TCC_SLICE_DEFINED

#include <tcc/TensorExpression.hpp>

#include <tcc/SliceOperation.hpp>
#include <tcc/SliceIntoOperation.hpp>
#include <util/SharedPointer.hpp>
#include <vector>

namespace cc4s {
  /**
   * \brief 
   **/
  template <typename F, typename TE>
  class Slice: public TensorExpression<F,TE> {
  public:
    Slice(
      const Ptr<TensorExpression<F,TE>> &source_,
      const std::vector<size_t> &begins_,
      const std::vector<size_t> &ends_,
      const typename Expression<TE>::ProtectedToken &
    ): source(source_), begins(begins_), ends(ends_) {
    }

    static Ptr<Slice<F,TE>> create(
      const Ptr<TensorExpression<F,TE>> &source,
      const std::vector<size_t> &begins,
      const std::vector<size_t> &ends
    ) {
      return New<Slice<F,TE>>(
        source, begins, ends, typename Expression<TE>::ProtectedToken()
      );
    }

    Ptr<Operation<TE>> compile(Scope &scope) override {
      Scope sourceScope(scope.file, scope.line);
      auto sourceOperation(
        dynamicPtrCast<TensorOperation<F,TE>>(source->compile(sourceScope))
      );
      return SliceOperation<F,TE>::create(
        sourceOperation,
        Tensor<F,TE>::create(
          getLens(), sourceOperation->getResult()->getName()+"$"
        ),
        begins, ends,
        scope
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    Ptr<TensorOperation<F,TE>> lhsCompile(
      const Ptr<TensorOperation<F,TE>> &rhsOperation
    ) override {
      auto lhsTensor(dynamicPtrCast<Tensor<F,TE>>(source));
      if (!lhsTensor) {
        throw New<Exception>(
          "Expecting tensor for slice operation on left-hand-side.",
          SourceLocation(rhsOperation->file, rhsOperation->line)
        );
      }
      if (!rhsOperation->getResult()->assumedShape) {
        // create intermediate tensor as result tensor for rhsOperation
        auto resultLens(ends);
        for (size_t i(0); i < resultLens.size(); ++i) resultLens[i] -= begins[i];
        auto intermediateTensor(
          Tensor<F,TE>::create(resultLens, lhsTensor->getName() + "'")
        );
        rhsOperation->result = intermediateTensor;
        WARNING_LOCATION(
          SourceLocation(rhsOperation->file, rhsOperation->line)
        ) <<
          "updating parts of a slice on the left-hand-side "
          "currently requires the entire slice as intermediate tensor. "
          "This is less efficient than using write or read "
          "for the intended indices." << std::endl;
        if (rhsOperation->beta != F(1)) {
          WARNING_LOCATION(
            SourceLocation(rhsOperation->file, rhsOperation->line)
          ) <<
            "updating slice parts with operations other than += or -= "
            "may currently give wrong results as the entire slice is not "
            "read from the left-hand-side tensor before doing the update "
            "with the right-hand-side." << std::endl;
        }
      }
      auto sliceIntoOperation(
        SliceIntoOperation<F,TE>::create(
          rhsOperation,
          lhsTensor,
          begins, ends,
          Scope(rhsOperation->file, rhsOperation->line)
        )
      );
      // transfer beta from inner to outer operation
      sliceIntoOperation->beta = rhsOperation->beta;
      rhsOperation->beta = F(0);
      // TODO: transfer alpha in case of moves or contractions
      return sliceIntoOperation;
    }

    operator std::string () const override {
      return std::string(*source) + "( " +
        SliceOperation<F,TE>::coordinateString(begins) + "-" +
        SliceOperation<F,TE>::coordinateString(ends) + " )";
    }

  protected:
    std::vector<size_t> getLens() {
      std::vector<size_t> lens(ends);
      for (size_t i(0); i < lens.size(); ++i) {
        lens[i] -= begins[i];
      }
      return lens;
    }

    Ptr<TensorExpression<F,TE>> source;
    std::vector<size_t> begins, ends;
  };
}

#endif

