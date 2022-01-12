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

#ifndef TCC_EXPRESSION_DEFINED
#define TCC_EXPRESSION_DEFINED

#include <Object.hpp>

#include <tcc/Operation.hpp>
#include <tcc/Scope.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>
#include <Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class TensorRecipe;
  template <typename F, typename TE>
  Ptr<TensorRecipe<F,TE>> createTensorRecipe(
    const Ptr<Tensor<F,TE>> &result, const Ptr<Operation<TE>> &recipe
  );

  template <typename TE>
  class Expression: public Object {
  public:
    /**
     * \brief Allow inferring the tensor engine TE given any expression types.
     **/
    typedef TE TensorEngine;

    virtual ~Expression() {
    }

    /**
     * \brief Compiles this expression and its subexpressions and returns
     * the resulting operation.
     **/
    virtual Ptr<Operation<TE>> compile(Scope &scope) {
      throw New<Exception>(
        "Sequence (,) of move operation (<<=, +=, -=) expected.",
        SourceLocation(scope.file, scope.line)
      );
    }

    virtual Ptr<Operation<TE>> compile(
      const std::string &file, const size_t line
    ) {
      Scope scope(file, line);
      return this->compile(scope);
    }

    template <typename F>
    Ptr<TensorRecipe<F,TE>> compileRecipe(
      const Ptr<Tensor<F,TE>> &result,
      const std::string &file,
      const size_t line
    ) {
      return createTensorRecipe(result, compile(file, line));
    }

    // TODO: should be protected
    /**
     * \brief Count all indices occurring in this expression and its
     * subexpressions.
     **/
    virtual void countIndices(Scope &) = 0;

    virtual operator std::string () const = 0;

  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };
  };
}

#endif

