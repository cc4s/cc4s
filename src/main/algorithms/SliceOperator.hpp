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

#ifndef SLICE_OPERATOR_DEFINED
#define SLICE_OPERATOR_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/TensorExpression.hpp>
#include <TensorSet.hpp>

namespace cc4s {
  class SliceOperator: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(SliceOperator)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);
    template <typename F, typename TE>
    void slice(
      const Ptr<TensorExpression<F,TE>> &tensor,
      const std::string &parts,
      const Ptr<TensorSet<F,TE>> &slices
    );

    std::vector<Natural<>> dims;
    Natural<> No, Nv;
  };
}

#endif

