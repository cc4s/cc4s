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

#ifndef UEG_VERTEX_GENERATOR_DEFINED
#define UEG_VERTEX_GENERATOR_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/TensorExpression.hpp>
#include <TensorSet.hpp>

using ivec  = array<int,3>;
using dvec  = array<double,4>;


namespace cc4s {
  class UegVertexGenerator: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UegVertexGenerator)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    Ptr<MapNode> run(const Ptr<MapNode> &arguments);

    double evalMadelung(double volume);
    double Vijji(const dvec a, const dvec b, const double v);

    int halfGrid;
    size_t No, Nv, NF;
    double rs, madelung;
  };
}

#endif

