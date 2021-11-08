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

#ifndef BASIS_SET_CORRECTION
#define BASIS_SET_CORRECTION

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  /**
   * \brief Caclulates basis set correction scheme called ps-ppl
   */
  class BasisSetCorrection: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(BasisSetCorrection)
    /**
     * \brief run routine 
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    bool run(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);
  
  };
}

#endif

