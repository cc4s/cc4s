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

#ifndef CTF_TENSOR_ENGINE_DEFINED
#define CTF_TENSOR_ENGINE_DEFINED

#include <engines/CtfMachineTensor.hpp>
#include <tcc/Costs.hpp>
#include <math/MathFunctions.hpp>

// TODO: create object of this class for runtime arguments of tensor engine
// such as MPI communicators
namespace cc4s {
  /**
   * \brief Tensor engine using the Cyclops tensor framework.
   **/
  class CtfTensorEngine {
  public:
    template <typename FieldType>
    using MachineTensor = CtfMachineTensor<FieldType>;

    /**
     * \brief Returns -1,0, or +1 depending on whether the given costs
     * satisfy l<r, l=r, or l>r, respectively.
     * The comparison depends on the tensor engine's estimate.
     **/
    template <typename F>
    static int compareCosts(const Costs &l, const Costs &r) {
      Natural<128> lTotal(
        1000 * l.maxStorageCount +
        10 * l.accessCount +
        sizeof(F) / sizeof(real<F>(F(0))) * l.multiplicationsCount +
        l.additionsCount
      );
      Natural<128> rTotal(
        1000 * l.maxStorageCount +
        10 * l.accessCount +
        sizeof(F) / sizeof(real<F>(F(0))) * r.multiplicationsCount +
        r.additionsCount
      );
      return (rTotal < lTotal) - (lTotal < rTotal);
    }
  };
}

#endif

