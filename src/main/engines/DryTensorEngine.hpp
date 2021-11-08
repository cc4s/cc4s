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

#ifndef DRY_TENSOR_ENGINE_DEFINED
#define DRY_TENSOR_ENGINE_DEFINED

#include <engines/DryMachineTensor.hpp>
#include <tcc/Costs.hpp>

// TODO: create object of this class for runtime arguments of tensor engine
// such as MPI communicators
namespace cc4s {
  /**
   * \brief Traits for inferring the respective DryMachineTensor types
   * from the respective tensor field types.
   * Tcc is given these traits upon compiling and execution.
   **/
  template <typename EmulatedTensorEngine>
  class DryTensorEngine {
  public:
    template <typename FieldType>
    using MachineTensor = DryMachineTensor<FieldType, EmulatedTensorEngine>;

    template <typename FieldType>
    static int64_t compareCosts(const Costs &l, const Costs &r) {
      return EmulatedTensorEngine::template compareCosts<FieldType>(l,r);
    }
  };


}

#endif

