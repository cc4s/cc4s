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

#include <tcc/Tensor.hpp>
#include <Cc4s.hpp>

namespace cc4s {
  /**
   * \brief The version to be given the next tensor which is updated.
   * This number is incremented on each update to allow for comparing
   * the 'age' of tensor data.
   **/
  static size_t nextTensorVersion = 0;

  size_t getNextTensorVersion() {
    // all processes have to come together to ensure the same version is given
    // to each of them
    Cc4s::world->barrier();
    return ++nextTensorVersion;
  }

  std::map<std::string, Ptr<TensorDimension>> TensorDimension::dimensions;
}
