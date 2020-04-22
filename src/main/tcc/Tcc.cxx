/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/

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
}
