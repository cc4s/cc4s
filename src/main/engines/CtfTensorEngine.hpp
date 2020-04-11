/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CTF_TENSOR_ENGINE_DEFINED
#define CTF_TENSOR_ENGINE_DEFINED

#include <engines/CtfMachineTensor.hpp>

// TODO: create object of this class for runtime arguments of tensor engine
// such as MPI communicators
namespace cc4s {
  /**
   * \brief Traits for inferring the respective DryMachineTensor types
   * from the respective tensor field types.
   * Tcc is given these traits upon compiling and execution.
   **/
  class CtfTensorEngine {
  public:
    template <typename FieldType>
    using MachineTensor = CtfMachineTensor<FieldType>;
  };
}

#endif

