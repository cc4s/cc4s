/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_ENGINE_DEFINED
#define DRY_TENSOR_ENGINE_DEFINED

#include <engines/DryMachineTensor.hpp>

// TODO: create object of this class for runtime arguments of tensor engine
// such as MPI communicators
namespace cc4s {
  /**
   * \brief Traits for inferring the respective DryMachineTensor types
   * from the respective tensor field types.
   * Tcc is given these traits upon compiling and execution.
   **/
  class DryTensorEngine {
  public:
    template <typename FieldType>
    using MachineTensor = DryMachineTensor<FieldType>;
  };
}

#endif

