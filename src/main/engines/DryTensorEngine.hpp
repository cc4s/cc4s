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

