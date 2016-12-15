/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/Tensor.hpp>
#include <tcc/Move.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/DryTensor.hpp>
#include <tcc/DryMachineTensor.hpp>

#include <vector>
#include <string>
#include <memory>

namespace tcc {
  template <typename F>
  class Tcc: public std::enable_shared_from_this<Tcc<F>> {
  protected:
    class ProtectedToken {
    };

  public:
    Tcc(
      const std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory_,
      const ProtectedToken &
    ): machineTensorFactory(machineTensorFactory_) {
    }

    static std::shared_ptr<Tcc<F>> create(
      const std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory_
    ) {
      return std::make_shared<Tcc<F>>(machineTensorFactory_, ProtectedToken());
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Symmetryies are not supported at the moment.
     * Note that tensor objects should only be created by the Tcc object
     * which specifies the environment the tensor lives in.
     */
    std::shared_ptr<Tensor<F>> createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::make_shared<Tensor<F>>(
        lens, name, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    std::shared_ptr<Tensor<F>> createTensor(
      const std::shared_ptr<MachineTensor<F>> &machineTensor
    ) {
      return std::make_shared<Tensor<F>>(
        machineTensor, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    std::shared_ptr<MachineTensor<F>> createMachineTensor(
      const std::shared_ptr<Tensor<F>> &tensor
    ) {
      return machineTensorFactory->createTensor(tensor->lens, tensor->name);
    }

  protected:
    std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory;
  };
}

#endif

