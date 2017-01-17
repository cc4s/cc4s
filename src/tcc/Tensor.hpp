/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/IndexedTensor.hpp>
#include <tcc/MachineTensor.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace tcc {
  template <typename F=double>
  class Tcc;

  /**
   * \brief 
   **/
  template <typename F>
  class Tensor: public std::enable_shared_from_this<Tensor<F>> {
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };

  public:
    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::vector<int> &lens_,
      const std::string &name_,
      const std::shared_ptr<Tcc<F>> &tcc_,
      const ProtectedToken &
    ): lens(lens_), name(name_), tcc(tcc_) {
      // the machine tensor is not allocated initially
    }

    /**
     * \brief Create a tensor from a given machine tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::shared_ptr<MachineTensor<F>> &machineTensor_,
      const std::shared_ptr<Tcc<F>> &tcc_,
      const ProtectedToken &
    ):
      lens(machineTensor_->getLens()), name(machineTensor_->getName()),
      machineTensor(machineTensor_),
      tcc(tcc_)
    {
    }

    void setName(const std::string &name_) {
      name = name_;
    }
    const std::string &getName() const {
      return name;
    }

    std::shared_ptr<Tcc<F>> &getTcc() {
      return tcc;
    }

    template <class ActualMachineTensor=MachineTensor<F>>
    std::shared_ptr<ActualMachineTensor> getMachineTensor() {
      if (!machineTensor) {
        // allocate the implementation specific machine tensor upon request
        machineTensor = tcc->createMachineTensor(this->shared_from_this());
      }
      return std::dynamic_pointer_cast<ActualMachineTensor>(machineTensor);
    }

    /**
     * \brief Returns the number of elements contained in this tensor.
     **/
    int64_t getElementsCount() const {
      int64_t elementsCount(1);
      for (unsigned int i(0); i < lens.size(); ++i) {
        elementsCount *= lens[i];
      }
      return elementsCount;
    }

    /**
     * \brief Specify named indices of this tensor to be used in a
     * tensor expression. Indexed tensors are atomic types of tensor
     * expressions.
     **/
    std::shared_ptr<IndexedTensor<F>> operator[](const std::string &indices) {
      return IndexedTensor<F>::create(
        this->shared_from_this(), indices
      );
    }

    std::vector<int> lens;

  protected:
    /**
     * \brief The tensor name.
     **/
    std::string name;

    std::shared_ptr<MachineTensor<F>> machineTensor;

    /**
     * \brief The pointer to the tcc object that created this tensor.
     **/
    std::shared_ptr<Tcc<F>> tcc;

    friend class Tcc<F>;
  };
}

#endif

