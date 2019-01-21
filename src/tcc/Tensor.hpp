/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/IndexedTensor.hpp>
//#include <tcc/MachineTensor.hpp>

#include <util/SharedPointer.hpp>

#include <cstdint>
#include <vector>
#include <string>

namespace tcc {
  /**
   * \brief 
   **/
  template <typename F, typename TE>
  class Tensor: public THISABLE(ESC(Tensor<F,TE>)) {
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };

  public:
    typedef typename TE::template MachineTensor<F> MT;
    /**
     * \brief Create a tcc tensor of yet unknown shape.
     * It must be first on the left-hand-side of an assignment.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::string &name_,
      const ProtectedToken &
    ): assumedShape(false), name(name_) {
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::vector<size_t> &lens_,
      const std::string &name_,
      const ProtectedToken &
    ): lens(lens_), assumedShape(true), name(name_) {
      // the machine tensor is not allocated initially
    }

    /**
     * \brief Create a tensor from a given machine tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const typename MT::T &unadaptedTensor_,
      const ProtectedToken &
    ): assumedShape(true) {
      auto mt(MT::create(unadaptedTensor_));
      lens = mt->getLens();
      name = mt->getName();
      machineTensor = mt;
    }

    static PTR(ESC(Tensor<F,TE>)) create(const std::string &name) {
      return NEW(ESC(Tensor<F,TE>), name, ProtectedToken());
    }

    static PTR(ESC(Tensor<F,TE>)) create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return NEW(ESC(Tensor<F,TE>), lens, name, ProtectedToken());
    }

    /**
     * \brief Create an empty tensor of identical shape as the given tensor.
     * The name, however, should differ.
     **/
    static PTR(ESC(Tensor<F,TE>)) create(
      const PTR(ESC(Tensor<F,TE>)) &tensor,
      const std::string &name
    ) {
      return create(tensor->lens, name);
    }

    static PTR(ESC(Tensor<F,TE>)) create(
      const typename MT::T &unadaptedTensor
    ) {
      return NEW(ESC(Tensor<F,TE>), unadaptedTensor, ProtectedToken());
    }

    void setName(const std::string &name_) {
      name = name_;
    }
    const std::string &getName() const {
      return name;
    }

    PTR(MT) getMachineTensor() {
      if (!machineTensor) {
        if (!assumedShape) {
          throw new EXCEPTION(
            "Tried to execute operation on tensor " + name +
            " before its shape has been assumed."
          );
        }
        // allocate the implementation specific machine tensor upon request
        machineTensor = MT::create(lens, name);
      }
      return machineTensor;
    }

    const std::vector<size_t> &getLens() const {
      return lens;
    }

    /**
     * \brief Returns the number of elements contained in this tensor.
     **/
    size_t getElementsCount() const {
      size_t elementsCount(1);
      for (auto len: lens) {
        elementsCount *= len;
      }
      return elementsCount;
    }

    /**
     * \brief Specify named indices of this tensor to be used in a
     * tensor expression. Indexed tensors are atomic types of tensor
     * expressions.
     **/
    PTR(ESC(IndexedTensor<F,TE>)) operator[](const std::string &indices) {
      return IndexedTensor<F,TE>::create(THIS, indices);
    }

    std::vector<size_t> lens;    
    bool assumedShape;

  protected:
    /**
     * \brief The tensor name.
     **/
    std::string name;

    PTR(MT) machineTensor;
  };
}

#endif

