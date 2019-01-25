/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/ClosedTensorExpression.hpp>

#include <tcc/IndexedTensorExpression.hpp>
#include <util/SharedPointer.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace tcc {
  /**
   * \brief 
   **/
  template <typename F, typename TE>
  class Tensor: public ClosedTensorExpression<F,TE> {
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

    Tensor(
      const std::vector<size_t> &lens_,
      const std::string &name_,
      const bool assumedShape_,
      const ProtectedToken &
    ): lens(lens_), assumedShape(assumedShape_), name(name_) {
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
      return NEW(ESC(Tensor<F,TE>),
        tensor->lens, name, tensor->assumedShape, ProtectedToken()
      );
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

    std::vector<size_t> lens;    
    bool assumedShape;

    virtual PTR(Operation<TE>) compile(Scope &) {
      return TensorOperation<F,TE>::create(
        DYNAMIC_PTR_CAST(ESC(Tensor<F,TE>), THIS),
        Costs(getElementsCount())
      );
    }

    virtual PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    ) {
      if (!assumedShape) {
        // let this tensor assume the shape of the rhs result
        lens = rhsOperation->getResult()->getLens();
        assumedShape = true;
        // otherwise check shape
      } else if (rhsOperation->getResult()->getLens() != lens) {
        std::stringstream lhsShape;
        for (auto i: getLens()) { lhsShape << " " << i; }
        std::stringstream rhsShape;
        for (auto i: lens) { rhsShape << " " << i; }
        throw new EXCEPTION(
          "Shape of left-hand-side tensor " + getName() +
          " (" + lhsShape.str() + ") "
          " must match the shape of the result tensor " +
          rhsOperation->getResult()->getName() +
          " (" + rhsShape.str() + ")"
        );
      }
      // make the rhs operation directly operate on this tensor
      rhsOperation->result = DYNAMIC_PTR_CAST(ESC(Tensor<F,TE>), THIS);
      return rhsOperation;
    }

  protected:
    /**
     * \brief The tensor name.
     **/
    std::string name;

    PTR(MT) machineTensor;
  };
}

#endif

