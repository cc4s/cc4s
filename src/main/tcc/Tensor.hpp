/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/ClosedTensorExpression.hpp>

#include <tcc/TensorLoadOperation.hpp>
#include <util/SharedPointer.hpp>
#include <Cc4s.hpp>
#include <cstdint>
#include <vector>
#include <string>
#include <mpi.h>

namespace cc4s {
  size_t getNextTensorVersion();

  /**
   * \brief 
   **/
  template <typename F, typename TE>
  class Tensor: public ClosedTensorExpression<F,TE> {
  public:
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared aka New.
     **/
    class ProtectedToken {
    };

  public:
    // TODO: should be protected
    /**
     * \brief Dimensions of the tensor if assumedShape is true, otherwise empty.
     **/
    std::vector<size_t> lens;

    /**
     * \brief Whether this tensor has already well defined dimensions.
     * A tensor may be created without shape. Such a tensor assumes shape
     * from the result of an assignment.
     **/
    bool assumedShape;

  protected:
    /**
     * \brief The version of the data in this tensor. It is updated everytime
     * the tensor is updated.
     **/
    size_t version;

    /**
     * \brief The tensor name for in verbose info and error messages.
     **/
    std::string name;

    /**
     * \brief Pointer to the machine tensor, wrapping the underlying
     * tensor representation. The machine tensor is not allocated until
     * requested during execution or by manually calling getMachineTensor().
     **/
    typedef typename TE::template MachineTensor<F> MT;
    PTR(MT) machineTensor;

  public:
    /**
     * \brief Create a tcc tensor of yet unknown shape.
     * It must be first on the left-hand-side of an assignment.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::string &name_,
      const ProtectedToken &
    ): assumedShape(false), version(0), name(name_) {
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
    ):
      lens(lens_), assumedShape(true),
      version(0), name(name_)
    {
      // the machine tensor is not allocated initially
    }

    Tensor(
      const std::vector<size_t> &lens_,
      const std::string &name_,
      const bool assumedShape_,
      const ProtectedToken &
    ):
      lens(lens_), assumedShape(assumedShape_),
      version(0), name(name_)
    {
      // the machine tensor is not allocated initially
    }

    /**
     * \brief Create a tensor from a given machine tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const typename MT::T &unadaptedTensor_,
      const ProtectedToken &
    ): assumedShape(true), version(0) {
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
        ASSERT_LOCATION(assumedShape,
          "Tried to execute operation on tensor " + name +
          " before its shape has been assumed.",
          SOURCE_LOCATION
        );
        // allocate the implementation specific machine tensor upon request
        machineTensor = MT::create(lens, name);
        // wait until allocation is done on all processes
        Cc4s::world->barrier();
        LOG() << "Allocated tensor with " <<
          getElementsCount() << " elements" << std::endl;
      }
      return machineTensor;
    }

    bool allocated() {
      return machineTensor != nullptr;
    }

    size_t getVersion() const {
      return version;
    }

    void updated() {
      version = getNextTensorVersion();
    }

    // read tensor elements to buffer
    // each rank must specify its range of indices to read
    // behavior is undefined if ranges overlap
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
      getMachineTensor()->read(elementsCount, indexData, valueData);
    }
    F read(const size_t index = 0) {
      std::vector<size_t> indices(1);
      std::vector<F> values(1);
      if (Cc4s::world->getRank() == 0) {
        indices[0] = index;
        read(1, indices.data(), values.data());
      } else {
        read(0, indices.data(), values.data());
      }
      // broadcast value read on root
      Cc4s::world->broadcast(values);
      return values[0];
    }
    void readToFile(MPI_File &file, const size_t offset = 0) {
      getMachineTensor()->readToFile(file, offset);
    }

    // write tensor elements from buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
      getMachineTensor()->write(elementsCount, indexData, valueData);
      updated();
    }
    void write(F value, size_t index = 0) {
      write(1, &index, &value);
    }
    void writeFromFile(MPI_File &file, const size_t offset = 0) {
      getMachineTensor()->writeFromFile(file, offset);
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

    PTR(Operation<TE>) compile(Scope &scope) override {
      return TensorLoadOperation<F,TE>::create(
        this->template toPtr<Tensor<F,TE>>(),
        Costs(getElementsCount()),
        scope
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    ) override {
      // shape assuming:
      if (!assumedShape) {
        // lhs has not yet assumed shape
        if (!rhsOperation->getResult()->assumedShape) {
          // rhs has not yet assumed shape either
          ASSERT_LOCATION(false,
            "Neither left-hand-side tensor " + getName() +
            " nor right-hand-side result " +
            rhsOperation->getResult()->getName() +
            " have known shape. "
            "Try specifying left-hand-side tensor dimension manually.",
            SourceLocation(rhsOperation->file, rhsOperation->line)
          );
        } else {
          // let this lhs tensor assume the shape of the rhs result
          lens = rhsOperation->getResult()->getLens();
          assumedShape = true;
          // also assume name if empy
          if (name == "") name = rhsOperation->getResult()->getName();
        }
      } else {
        if (!rhsOperation->getResult()->assumedShape) {
          // let rhs tensor assume the shape of this lhs tensor
          rhsOperation->getResult()->lens = getLens();
          rhsOperation->getResult()->assumedShape = true;
        } else {
          // both have assumed shape: shapes must match
          if (rhsOperation->getResult()->getLens() != lens) {
            std::stringstream lhsShape;
            for (auto i: getLens()) { lhsShape << " " << i; }
            std::stringstream rhsShape;
            for (auto i: rhsOperation->getResult()->getLens()) { rhsShape << " " << i; }
            ASSERT_LOCATION(false,
              "Shape of left-hand-side tensor " + getName() +
              " (" + lhsShape.str() + ") "
              " must match the shape of the result tensor " +
              rhsOperation->getResult()->getName() +
              " (" + rhsShape.str() + ")",
              SourceLocation(rhsOperation->file, rhsOperation->line)
            );
          }
        }
      }

      // make the rhs operation directly operate on this tensor
      rhsOperation->result = this->template toPtr<Tensor<F,TE>>();
      return rhsOperation;
    }

    operator std::string () const override {
      return getName();
    }
  };
}

#endif

