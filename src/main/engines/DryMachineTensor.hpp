/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <engines/DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <string>
#include <mpi.h>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  class DryTensorEngine;

  /**
   * \brief MachineTensor adapter for a DryTensor
   **/
  template <typename F>
  class DryMachineTensor {
  protected:
    class ProtectedToken {
    };

  public:
    typedef DryTensor<F> T;
    typedef DryTensorEngine TensorEngine;

    // constructors called by factory
    DryMachineTensor(
      const std::vector<size_t> &lens,
      const std::string &name,
      const ProtectedToken &
    ):
      tensor(
        lens.size(), lens.data(),
        std::vector<int>(lens.size(), 0).data()
      )
    {
      tensor.set_name(name);
    }

    // copy constructor from DryTensor
    DryMachineTensor(const T &t, const ProtectedToken &): tensor(t) {
    }

    ~DryMachineTensor() {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<DryMachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[bIndices] = f(alpha * A[aIndices]) + beta * this[bIndices]
    template <typename Domain>
    void sum(
      Domain alpha,
      const Ptr<DryMachineTensor<Domain>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const Domain)> &f
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<DryMachineTensor<F>> &A,
      const std::string &aIndices,
      const Ptr<DryMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<DryMachineTensor<F>> &A,
      const std::string &aIndices,
      const Ptr<DryMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    void slice(
      F alpha,
      const Ptr<DryMachineTensor<F>> &A,
      const std::vector<size_t> aBegins,
      const std::vector<size_t> aEnds,
      F beta,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      // allocate tensor for A assuming index reordering
      DryTensor<F> intermediateA(A->tensor, SOURCE_LOCATION);
      // allocate tensor for result assuming index reordering
      DryTensor<F> intermediateResult(this->tensor, SOURCE_LOCATION);
    }

    // TODO: interfaces to be defined: permute, transform

    // read tensor elements to buffer
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
    }

    void readToFile(MPI_File &file, const size_t offset = 0) {
    }

    // write tensor elements from buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
    }

    void writeFromFile(MPI_File &file, const size_t offset = 0) {
    }

    std::vector<size_t> getLens() const {
      return tensor.lens;
    }

    std::string getName() const {
      return tensor.name;
    }

    // adapted DryTensor
    T tensor;

    // TODO: protect
    // create adapater from given DryTensor
    static Ptr<DryMachineTensor<F>> create(const T &t) {
      return NEW(DryMachineTensor<F>, t, ProtectedToken());
    }

    // create adapter from shape and name
    static Ptr<DryMachineTensor<F>> create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return NEW(DryMachineTensor<F>, lens, name, ProtectedToken());
    }
  protected:
    friend class Tensor<F,DryTensorEngine>;
  };
}

#endif

