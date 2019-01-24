/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <tcc/DryTensor.hpp>

#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  class DryEngine;

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
    typedef DryEngine TensorEngine;

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
    void move(
      F alpha,
      const PTR(DryMachineTensor<F>) &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      LOG(1, "TCC") << "move " <<
        getName() << "[" << bIndices << "] <<= " <<
        alpha << " * " << A->getName() << "[" << aIndices << "] + " <<
        beta << " * " << getName() << "[" << bIndices << "]" << std::endl;
    }

    // this[bIndices] = f(alpha * A[aIndices]) + beta * this[bIndices]
    template <typename Domain>
    void move(
      Domain alpha,
      const PTR(DryMachineTensor<Domain>) &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const Domain)> &f
    ) {
      LOG(1, "TCC") << "move " <<
        getName() << "[" << bIndices << "] <<= f(" <<
        alpha << " * " << A->getName() << "[" << aIndices << "]) + " <<
        beta << " * " << getName() << "[" << bIndices << "]" << std::endl;
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const PTR(DryMachineTensor<F>) &A,
      const std::string &aIndices,
      const PTR(DryMachineTensor<F>) &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      LOG(1, "TCC") << "contract " <<
        getName() << "[" << cIndices << "] <<= " <<
        alpha << " * " << A->getName() << "[" << aIndices << "] * " <<
        B->getName() << "[" << bIndices << "] + " <<
        beta << " * " << getName() << "[" << cIndices << "]" << std::endl;
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contract(
      F alpha,
      const PTR(DryMachineTensor<F>) &A,
      const std::string &aIndices,
      const PTR(DryMachineTensor<F>) &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      LOG(1, "TCC") << "contract " <<
        getName() << "[" << cIndices << "] <<= g(" <<
        alpha << " * " << A->getName() << "[" << aIndices << "], " <<
        B->getName() << "[" << bIndices << "]) + " <<
        beta << " * " << getName() << "[" << cIndices << "]" << std::endl;
    }

    void slice(
      F alpha,
      const PTR(DryMachineTensor<F>) &A,
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
/*
      LOG(2, "TCC") << "slice " <<
        getName() << "[" << begins << "," ends << ") <<= " <<
        alpha << " * " <<
        ctfA->getName() << "[" << aBegins << "," << aEnds << "] + " <<
        beta << " * " <<
        getName() << "[" << begins << "," << ends << "]" << std::endl;
*/
    }

    // TODO: interfaces to be defined: permute, transform

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
    static PTR(DryMachineTensor<F>) create(const T &t) {
      return NEW(DryMachineTensor<F>, t, ProtectedToken());
    }

    // create adapter from shape and name
    static PTR(DryMachineTensor<F>) create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return NEW(DryMachineTensor<F>, lens, name, ProtectedToken());
    }
  protected:
    friend class Tensor<F,DryEngine>;
  };

  /**
   * \brief Traits for inferring the respective DryMachineTensor types
   * from the respective tensor field types.
   * Tcc is given these traits upon compiling and execution.
   **/
  class DryEngine {
  public:
    template <typename FieldType>
    using MachineTensor = DryMachineTensor<FieldType>;
  };
}

#endif

