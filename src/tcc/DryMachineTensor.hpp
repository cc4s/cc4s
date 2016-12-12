/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <tcc/MachineTensor.hpp>
#include <string>
#include <memory>

namespace cc4s {
  template <typename F=double>
  class DryMachineTensorFactory;

  template <typename F>
  class DryMachineTensor: public tcc::MachineTensor<F> {
  protected:
    class ProtectedToken {
    };

  public:
    // constructors called by factory
    DryMachineTensor(
      const std::vector<int> &lens_,
      const std::string &name_,
      const ProtectedToken &
    ): lens(lens_), name(name_) {
    }

    virtual ~DryMachineTensor() {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    virtual void move(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      std::shared_ptr<DryMachineTensor<F>> dryA(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(A)
      );
      if (!dryA) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      LOG(1, "TCC") << "move " <<
        getName() << "[" << bIndices << "] <<= " <<
        dryA->getName() << "[" << aIndices << "]" << std::endl;
    }

    // this[bIndices] = alpha * f(A[aIndices]) + beta*this[bIndices]
    void move(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const F)> &f
    ) {
      std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A)
      );
      if (!ctfA) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      ctfTensor.sum(
        alpha,
        ctfA->ctfTensor, aIndices.c_str(),
        beta,
        bIndices.c_str(),
        CTF::Univar_Function<F>(f)
      );
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contraction(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      const std::shared_ptr<tcc::MachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      std::string &cIndices
    ) {
      std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A)
      );
      std::shared_ptr<CtfMachineTensor<F>> ctfB(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(B)
      );
      if (!ctfA || !ctfB) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      ctfTensor.contract(
        alpha,
        ctfA->ctfTensor, aIndices.c_str(),
        ctfB->ctfTensor, bIndices.c_str(),
        beta,
        cIndices.c_str()
      );
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contraction(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      const std::shared_ptr<tcc::MachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      std::shared_ptr<CtfMachineTensor<F>> ctfA(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(A)
      );
      std::shared_ptr<CtfMachineTensor<F>> ctfB(
        std::dynamic_pointer_cast<CtfMachineTensor<F>>(B)
      );
      if (!ctfA || !ctfB) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      ctfTensor.contract(
        alpha,
        ctfA->ctfTensor, aIndices.c_str(),
        ctfB->ctfTensor, bIndices.c_str(),
        beta,
        cIndices.c_str(),
        CTF::Bivar_Function<F>(g)
      );
    }

    // TODO: interfaces to be defined: slice, permute, transform

    virtual std::vector<int> getLens() const {
      return std::vector<int>(ctfTensor.lens, ctfTensor.lens+ctfTensor.order);
    }

    virtual std::string getName() const {
      return std::string(ctfTensor.get_name());
    }

  protected:
    /**
     * \brief The adapted CTF tensor
     **/
    CTF::Tensor<F> ctfTensor;

    friend class CtfMachineTensorFactory<F>;
  };

  template <typename F>
  class CtfMachineTensorFactory: public tcc::MachineTensorFactory<F> {
  protected:
    class ProtectedToken {
    };

  public:
    CtfMachineTensorFactory(
      CTF::World *world_, const ProtectedToken &
    ): world(world_) {
    }

    virtual ~CtfMachineTensorFactory() {
    }

    virtual std::shared_ptr<tcc::MachineTensor<F>> createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::shared_ptr<typename tcc::MachineTensor<F>>(
        std::make_shared<CtfMachineTensor<F>>(
          lens, name, world, typename CtfMachineTensor<F>::ProtectedToken()
        )
      );
    }

    static std::shared_ptr<CtfMachineTensorFactory<F>> create(
      CTF::World *world
    ) {
      return std::make_shared<CtfMachineTensorFactory<F>>(
        world, ProtectedToken()
      );
    }
  protected:
    CTF::World *world;
  };
}

#endif

