/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <tcc/MachineTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

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
    // required by templates to infer corresponding Factory type
    typedef DryMachineTensorFactory<F> Factory;
    typedef DryTensor<F> Tensor;

    // constructors called by factory
    DryMachineTensor(
      const std::vector<int> &lens,
      const std::string &name,
      const ProtectedToken &
    ):
      tensor(
        static_cast<int>(lens.size()), lens.data(),
        std::vector<int>(lens.size(), 0).data()
      )
    {
      tensor.set_name(name);
    }

    // copy constructor from DryTensor, for compatibility
    DryMachineTensor(const Tensor &T, const ProtectedToken &): tensor(T) {
    }

    static std::shared_ptr<DryMachineTensor<F>> create(const Tensor &T) {
      return std::make_shared<DryMachineTensor<F>>(T, ProtectedToken());
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
        alpha << " * " << dryA->getName() << "[" << aIndices << "] + " <<
        beta << " * " << getName() << "[" << bIndices << "]" << std::endl;
    }

    // this[bIndices] = f(alpha * A[aIndices]) + beta * this[bIndices]
    virtual void move(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const F)> &f
    ) {
      std::shared_ptr<DryMachineTensor<F>> dryA(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(A)
      );
      if (!dryA) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      LOG(1, "TCC") << "move " <<
        getName() << "[" << bIndices << "] <<= f(" <<
        alpha << " * " << dryA->getName() << "[" << aIndices << "]) + " <<
        beta << " * " << getName() << "[" << bIndices << "]" << std::endl;
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    virtual void contract(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      const std::shared_ptr<tcc::MachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      std::shared_ptr<DryMachineTensor<F>> dryA(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(A)
      );
      std::shared_ptr<DryMachineTensor<F>> dryB(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(B)
      );
      if (!dryA || !dryB) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      LOG(1, "TCC") << "contract " <<
        getName() << "[" << cIndices << "] <<= " <<
        alpha << " * " << dryA->getName() << "[" << aIndices << "] * " <<
        dryB->getName() << "[" << bIndices << "] + " <<
        beta << " * " << getName() << "[" << cIndices << "]" << std::endl;
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    virtual void contract(
      F alpha,
      const std::shared_ptr<tcc::MachineTensor<F>> &A,
      const std::string &aIndices,
      const std::shared_ptr<tcc::MachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      std::shared_ptr<DryMachineTensor<F>> dryA(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(A)
      );
      std::shared_ptr<DryMachineTensor<F>> dryB(
        std::dynamic_pointer_cast<DryMachineTensor<F>>(B)
      );
      if (!dryA || !dryB) {
        throw new EXCEPTION("Passed machine tensor of wrong implementation.");
      }
      LOG(1, "TCC") << "contract " <<
        getName() << "[" << cIndices << "] <<= g(" <<
        alpha << " * " << dryA->getName() << "[" << aIndices << "], " <<
        dryB->getName() << "[" << bIndices << "]) + " <<
        beta << " * " << getName() << "[" << cIndices << "]" << std::endl;
    }

    // TODO: interfaces to be defined: slice, permute, transform

    virtual std::vector<int> getLens() const {
      return tensor.lens;
    }

    virtual std::string getName() const {
      return tensor.name;
    }

    // adapted DryTensor
    Tensor tensor;

    friend class DryMachineTensorFactory<F>;
  };

  template <typename F>
  class DryMachineTensorFactory: public tcc::MachineTensorFactory<F> {
  protected:
    class ProtectedToken {
    };

  public:
    DryMachineTensorFactory(const ProtectedToken &) {
    }

    virtual ~DryMachineTensorFactory() {
    }

    virtual std::shared_ptr<tcc::MachineTensor<F>> createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::shared_ptr<tcc::MachineTensor<F>>(
        std::make_shared<DryMachineTensor<F>>(
          lens, name, typename DryMachineTensor<F>::ProtectedToken()
        )
      );
    }

    static std::shared_ptr<DryMachineTensorFactory<F>> create(
    ) {
      return std::make_shared<DryMachineTensorFactory<F>>(
        ProtectedToken()
      );
    }
  };
}

#endif

