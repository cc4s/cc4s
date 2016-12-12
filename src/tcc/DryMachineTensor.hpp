/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <tcc/MachineTensor.hpp>
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
      return lens;
    }

    virtual std::string getName() const {
      return name;
    }

  protected:
    std::vector<int> lens;
    std::string name;

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

