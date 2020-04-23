/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CTF_MACHINE_TENSOR_DEFINED
#define CTF_MACHINE_TENSOR_DEFINED

#include <util/SharedPointer.hpp>
#include <ctf.hpp>
#include <string>

// TODO: specify MPI communicator when creating CtfTensorEngine
namespace cc4s {
  template <typename F,typename TE> class Tensor;
  class CtfTensorEngine;

  /**
   * \brief MachineTensor adapter for a CTF::Tensor
   **/
  template <typename F>
  class CtfMachineTensor {
  protected:
    class ProtectedToken {
    };

  public:
    typedef CTF::Tensor<F> T;
    typedef CtfTensorEngine TensorEngine;

    // constructors
    CtfMachineTensor(
      const std::vector<size_t> &lens,
      const std::string &name,
      const ProtectedToken &
    ):
      tensor(
        static_cast<int>(lens.size()),
        std::vector<int64_t>(lens.begin(), lens.end()).data(),
        std::vector<int>(0, lens.size()).data(),
        CTF::get_universe(), name.c_str()
      )
    {
    }

    // copy constructor from CTF tensor
    CtfMachineTensor(const T &t, const ProtectedToken &): tensor(t) {
    }

    ~CtfMachineTensor() {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      tensor.sum(
        alpha,
        A->tensor, aIndices.c_str(),
        beta,
        bIndices.c_str()
      );
    }

    // this[bIndices] = alpha * f(A[aIndices]) + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const F)> &f
    ) {
      tensor.sum(
        alpha,
        A->tensor, aIndices.c_str(),
        beta,
        bIndices.c_str(),
        CTF::Univar_Function<F>(f)
      );
    }

    template <typename G>
    void sum(
      G alpha,
      const Ptr<CtfMachineTensor<G>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const G)> &f
    ) {
      CTF::Transform<G,F>(
        std::function<void(const G, F &)>(
          [f,alpha,beta](const G x, F &y) { y = f(alpha*x) + beta*y; }
        )
      ) (
        A->tensor[aIndices.c_str()], tensor[bIndices.c_str()]
      );
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      const Ptr<CtfMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      tensor.contract(
        alpha,
        A->tensor, aIndices.c_str(),
        B->tensor, bIndices.c_str(),
        beta,
        cIndices.c_str()
      );
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contract(
      F alpha,
      const std::shared_ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      const std::shared_ptr<CtfMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      tensor.contract(
        alpha,
        A->tensor, aIndices.c_str(),
        B->tensor, bIndices.c_str(),
        beta,
        cIndices.c_str(),
        CTF::Bivar_Function<F>(g)
      );
    }

    void slice(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::vector<size_t> aBegins,
      const std::vector<size_t> aEnds,
      F beta,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      tensor.slice(
        std::vector<int>(begins.begin(), begins.end()).data(),
        std::vector<int>(ends.begin(), ends.end()).data(),
        beta,
        A->tensor,
        std::vector<int>(aBegins.begin(), aBegins.end()).data(),
        std::vector<int>(aEnds.begin(), aEnds.end()).data(),
        alpha
      );
    }

    // TODO: interfaces to be defined: permute, transform

    // read tensor elements to buffer
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
      tensor.read(
        elementsCount,
        reinterpret_cast<const int64_t *>(indexData),
        valueData
      );
    }

    // write tensor elements to buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
      tensor.write(
        elementsCount,
        reinterpret_cast<const int64_t *>(indexData),
        valueData
      );
    }

    std::vector<size_t> getLens() const {
      return std::vector<size_t>(tensor.lens, tensor.lens+tensor.order);
    }

    std::string getName() const {
      return std::string(tensor.get_name());
    }

    /**
     * \brief The adapted CTF tensor
     **/
    T tensor;

  protected:
    static Ptr<CtfMachineTensor<F>> create(const T &t) {
      return NEW(CtfMachineTensor<F>, t, ProtectedToken());
    }

    static Ptr<CtfMachineTensor<F>> create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return NEW(CtfMachineTensor<F>, lens, name, ProtectedToken());
    }

    friend class Tensor<F,CtfTensorEngine>;
  };
}

#endif
