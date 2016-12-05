/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CTF_TENSOR_ADAPTER_DEFINED
#define CTF_TENSOR_ADAPTER_DEFINED

#include <ctf.hpp>
#include <string>
#include <memory>

namespace cc4s {
  template <typename F>
  class CtfTensorFactory {
  public:
    // zero constructor, no symmetries supported 
    std::shared_ptr<Tensor<F>> create(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::make_shared<F>(lens, name, world)
    }

  protected:
    CTF::World *world;
  };

  template <typename F>
  class CtfTensorAdapter {
    // constructors called by factory
    Tensor(
      const std::vector<int> &lens,
      const std::string &name,
      CTF::World *world
    ):
      ctfTensor(
        lens.size(), lens.data(), std::vector<int>(0, lens.size()),
        *world, name
      )
    {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    void move(
      F alpha,
      Tensor<F> &A, const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      ctfTensor.sum(
        alpha,
        A.ctfTensor, aIndices.c_str(),
        beta,
        bIndices.c_str()
      );
    }

    // this[bIndices] = alpha * f(A[aIndices]) + beta*this[bIndices]
    void move(
      F alpha,
      Tensor<F> &A, const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const F)> &f
    ) {
      ctfTensor.sum(
        alpha
        A.ctfTensor, aIndices.c_str(),
        beta,
        bIndices.c_str(),
        CTF::Univar_Function<F>(f)
      );
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contraction(
      F alpha,
      Tensor &A, const std::string &aIndices,
      Tensor &B, const std::string &bIndices,
      F beta,
      std::string &cIndices
    ) {
      ctfTensor.contract(
        alpha,
        A.ctfTensor, aIndices.c_str(),
        B.ctfTensor, bIndices.c_str(),
        beta,
        cIndices.c_str()
      );
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contraction(
      F alpha,
      Tensor &A, const std::string &aIndices,
      Tensor &B, const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      ctfTensor.contract(
        alpha,
        A.ctfTensor, aIndices.c_str(),
        B.ctfTensor, bIndices.c_str(),
        beta,
        cIndices.c_str(),
        CTF::Bivar_Function(g)
      );
    }

    // TODO: interfaces to be defined: slice, permute, transform

    CTF::Tensor<F> ctfTensor;
  };
}

