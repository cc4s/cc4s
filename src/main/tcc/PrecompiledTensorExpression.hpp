/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PRECOMPILED_TENSOR_EXPRESSION_DEFINED
#define PRECOMPILED_TENSOR_EXPRESSION_DEFINED

#include <tcc/ClosedTensorExpression.hpp>
#include <tcc/TensorOperation.hpp>

#include <tcc/IndexedTensorExpression.hpp>
#include <util/SharedPointer.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace cc4s {
  /**
   * \brief
   **/
  template <typename F, typename TE>
  class PrecompiledTensorExpression:
    public ClosedTensorExpression<F,TE>,
    public TensorOperation<F,TE>
  {
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared.
     **/
    class ProtectedToken {
    };
    /**
     * \brief precompiled operation to compute the result of this expression
     **/
    Ptr<TensorOperation<F,TE>> source;

  public:
    /**
     * \brief Create a tcc tensor of yet unknown shape.
     * It must be first on the left-hand-side of an assignment.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    PrecompiledTensorExpression(
      const std::string &lhsIndices,
      const Ptr<IndexedTensorExpression<F,TE>> &rhs,
      const std::string &name,
      const ProtectedToken &
    ): ClosedTensorExpression<F,TE>(),
    TensorOperation<F,TE>(
      nullptr, Costs(0), typename Operation<TE>::ProtectedToken()
    ),
    source(
      dynamicPtrCast<TensorOperation<F,TE>>(
        ((*Tensor<F,TE>::create(""))[lhsIndices] <<= rhs)->compile()
      )
    ) {
      Assert(source,
        "Expected operation returning a tensor as right-hand-side, got " +
        static_cast<std::string>(*rhs)
      );
      if (name != "") {
        source->getResult()->setName(name);
      }
    }

    static Ptr<PrecompiledTensorExpression<F,TE>> create(
      const std::string &lhsIndices,
      const Ptr<IndexedTensorExpression<F,TE>> &rhs,
      const std::string &name = ""
    ) {
      return New<PrecompiledTensorExpression<F,TE>>(
        lhsIndices, rhs, name, ProtectedToken()
      );
    }

    Ptr<Operation<TE>> compile(Scope &) override {
      // use this object for execution
      return this->template toPtr<Operation<TE>>();
    }

    Ptr<Tensor<F,TE>> getResult() override {
      return source->getResult();
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    operator std::string () const override {
      return source->getResult()->getName();
    }

    void execute() override {
      if (!source->getResult()->allocated()) source->execute();
    }
  };
}

#endif

