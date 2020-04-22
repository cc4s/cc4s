/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_LOAD_OPERATION_DEFINED
#define TCC_TENSOR_LOAD_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Tensor;

  template <typename F, typename TE>
  class TensorLoadOperation: public TensorOperation<F,TE> {
  public:
    TensorLoadOperation(
      const PTR(ESC(Tensor<F,TE>)) &source_,
      const Costs &costs_,
      const typename Operation<TE>::ProtectedToken &
    ):
      // the source tensor is also the result tensor, unless altered
      TensorOperation<F,TE>(
        source_,
        costs_,
        typename Operation<TE>::ProtectedToken()
      ),
      source(source_)
    {
    }

    void execute(const size_t targetVersion) override {
      if (this->getResult() != source) {
        // move the data only if source and result tensors are different
        *this->getResult() = *source;
        LOG(2,"TCC") << "move " << this->getResult()->getName() << " <<= " <<
          source->getName() << std::endl;
      }
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Load( " <<
        std::string(*source) <<
        " from " << source.get() <<
        " into " << this->result.get() << " )";
      return stream.str();
    }

  protected:
    PTR(ESC(Tensor<F,TE>)) source;

    static PTR(ESC(TensorOperation<F,TE>)) create(
      const PTR(ESC(Tensor<F,TE>)) &source,
      const Costs &costs
    ) {
      return NEW(ESC(TensorLoadOperation<F,TE>),
        source, costs, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Tensor<F,TE>;
  };
}

#endif

