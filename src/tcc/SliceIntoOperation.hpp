/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SLICE_INTO_OPERATION_DEFINED
#define TCC_SLICE_INTO_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace tcc {
  template <typename F, typename TE> Slice;

  template <typename F, typename TE>
  class SliceIntoOperation: public TensorOperation<F,TE> {
  public:
    SliceIntoOperation(
      const PTR(ESC(TensorOperation<F,TE>)) &target_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const std::vector<size_t> begins_,
      const std::vector<size_t> ends_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<TE>(target->costs),
      source(source_), begins(begins_), ends(ends_)
    {
    }

    virtual ~SliceOperation() {
    }

    virtual void execute() {
      this->getResult()->getMachineTensor()->(
        F(1),
        getResult()->getMachineTensor(),
        // read from the entire result tensor
        std::vector<size_t>(this->getResult()->getLens()->size()),
        this->getResult()->getLens()
        F(0),
        begins, ends,
      );
      target->execute();
    }

  protected:
    PTR(ESC(TensorOperation<F,TE>)) source;
    std::vector<size_t> begins, ends;

    static PTR(ESC(SliceOperation<F,TE>)) create(
      const PTR(ESC(TensorOperation<F,TE>)) &source,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      return NEW(ESC(SliceOperation<F,TE>),
        source, begins, ends, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

