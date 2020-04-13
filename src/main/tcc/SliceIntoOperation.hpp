/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SLICE_INTO_OPERATION_DEFINED
#define TCC_SLICE_INTO_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/SliceOperation.hpp>
#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Slice;

  template <typename F, typename TE>
  class SliceIntoOperation: public TensorOperation<F,TE> {
  public:
    SliceIntoOperation(
      const PTR(ESC(TensorOperation<F,TE>)) &source_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const std::vector<size_t> begins_,
      const std::vector<size_t> ends_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, source_->costs, typename Operation<TE>::ProtectedToken()
      ),
      source(source_), begins(begins_), ends(ends_)
    {
    }

    virtual ~SliceIntoOperation() {
    }

    virtual void execute() {
      source->execute();
      this->getResult()->getMachineTensor()->slice(
        F(1),
        source->getResult()->getMachineTensor(),
        // read from the entire result tensor of the source
        std::vector<size_t>(source->getResult()->getLens().size()),
        source->getResult()->getLens(),
        F(this->beta),
        begins, ends
      );
    }

    virtual operator std::string () const {
      return "SliceInto( " + std::string(*source) + ", " +
        SliceOperation<F,TE>::coordinateString(begins) + "-" +
        SliceOperation<F,TE>::coordinateString(ends) + " )";
    }

  protected:
    PTR(ESC(TensorOperation<F,TE>)) source;
    std::vector<size_t> begins, ends;

    static PTR(ESC(SliceIntoOperation<F,TE>)) create(
      const PTR(ESC(TensorOperation<F,TE>)) &source,
      const PTR(ESC(Tensor<F,TE>)) &result,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      return NEW(ESC(SliceIntoOperation<F,TE>),
        source, result, begins, ends, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

