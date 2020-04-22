/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SLICE_OPERATION_DEFINED
#define TCC_SLICE_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Slice;

  template <typename F, typename TE>
  class SliceOperation: public TensorOperation<F,TE> {
  public:
    SliceOperation(
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

    void execute() override {
      source->execute();
      if (this->template isOlderThan<F>(source)) {
        this->getResult()->getMachineTensor()->slice(
          F(1),
          source->getResult()->getMachineTensor(),
          begins, ends,
          F(0),
          // write into entire result tensor
          std::vector<size_t>(this->getResult()->getLens().size()),
          this->getResult()->getLens()
        );
        this->updated();
      }
    }

    operator std::string () const override {
      return "Slice( " + std::string(*source) + ", " +
        SliceOperation<F,TE>::coordinateString(begins) + "-" +
        SliceOperation<F,TE>::coordinateString(ends) + " )";
    }

    static std::string coordinateString(std::vector<size_t> const &coordinates) {
      std::stringstream stream;
      stream << "(";
      std::string delimiter("");
      for (auto i: coordinates) {
        stream << delimiter << i;
        delimiter = ",";
      }
      stream << ")";
      return stream.str();
    }

  protected:
    PTR(ESC(TensorOperation<F,TE>)) source;
    std::vector<size_t> begins, ends;

    static PTR(ESC(SliceOperation<F,TE>)) create(
      const PTR(ESC(TensorOperation<F,TE>)) &source,
      const PTR(ESC(Tensor<F,TE>)) &result,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      return NEW(ESC(SliceOperation<F,TE>),
        source, result, begins, ends, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

