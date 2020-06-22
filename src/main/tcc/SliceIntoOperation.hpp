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
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, source_->costs,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      source(source_), begins(begins_), ends(ends_)
    {
    }

    void execute() override {
      source->execute();
      if (this->template isOlderThan<F>(source)) {
        // read from the entire result tensor of the source
        auto aEnds(source->getResult()->getLens());
        auto aBegins(std::vector<size_t>(aEnds.size()));

        std::stringstream beginsStream, endsStream, aBeginsStream, aEndsStream;
        for (auto d: begins) { beginsStream << " " << d; }
        for (auto d: ends) { endsStream << " " << d; }
        for (auto d: aBegins) { aBeginsStream << " " << d; }
        for (auto d: aEnds) { aEndsStream << " " << d; }

        LOG_FILE_LINE(2, this->file, this->line) << "executing: slice " <<
          this->getName() << "(" <<
            beginsStream.str() << "," << endsStream.str() <<
          ") <<= " << this->alpha << " * " << source->getName() << "(" <<
            aBeginsStream.str() << "," << aEndsStream.str() << ") + " <<
          this->beta << " * " << this->getName() << "(" <<
            beginsStream.str() << "," << endsStream.str() <<
          ")" << std::endl;

        this->getResult()->getMachineTensor()->slice(
          F(1), source->getResult()->getMachineTensor(), aBegins, aEnds,
          F(this->beta), begins, ends
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_FILE_LINE(3,this->file, this->line) << this->getName() <<
          " up-to-date with " << source->getName() << std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return source->getLatestSourceVersion();
    }

    operator std::string () const override {
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
      const std::vector<size_t> ends,
      const Scope &scope
    ) {
      return NEW(ESC(SliceIntoOperation<F,TE>),
        source, result, begins, ends,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

