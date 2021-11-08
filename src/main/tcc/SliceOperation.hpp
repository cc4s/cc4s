#ifndef TCC_SLICE_OPERATION_DEFINED
#define TCC_SLICE_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>
#include <vector>

namespace cc4s {
  template <typename F, typename TE> class Slice;

  template <typename F, typename TE>
  class SliceOperation: public TensorOperation<F,TE> {
  public:
    SliceOperation(
      const Ptr<TensorOperation<F,TE>> &source_,
      const Ptr<Tensor<F,TE>> &result_,
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
        // write into entire result tensor
        auto bEnds(this->getResult()->getLens());
        auto bBegins(std::vector<size_t>(bEnds.size()));

        std::stringstream beginsStream, endsStream, bBeginsStream, bEndsStream;
        for (auto d: begins) { beginsStream << " " << d; }
        for (auto d: ends) { endsStream << " " << d; }
        for (auto d: bBegins) { bBeginsStream << " " << d; }
        for (auto d: bEnds) { bEndsStream << " " << d; }

        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: slice " <<
          this->getName() << "(" <<
            bBeginsStream.str() << "," << bEndsStream.str() <<
          ") <<= " << this->alpha << " * " << source->getName() << "(" <<
            beginsStream.str() << "," << endsStream.str() << ") + " <<
          this->beta << " * " << this->getName() << "(" <<
            bBeginsStream.str() << "," << bEndsStream.str() <<
          ")" << std::endl;

        this->getResult()->getMachineTensor()->slice(
          F(1), source->getResult()->getMachineTensor(), begins, ends,
          F(0), bBegins, bEnds
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() <<
          " up-to-date with " << source->getName() << std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return source->getLatestSourceVersion();
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
    Ptr<TensorOperation<F,TE>> source;
    std::vector<size_t> begins, ends;

    static Ptr<SliceOperation<F,TE>> create(
      const Ptr<TensorOperation<F,TE>> &source,
      const Ptr<Tensor<F,TE>> &result,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends,
      const Scope &scope
    ) {
      return New<SliceOperation<F,TE>>(
        source, result, begins, ends,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Slice<F,TE>;
  };
}

#endif

