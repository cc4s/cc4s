#ifndef TCC_OPERATION_SEQUENCE_DEFINED
#define TCC_OPERATION_SEQUENCE_DEFINED

#include <tcc/Operation.hpp>

#include <util/SharedPointer.hpp>
#include <vector>

namespace cc4s {
  template <typename TE> class Sequence;

  template <typename TE>
  class SequenceOperation: public Operation<TE> {
  public:
    SequenceOperation(
      const std::vector<Ptr<Operation<TE>>> &operations_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ): Operation<TE>(
      operations_[0]->costs, file_, line_
    ), operations(
      operations_
    ) {
      for (size_t i(1); i < operations_.size(); ++i) {
        this->costs += operations_[i]->costs;
      }
      // no elements are required by the result of a sequence
      this->costs.elementsCount = 0;
    }

    void execute() override {
      // execute each operation in turn
      for (auto &operation: operations) {
        operation->execute();
      }
    }

    size_t getLatestSourceVersion() override {
      size_t latestVersion(0);
      for (auto &operation: operations) {
        latestVersion = std::max(
          latestVersion,
          operation->getLatestSourceVersion()
        );
      }
      return latestVersion;
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Sequence( ";
      std::string delimiter("");
      for (auto const &operation: operations) {
        stream << delimiter << std::string(*operation);
        delimiter = ", ";
      }
      stream << " )";
      return stream.str();
    }

  protected:
    static Ptr<SequenceOperation<TE>> create(
      const std::vector<Ptr<Operation<TE>>> &operations_,
      const Scope &scope
    ) {
      return New<SequenceOperation<TE>>(
        operations_,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    std::vector<Ptr<Operation<TE>>> operations;

    friend class Sequence<TE>;
  };
}

#endif

