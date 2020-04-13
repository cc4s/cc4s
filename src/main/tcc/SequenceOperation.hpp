/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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
      const std::vector<PTR(Operation<TE>)> &operations_,
      const typename Operation<TE>::ProtectedToken &
    ): Operation<TE>(operations_[0]->costs), operations(operations_) {
      for (size_t i(1); i < operations_.size(); ++i) {
        this->costs += operations_[i]->costs;
      }
      // no elements are required by the result of a sequence
      this->costs.elementsCount = 0;
    }
    virtual ~SequenceOperation() {
    }

    virtual void execute() {
      // execute each operation in turn
      for (auto &operation: operations) {
        operation->execute();
      }
    }

    virtual operator std::string () const {
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
    static PTR(SequenceOperation<TE>) create(
      const std::vector<PTR(Operation<TE>)> &operations_
    ) {
      return NEW(SequenceOperation<TE>,
        operations_, typename Operation<TE>::ProtectedToken()
      );
    }

    std::vector<PTR(Operation<TE>)> operations;

    friend class Sequence<TE>;
  };
}

#endif

