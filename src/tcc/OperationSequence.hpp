/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_OPERATION_SEQUENCE_DEFINED
#define TCC_OPERATION_SEQUENCE_DEFINED

#include <tcc/Operation.hpp>

#include <vector>
#include <memory>

namespace tcc {
  template <typename F>
  class OperationSequence: public Operation<F> {
  public:
    OperationSequence(
      const std::vector<std::shared_ptr<Operation<F>>> &operations_,
      const typename Operation<F>::ProtectedToken &
    ): Operation<F>(operations_[0]->costs), operations(operations_) {
      for (unsigned int i(1); i < operations_.size(); ++i) {
        this->costs += operations_[i]->costs;
      }
      // no elements are required by the result of a sequence
      this->costs.elementsCount = 0;
    }
    virtual ~OperationSequence() {
    }

    virtual void execute() {
      // execute each operation in turn
      for (unsigned int i(0); i < operations.size(); ++i) {
        operations[i]->execute();
      }
    };

    // this operation returns void
    virtual std::shared_ptr<Tensor<F>> getResult() {
      return std::shared_ptr<Tensor<F>>();
    }
    virtual const std::string &getResultIndices() {
      // TODO: think of better way
      return emptyIndices;
    }

  protected:
    static std::shared_ptr<OperationSequence<F>> create(
      const std::vector<std::shared_ptr<Operation<F>>> &operations_
    ) {
      return std::make_shared<OperationSequence<F>>(
        operations_, typename Operation<F>::ProtectedToken()
      );
    }

    std::vector<std::shared_ptr<Operation<F>>> operations;
    std::string emptyIndices;

    friend class Tcc<F>;
  };
}

#endif
