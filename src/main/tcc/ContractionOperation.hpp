/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_OPERATION_DEFINED
#define TCC_CONTRACTION_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <tcc/Tensor.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class ContractionOperation: public IndexedTensorOperation<F,TE> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    ContractionOperation(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &left_,
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &right_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      Costs contractionCosts,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_, left_->costs + right_->costs,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      left(left_), right(right_)
    {
      // so far, costs contains costs involved to get left and right factors
      // during contraction all elements of left,right and result are present
      contractionCosts.maxElementsCount =
        contractionCosts.elementsCount + this->costs.elementsCount;
      // the intermediate results are, however, no longer needed afterwards
      this->costs.elementsCount = 0;
      this->costs += contractionCosts;
    }

    void execute() override {
      left->execute();
      right->execute();
      if (
        this->template isOlderThan<F>(left) ||
        this->template isOlderThan<F>(right)
      ) {
        LOG_FILE_LINE(2, this->file, this->line) << "executing: contraction " <<
          this->getName() << " <<= " <<
          this->alpha << " * " <<
          left->getName() << " * " << right->getName() << " + " <<
          this->beta << " * " << this->getName() << std::endl;

        this->getResult()->getMachineTensor()->contract(
          this->alpha,
          left->getResult()->getMachineTensor(), left->getResultIndices(),
          right->getResult()->getMachineTensor(), right->getResultIndices(),
          this->beta,
          this->getResultIndices()
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_FILE_LINE(3, this->file, this->line) <<
          this->getResult()->getName() << " up-to-date with " <<
          "(" << left->getName() << ", " << right->getName() << ")" <<
          std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return std::max(
        left->getLatestSourceVersion(),
        right->getLatestSourceVersion()
      );
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Contraction( " << this->alpha << ", " <<
        std::string(*left) << ", " << std::string(*right) << ", " <<
        this->beta << " )";
      return stream.str();
    }

  protected:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     **/
    static PTR(ESC(ContractionOperation<F,TE>)) create(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &left_,
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &right_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      const Costs &contractionCosts,
      const Scope &scope
    ) {
      return NEW(ESC(ContractionOperation<F,TE>),
        left_, right_,
        result_, resultIndices_,
        contractionCosts,
        scope.file, scope.line,
        typename Operation<TE>::ProtectedToken()
      );
    }

    PTR(ESC(IndexedTensorOperation<F,TE>)) left;
    PTR(ESC(IndexedTensorOperation<F,TE>)) right;

    friend class Contraction<F,TE>;
  };
}

#endif

