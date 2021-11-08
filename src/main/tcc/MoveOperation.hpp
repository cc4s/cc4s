#ifndef TCC_MOVE_OPERATION_DEFINED
#define TCC_MOVE_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class MoveOperation: public IndexedTensorOperation<F,TE> {
  public:
    /**
     * \brief Creates a move operation moving the results of
     * the right hand side operation into given tensor of the
     * left hand side after applying the function f.
     * The function f defaults to the identity operation.
     * Not intended for direct invocation. Use expression->compile()
     * to generate operations.
     **/
    MoveOperation(
      const Ptr<IndexedTensorOperation<F,TE>> &rhs_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      Costs moveCosts,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_,
        rhs_->costs,
        file_, line_,
        typename Operation<TE>::ProtectedToken()
      ),
      rhs(rhs_)
    {
      this->costs += moveCosts;
    }

    void execute() override {
      rhs->execute();
      if (this->template isOlderThan<F>(rhs)) {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: sum " <<
          this->getName() << " <<= " <<
          this->alpha << " * " << rhs->getName() << " + " <<
          this->beta << " * " << this->getName() << std::endl;

        this->getResult()->getMachineTensor()->sum(
          this->alpha,
          rhs->getResult()->getMachineTensor(), rhs->getResultIndices(),
          this->beta,
          this->resultIndices
        );
        this->updated();
        this->accountFlops();
      } else {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          this->getName() << " up-to-date with " << rhs->getName() << std::endl;
      }
    }

    size_t getLatestSourceVersion() override {
      return rhs->getLatestSourceVersion();
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Move( " << this->alpha << ", " <<
        std::string(*this->result) << ", " << std::string(*rhs) << ", " <<
        this->beta << " )";
      return stream.str();
    }

  protected:
    static Ptr<MoveOperation<F,TE>> create(
      const Ptr<IndexedTensorOperation<F,TE>> &rhs_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &moveCosts,
      const Scope &scope
    ) {
      return New<MoveOperation<F,TE>>(
        rhs_,
        result_, resultIndices_, moveCosts,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    Ptr<IndexedTensorOperation<F,TE>> rhs;
//    F alpha, beta;

    friend class Contraction<F,TE>;
    friend class Indexing<F,TE>;
  };
}

#endif

