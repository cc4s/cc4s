/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MAP_OPERATION_DEFINED
#define TCC_MAP_OPERATION_DEFINED

#include <tcc/Operation.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename F> class Map;

  template <typename F>
  class MapOperation: public Operation<F> {
  public:
    MapOperation(
      const std::function<F(const F)> &f_,
      const PTR(Operation<F>) &source_,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(source_->costs), f(f_), source(source_),
      result(
        source_->getResult()->getTcc()->createTensor(
          source_->getResult(), "f(" + source->getResult()->getName() + ")"
        )
      )
    {
      // TODO: assess map operation costs
    }
    virtual ~MapOperation() {
    }

    virtual void execute() {
      source->execute();
      // execute machine tensor's move with custom map
      result->getMachineTensor()->move(
        F(1),
        source->getResult()->getMachineTensor(), source->getResultIndices(),
        F(0),
        getResultIndices(),
        f
      );
    }

    virtual PTR(Tensor<F>) getResult() {
      return result;
    }
    virtual const std::string &getResultIndices() {
      return source->getResultIndices();
    }

  protected:
    static PTR(MapOperation<F>) create(
      const std::function<F(const F)> &f_,
      const PTR(Operation<F>) &source_
    ) {
      return NEW(MapOperation<F>,
        f_, source_, typename Operation<F>::ProtectedToken()
      );
    }

    std::function<F(const F)> f;
    PTR(Operation<F>) source;
    PTR(Tensor<F>) result;

    friend class Map<F>;
  };
}

#endif

