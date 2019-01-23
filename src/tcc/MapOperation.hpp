/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MAP_OPERATION_DEFINED
#define TCC_MAP_OPERATION_DEFINED

#include <tcc/Operation.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename Target, typename Domain, typename TE> class Map;

  template <typename Target, typename Domain, typename TE>
  class MapOperation: public IndexedTensorOperation<Target,TE> {
  public:
    MapOperation(
      const std::function<Target(const Domain)> &f_,
      const PTR(ESC(IndexedTensorOperation<Domain,TE>)) &source_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<Target,TE>(
        Tensor<Target,TE>::create(
          source_->getResult()->getLens(),  // target tensor has identical lens
          "f(" + source_->getResult()->getName() + ")"
        ),
        source_->getResultIndices().c_str(), // target has identical indices
        source_->costs,
        typename Operation<TE>::ProtectedToken()
      ),
      f(f_), source(source_)
    {
      // TODO: assess map operation costs
    }
    virtual ~MapOperation() {
    }

    virtual void execute() {
      source->execute();
      // execute machine tensor's move with custom map
      this->getResult()->getMachineTensor()->move(
        Domain(1),
        source->getResult()->getMachineTensor(), source->getResultIndices(),
        Target(0),
        this->getResultIndices(),
        f
      );
    }

  protected:
    static PTR(ESC(MapOperation<Target,Domain,TE>)) create(
      const std::function<Target(const Domain)> &f_,
      const PTR(ESC(IndexedTensorOperation<Domain,TE>)) &source_
    ) {
      return NEW(ESC(MapOperation<Target,Domain,TE>),
        f_, source_, typename Operation<TE>::ProtectedToken()
      );
    }

    std::function<Target(const Domain)> f;
    PTR(ESC(IndexedTensorOperation<Domain,TE>)) source;

    friend class Map<Target,Domain,TE>;
  };
}

#endif

