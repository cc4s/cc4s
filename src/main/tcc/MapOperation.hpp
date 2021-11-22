/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TCC_MAP_OPERATION_DEFINED
#define TCC_MAP_OPERATION_DEFINED

#include <tcc/Operation.hpp>

#include <util/SharedPointer.hpp>
#include <functional>

namespace cc4s {
  template <typename Target, typename Domain, typename TE> class Map;

  template <typename Target, typename Domain, typename TE>
  class MapOperation: public IndexedTensorOperation<Target,TE> {
  public:
    MapOperation(
      const std::function<Target(Domain)> &f_,
      const Ptr<IndexedTensorOperation<Domain,TE>> &source_,
      const Costs &mapCosts_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<Target,TE>(
        Tensor<Target,TE>::create(
          source_->getResult()->getLens(),  // target tensor has identical lens
          "f(" + source_->getResult()->getName() + ")"
        ),
        source_->getResultIndices().c_str(), // target has identical indices
        mapCosts_, source_->costs,
        file_, line_,
        typename Operation<TE>::ProtectedToken()
      ),
      f(f_), source(source_)
    {
      // TODO: assess map operation costs
    }

    virtual ~MapOperation() {
    }

    void execute() override {
      source->execute();
      if (this->template isOlderThan<Domain>(source)) {
        LOG_LOCATION(SourceLocation(this->file, this->line)) <<
          "executing: unary map " <<
          this->getName() << " <<= "<<
          "f(" << this->alpha << " * " << source->getName() << ") + " <<
          this->beta << " * " << this->getName() << std::endl;

        // execute machine tensor's sum with custom map
        this->getResult()->getMachineTensor()->sum(
          Domain(1),
          source->getResult()->getMachineTensor(), source->getResultIndices(),
          Target(0), this->getResultIndices(),
          f
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
      return "map( f, " + std::string(*source) + " )";
    }

 protected:
    static Ptr<MapOperation<Target,Domain,TE>> create(
      const std::function<Target(Domain)> &f_,
      const Ptr<IndexedTensorOperation<Domain,TE>> &source_,
      const Scope &scope
    ) {
      auto elementsCount(source_->getResult()->getElementsCount());
      return New<MapOperation<Target,Domain,TE>>(
        f_, source_,
        // FIXME: costs of map assumed 10 times costs of addition
        // but depends on actual map
        Costs(
          elementsCount, elementsCount,
          0, 10*elementsCount
        ),
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    std::function<Target(Domain)> f;
    Ptr<IndexedTensorOperation<Domain,TE>> source;

    friend class Map<Target,Domain,TE>;
  };
}

#endif

