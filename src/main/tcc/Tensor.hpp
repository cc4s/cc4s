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

#ifndef TCC_TENSOR_DEFINED
#define TCC_TENSOR_DEFINED

#include <tcc/TensorExpression.hpp>

#include <tcc/TensorLoadOperation.hpp>
#include <util/SharedPointer.hpp>
#include <math/Integer.hpp>
#include <Cc4s.hpp>
#include <cstdint>
#include <vector>
#include <string>
#include <mpi.h>

namespace cc4s {
  size_t getNextTensorVersion();

  class TensorDimensionProperty {
  public:
    // property name for reference in tensor expressions,
    // e.g. "spin" for the notation "i.spin"
    std::string name;
    // map: index -> property
    std::map<Natural<>,Natural<>> propertyOfIndex;
    // index for reverse-lookup of the map index -> property
    std::map<Natural<>,std::set<Natural<>>> indicesOfProperty;
  };

  class TensorDimension {
  public:
    // dimension name for reference in tensor expressions,
    // e.g. "State" for indicating that a dimension referes to a state p
    std::string name;
    std::map<std::string, Ptr<TensorDimensionProperty>> properties;

    // table of tensor dimension information
    static std::map<std::string, Ptr<TensorDimension>> dimensions;
  };

  class TensorDimensionPropertyReference {
  public:
    Natural<> dimension;
    Ptr<TensorDimensionProperty> property;
  };

  class TensorNonZeroCondition {
  public:
    // name of condition for reference
    std::string name;
    std::vector<TensorDimensionPropertyReference> dimensionPropertyReferences;
    // list tuples of dimension properties
    std::vector<std::vector<Natural<>>> tuples;
  };

  class TensorNonZeroConditions {
  public:
    // NOTE: currently, only "all" is supported, i.e. all of the listed
    // conditions must be met for a tensor entry to be non-zero.
    std::vector<Ptr<TensorNonZeroCondition>> all;
  };

  class TensorNonZeroBlockElementIterator {
  public:
    std::vector<Natural<>> lens;
    std::vector<std::vector<Natural<>>> indices;
    std::vector<Natural<>> indexPositions;

    TensorNonZeroBlockElementIterator(
      const std::vector<Natural<>> &lens_,
      std::vector<std::vector<Natural<>>> &&indices_
    ): lens(lens_), indices(indices_), indexPositions(indices_.size(), 0) {
    }

    TensorNonZeroBlockElementIterator &operator ++() {
      if (!atEnd()) {
        ++(indexPositions[0]);
        Natural<> d(0);
        while (
          d+1 < indexPositions.size()  &&
          indexPositions[d] >= indices[d].size()
        ) {
          // carry?
          indexPositions[d] = 0;
          ++d;
          ++(indexPositions[d]);
        }
      }
      return *this;
    }

    bool atEnd() const {
      if (indexPositions.size() == 0) return true;
      const Natural<> lastIndex(indexPositions.size()-1);
      return indexPositions[lastIndex] >= indices[lastIndex].size();
    }

    Natural<> getGlobalIndex() const {
      Natural<> I(0), D(indexPositions.size()-1);
      for (Natural<> d(0); d < indexPositions.size(); ++d) {
        I *= lens[D-d];
        I += indices[D-d][indexPositions[D-d]];
      }
      return I;
    }

    Natural<> getCount() const {
      Natural<> count(1);
      for (Natural<> d(0); d < indices.size(); ++d) {
        count *= indices[d].size();
      }
      return count;
    }

    std::string print() {
      std::stringstream s;
      s << "(";
      for (Natural<> d(0); d < indexPositions.size(); ++d) {
        s << " " << indices[d][indexPositions[d]] << " ";
      }
      s << ")";
      return s.str();
    }
  };

  class TensorNonZeroBlockIterator {
  public:
    // the non-zero conditions over which this iterator iterates
    Ptr<TensorNonZeroConditions> nonZeroConditions;
    // references the position within the tuple list in each non-zero condition
    std::vector<Natural<>> tuplePositions;
    bool ended;

    TensorNonZeroBlockIterator(
      const Ptr<TensorNonZeroConditions> &nonZeroConditions_
    ):
      nonZeroConditions(nonZeroConditions_),
      tuplePositions(nonZeroConditions_->all.size(), 0),
      ended(false)
    {
    }

    TensorNonZeroBlockIterator(
      const Ptr<TensorNonZeroConditions> &nonZeroConditions_,
      std::vector<Natural<>> &&tuplePositions_
    ):
      nonZeroConditions(nonZeroConditions_),
      tuplePositions(tuplePositions_),
      ended(false)
    {
    }

    TensorNonZeroBlockIterator &operator ++() {
      if (ended) return *this;
      if (tuplePositions.size() > 0) {
        ++(tuplePositions[0]);
        Natural<> d(0);
        while (
          d+1 < tuplePositions.size()  &&
          tuplePositions[d] >= nonZeroConditions->all[d]->tuples.size()
        ) {
          // carry?
          tuplePositions[d] = 0;
          ++d;
          ++(tuplePositions[d]);
        }
        if (tuplePositions[d] >= nonZeroConditions->all[d]->tuples.size()) {
          ended = true;
        }
      } else {
        ended = true;
      }
      return *this;
    }

    bool atEnd() const {
      return ended;
    }

    TensorNonZeroBlockElementIterator getElementIterator(
      const std::vector<Natural<>> &lens
    ) {
      // build index map
      std::vector<Natural<>> conditionsCountOnDimension(lens.size());
      std::vector<std::vector<Natural<>>> matchCountOnDimensionIndex(lens.size());
      for (Natural<> d(0); d < lens.size(); ++d) {
        matchCountOnDimensionIndex[d].resize(lens[d]);
      }
      // go through all conditions
      Natural<> c(0);
      for (auto tp: tuplePositions) {
        auto nonZeroCondition(nonZeroConditions->all[c]);
        for (Natural<> dp(0); dp < nonZeroCondition->dimensionPropertyReferences.size(); ++dp) {
          auto dimensionPropertyReference(
            nonZeroCondition->dimensionPropertyReferences[dp]
          );
          auto d(dimensionPropertyReference.dimension);
          auto property(dimensionPropertyReference.property);
          // one more condition on dimension d
          ++conditionsCountOnDimension[d];
          // go through all selected indices of this dimension and mark them
          auto selectedProperty(nonZeroCondition->tuples[tp][dp]);
          auto indices(property->indicesOfProperty[selectedProperty]);
          for (auto index: indices) {
            ++matchCountOnDimensionIndex[d][index];
          }
        }
        ++c;
      }
      // finally, build list of indices, mathcing all conditions
      std::vector<std::vector<Natural<>>> indices(lens.size());
      for (Natural<> d(0); d < lens.size(); ++d) {
        for (Natural<> i(0); i < lens[d]; ++i) {
          if (matchCountOnDimensionIndex[d][i] == conditionsCountOnDimension[d]) {
            indices[d].push_back(i);
          }
        }
      }
      return TensorNonZeroBlockElementIterator(lens, std::move(indices));
    }

    std::string print() const {
      std::stringstream s;
      Natural<> d(0);
      for (auto tp: tuplePositions) {
        auto nonZeroCondition(nonZeroConditions->all[d]);
        for (Natural<> i(0); i < nonZeroCondition->dimensionPropertyReferences.size(); ++i) {
          auto dimensionPropertyReference(
            nonZeroCondition->dimensionPropertyReferences[i]
          );
          s << " " <<
            dimensionPropertyReference.dimension << "." <<
            dimensionPropertyReference.property->name << ": " <<
            nonZeroCondition->tuples[tp][i];
        }
        s << ";";
        ++d;
      }
      return s.str();
    }
  };


  /**
   * \brief
   **/
  template <typename F, typename TE>
  class Tensor: public TensorExpression<F,TE> {
  public:
  protected:
    /**
     * \brief Dummy objects of that type are used to guarantee that
     * constructors are only called from within the class although they
     * need to be declared public to work with make_shared aka New.
     **/
    class ProtectedToken {
    };

  public:
    // TODO: should be protected
    /**
     * \brief Dimensions of the tensor if assumedShape is true, otherwise empty.
     **/
    std::vector<size_t> lens;

    /**
     * \brief Whether this tensor has already well defined dimensions.
     * A tensor may be created without shape. Such a tensor assumes shape
     * from the result of an assignment.
     **/
    bool assumedShape;

    std::vector<Ptr<TensorDimension>> dimensions;
    Ptr<TensorNonZeroConditions> nonZeroConditions;

  protected:
    /**
     * \brief The version of the data in this tensor. It is updated everytime
     * the tensor is updated.
     **/
    size_t version;

    /**
     * \brief The tensor name for in verbose info and error messages.
     **/
    std::string name;

    /**
     * \brief Pointer to the machine tensor, wrapping the underlying
     * tensor representation. The machine tensor is not allocated until
     * requested during execution or by manually calling getMachineTensor().
     **/
    typedef typename TE::template MachineTensor<F> MT;
    Ptr<MT> machineTensor;

  public:
    /**
     * \brief Create a tcc tensor of yet unknown shape.
     * It must be first on the left-hand-side of an assignment.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::string &name_,
      const ProtectedToken &
    ): assumedShape(false), version(0), name(name_) {
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const std::vector<size_t> &lens_,
      const std::string &name_,
      const ProtectedToken &
    ):
      lens(lens_), assumedShape(true), dimensions(lens_.size()),
      version(0), name(name_)
    {
      // the machine tensor is not allocated initially
    }

    Tensor(
      const std::vector<size_t> &lens_,
      const std::string &name_,
      const bool assumedShape_,
      const ProtectedToken &
    ):
      lens(lens_), assumedShape(assumedShape_), dimensions(lens_.size()),
      version(0), name(name_)
    {
      // the machine tensor is not allocated initially
    }

    /**
     * \brief Create a tensor from a given machine tensor.
     * Not intended for direct invocation. Use Tcc::createTensor instead.
     **/
    Tensor(
      const typename MT::T &unadaptedTensor_,
      const ProtectedToken &
    ): assumedShape(true), version(0) {
      auto mt(MT::create(unadaptedTensor_));
      lens = mt->getLens();
      name = mt->getName();
      dimensions.resize(lens.size());
      machineTensor = mt;
    }

    static Ptr<Tensor<F,TE>> create(const std::string &name) {
      return New<Tensor<F,TE>>(name, ProtectedToken());
    }

    static Ptr<Tensor<F,TE>> create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return New<Tensor<F,TE>>(lens, name, ProtectedToken());
    }

    /**
     * \brief Create an empty tensor of identical shape as the given tensor.
     * The name, however, should differ.
     **/
    static Ptr<Tensor<F,TE>> create(
      const Ptr<Tensor<F,TE>> &tensor,
      const std::string &name
    ) {
      return New<Tensor<F,TE>>(
        tensor->lens, name, tensor->assumedShape, ProtectedToken()
      );
    }

    static Ptr<Tensor<F,TE>> create(
      const typename MT::T &unadaptedTensor
    ) {
      return New<Tensor<F,TE>>(unadaptedTensor, ProtectedToken());
    }

    void setName(const std::string &name_) {
      name = name_;
    }
    const std::string &getName() const {
      return name;
    }

    Ptr<MT> getMachineTensor() {
      if (!machineTensor) {
        ASSERT_LOCATION(assumedShape,
          "Tried to execute operation on tensor " + name +
          " before its shape has been assumed.",
          SOURCE_LOCATION
        );
        // allocate the implementation specific machine tensor upon request
        machineTensor = MT::create(lens, name);
        // wait until allocation is done on all processes
        Cc4s::world->barrier();
        LOG() << "Allocated tensor with " <<
          getElementsCount() << " elements" << std::endl;
      }
      return machineTensor;
    }

    bool allocated() {
      return machineTensor != nullptr;
    }

    size_t getVersion() const {
      return version;
    }

    void updated() {
      version = getNextTensorVersion();
    }

    // read tensor elements to buffer
    // each rank must specify its range of indices to read
    // behavior is undefined if ranges overlap
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
      getMachineTensor()->read(elementsCount, indexData, valueData);
    }
    F read(const size_t index = 0) {
      std::vector<size_t> indices(1);
      std::vector<F> values(1);
      if (Cc4s::world->getRank() == 0) {
        indices[0] = index;
        read(1, indices.data(), values.data());
      } else {
        read(0, indices.data(), values.data());
      }
      // broadcast value read on root
      Cc4s::world->broadcast(values);
      return values[0];
    }
    std::vector<F> readAll() {
      std::vector<F> values(getElementsCount());
      if (Cc4s::world->getRank() == 0) {
        std::vector<size_t> indices(values.size());
        for (size_t i(0); i < values.size(); ++i) {
          indices[i] = i;
        }
        read(indices.size(), indices.data(), values.data());
      } else {
        std::vector<size_t> indices(0);
        read(0, indices.data(), values.data());
      }
      // broadcast value read on root
      Cc4s::world->broadcast(values);
      return values;
    }
    void readToFile(MPI_File &file, const size_t offset = 0) {
      getMachineTensor()->readToFile(file, offset);
    }

    // write tensor elements from buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
      getMachineTensor()->write(elementsCount, indexData, valueData);
      updated();
    }
    void write(F value, size_t index = 0) {
      write(1, &index, &value);
    }
    void writeFromFile(MPI_File &file, const size_t offset = 0) {
      getMachineTensor()->writeFromFile(file, offset);
    }

    const std::vector<Natural<>> &getLens() const {
      return lens;
    }
    Natural<> getLen(Natural<> d) const {
      return lens[d];
    }


    /**
     * \brief Returns the number of elements contained in this tensor.
     **/
    size_t getElementsCount() const {
      size_t elementsCount(1);
      for (auto len: lens) {
        elementsCount *= len;
      }
      return elementsCount;
    }

    Ptr<Tensor<F,TE>> inspect() override {
      return this->template toPtr<Tensor<F,TE>>();
    }

    Ptr<Tensor<F,TE>> evaluate() override {
      return this->template toPtr<Tensor<F,TE>>();
    }

    Ptr<Operation<TE>> compile(Scope &scope) override {
      return TensorLoadOperation<F,TE>::create(
        this->template toPtr<Tensor<F,TE>>(),
        Costs(getElementsCount()),
        scope
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    Ptr<TensorOperation<F,TE>> lhsCompile(
      const Ptr<TensorOperation<F,TE>> &rhsOperation
    ) override {
      // shape assuming:
      if (!assumedShape) {
        // lhs has not yet assumed shape
        if (!rhsOperation->getResult()->assumedShape) {
          // rhs has not yet assumed shape either
          ASSERT_LOCATION(false,
            "Neither left-hand-side tensor " + getName() +
            " nor right-hand-side result " +
            rhsOperation->getResult()->getName() +
            " have known shape. "
            "Try specifying left-hand-side tensor dimension manually.",
            SourceLocation(rhsOperation->file, rhsOperation->line)
          );
        } else {
          // let this lhs tensor assume the shape of the rhs result
          lens = rhsOperation->getResult()->getLens();
          assumedShape = true;
          // also assume name if empy
          if (name == "") name = rhsOperation->getResult()->getName();
        }
      } else {
        if (!rhsOperation->getResult()->assumedShape) {
          // let rhs tensor assume the shape of this lhs tensor
          rhsOperation->getResult()->lens = getLens();
          rhsOperation->getResult()->assumedShape = true;
        } else {
          // both have assumed shape: shapes must match
          if (rhsOperation->getResult()->getLens() != lens) {
            std::stringstream lhsShape;
            for (auto i: getLens()) { lhsShape << " " << i; }
            std::stringstream rhsShape;
            for (auto i: rhsOperation->getResult()->getLens()) { rhsShape << " " << i; }
            ASSERT_LOCATION(false,
              "Shape of left-hand-side tensor " + getName() +
              " (" + lhsShape.str() + ") "
              " must match the shape of the result tensor " +
              rhsOperation->getResult()->getName() +
              " (" + rhsShape.str() + ")",
              SourceLocation(rhsOperation->file, rhsOperation->line)
            );
          }
        }
      }

      // make the rhs operation directly operate on this tensor
      rhsOperation->result = this->template toPtr<Tensor<F,TE>>();
      return rhsOperation;
    }

    operator std::string () const override {
      return getName();
    }
  };
}

#endif

