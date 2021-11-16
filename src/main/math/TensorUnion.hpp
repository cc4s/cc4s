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

#ifndef TENSOR_UNION_DEFINED
#define TENSOR_UNION_DEFINED

#include <math/MathFunctions.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <ostream>

namespace cc4s {
  /**
   * \brief Represents the union of tensors and provides
   * vector space operations: addition, scalar multiplication, inner product,
   * and complex conjugation.
   **/
  template <typename F, typename TE>
  class TensorUnion: public Object {
  public:
    typedef F FieldType;

    std::vector<Ptr<TensorExpression<F,TE>>> componentTensors;
    std::vector<std::string> componentIndices;

    /**
     * \brief Default constructor for an empty tensor union without elements.
     **/
    TensorUnion() {
    }

    /**
     * \brief Move constructor taking possession of the tensors owned by a.
     **/
    TensorUnion(
      TensorUnion &&a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices),
      indexEnds(a.componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Copy constructor copying the tensors owned by a.
     **/
    TensorUnion(
      const TensorUnion &a
    ):
      componentTensors(a.componentTensors.size()),
      componentIndices(a.componentIndices),
      indexEnds(a.componentTensors.size())
    {
      copyComponents(a.componentTensors);
      buildIndexTranslation();
    }

    /**
     * \brief Move constructor taking possession of the tensors given.
     **/
    TensorUnion(
      const std::vector<Ptr<TensorExpression<F,TE>>> &tensors,
      const std::vector<std::string> &indices
    ):
      componentTensors(tensors),
      componentIndices(indices),
      indexEnds(componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Move constructor taking possession of the tensors given
     * by the iterators.
     **/
    template <typename TensorsIterator, typename IndicesIterator>
    TensorUnion(
      TensorsIterator tensorsBegin, TensorsIterator tensorsEnd,
      IndicesIterator indicesBegin, IndicesIterator indicesEnd
    ):
      componentTensors(tensorsBegin, tensorsEnd),
      componentIndices(indicesBegin, indicesEnd),
      indexEnds(componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Retrieves the i-th component tensor. Note that
     * the Tensor is not const since rearrangement may be
     * required also in non-modifying tensor operations.
     **/
    const Ptr<TensorExpression<F,TE>> &get(const size_t i) const {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component tensor.
     **/
    Ptr<TensorExpression<F,TE>> &get(const size_t i) {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component indices.
     **/
    const std::string &getIndices(const size_t i) const {
      return componentIndices[i];
    }

    /**
     * \brief Retrieves the i-th component indices as modifiable string.
     **/
    std::string &getIndices(const size_t i) {
      return componentIndices[i];
    }

    /**
     * \brief Move assignment operator taking possession of the tensors
     * owned by a.
     **/
    TensorUnion &operator =(const TensorUnion &&a) {
      componentTensors = a.componentTensors;
      componentIndices = a.componentIndices;
      buildIndexTranslation();
      return *this;
    }

    /**
     * \brief Copy assignment operator copying the tensors owned by a.
     **/
    TensorUnion &operator =(const TensorUnion &a) {
      componentIndices = a.componentIndices;
      copyComponents(a.componentTensors);
      buildIndexTranslation();
      return *this;
    }

    /**
     * \brief Add-to assignment operator adding each component of a
     * to the respective component of this TensorUnion.
     **/
    TensorUnion &operator += (const TensorUnion &a) {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size(); ++i) {
//        const char *indices(componentIndices[i].c_str());
//        get(i)->sum(+1.0, *a.get(i), indices, 1.0, indices);
        COMPILE(
          (*get(i))[getIndices(i)] += (*a.get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Subtract-from assignment operator subtracting each component of a
     * from the respective component of this TensorUnion.
     **/
    TensorUnion &operator -= (const TensorUnion &a) {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size(); ++i) {
//        const char *indices(getIndices(i).c_str());
//        get(i)->sum(-1.0, *a.get(i), indices, 1.0, indices);
        COMPILE(
          (*get(i))[getIndices(i)] -= (*a.get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Multiply-by assignment operator scalar multiplying each component
     * each component of this TensorUnion by the given scalar.
     **/
    TensorUnion &operator *= (const F s) {
      for (size_t i(0); i < componentTensors.size(); ++i) {
//        const char *indices(getIndices(i).c_str());
//        get(i)->sum(s, *get(i), indices, 0.0, indices);
        COMPILE(
          (*get(i))[getIndices(i)] <<= s * (*get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Creates and returns the conjugate transpose of this TensorUnion.
     * The first and the second half of the inidices in each component are
     * swapped for the transposition. For real types F the conjugation
     * does nothing.
     **/
    // TOOD: precompile operations
    TensorUnion conjugateTranspose() const {
      TensorUnion result;
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t order(getIndices(i).length() / 2);
        std::vector<int> transposedLens(get(i)->lens, get(i)->lens + 2*order);
        std::rotate(
          transposedLens.begin(),
          transposedLens.begin() + order,
          transposedLens.begin() + 2*order
        );
        result.componentTensors.push_back(
          Tcc<TE>::template tensor<F>(
            transposedLens,
            std::string(get(i)->get_name()) + "*"
          )
        );
        result.componentIndices.push_back(
          getIndices(i).substr(order, 2*order) + getIndices(i).substr(0, order)
        );
        COMPILE(
          (*result.get(i))[result.getIndices(i)] <<= map<F>(
            cc4s::conj<F>,
            (*get(i))[getIndices(i)]
          )
        )->execute();
      }
      return std::move(result);
    }

    /**
     * \brief Returns the matrix product of this bra-TensorUnion with the
     * given dual ket-TensorUnion ket.
     **/
    F braket(const TensorUnion &ket) const {
      checkDualCompatibility(ket);
      auto result( Tcc<TE>::template tensor<F>("dot") );
      for (size_t i(0); i < componentTensors.size(); ++i) {
        // add to result
        COMPILE(
          (*result)[""] +=
            (*get(i))[getIndices(i)] * (*ket.get(i))[ket.getIndices(i)]
        )->execute();
      }
      return result->read();
    }

    /**
     * \brief Returns the inner product of this ket-TensorUnion with the
     * given ket-TensorUnion a. The elements of this TensorUnion are conjugated
     * in the inner product, i.e. this->dot(a) yields the same results as
     * this->conjugateTranspose().braket(a).
     **/
    F dot(const TensorUnion &a) const {
      checkCompatibilityTo(a);
      auto result( Tcc<TE>::template tensor<F>("dot") );
      for (size_t i(0); i < componentTensors.size(); ++i) {
        // add to result
        COMPILE(
          (*result)[""] += (*get(i))[getIndices(i)] *
            map<F>(cc4s::conj<F>, (*a.get(i))[getIndices(i)])
        )->execute();
      }
      return result->read();
    }

    /**
     * \brief Get the number of component tensors of this TensorUnion.
     */
    size_t getComponentsCount() const {
      return componentTensors.size();
    }

    /**
     * \brief Get the total number of degrees of freedom represented by this
     * TensorUnion, i.e. the total number of field values contained in all
     * component tensors. The indices used by read and write are between
     * 0 and getDimension()-1.
     */
    size_t getDimension() const {
      return indexEnds.back();
    }

    /**
     * \Brief Translates the given component and component index into
     * its element into an index between 0 and getDimension()-1.
     **/
    size_t getIndex(const size_t component, const size_t componentIndex) const {
      size_t base(component > 0 ? indexEnds[component-1] : 0);
      return base + componentIndex;
    }

    /**
     * \Brief Translates the given index between 0 and getDimension()-1
     * into a component number and component index into the corresponding
     * component tensor.
     **/
    void fromIndex(
      const size_t index, size_t &component, size_t &componentIndex
    ) const {
      component = 0;
      size_t base(0);
      while (component < indexEnds.size()) {
        if (index < indexEnds[component]) break;
        base = indexEnds[component++];
      }
      ASSERT(component < indexEnds.size(), "Index out bounds");
      componentIndex = index - base;
    }

    /**
     * \brief Reads out all locally stored values together with their
     * respective indices. The indices are between 0 and getDimension()-1.
     **/
    std::vector<std::pair<size_t,F>> readLocal() const {
      size_t elementsCount(0);
      std::vector<std::pair<size_t,F>> elements;
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t componentValuesCount;
        size_t *componentIndices;
        F *componentValues;
        get(i)->read_local(
          reinterpret_cast<int64_t *>(&componentValuesCount),
          reinterpret_cast<int64_t **>(&componentIndices),
          &componentValues
        );

        elements.resize(elementsCount+componentValuesCount);
        for (size_t k(0); k < componentValuesCount; ++k) {
          // translate index within component tensor to TensorUnion index
          elements[elementsCount+k].first = getIndex(i, componentIndices[k]);
          elements[elementsCount+k].second = componentValues[k];
        }
        elementsCount += componentValuesCount;
        free(componentIndices);
        free(componentValues);
      }
      return elements;
    }

    /**
     * \brief Writes the given values together with their
     * respective indices. The indices are between 0 and getDimension()-1.
     **/
    void write(const std::vector<std::pair<size_t,F>> &elements) {
      // vectors to contain indices and values for each component tensor
      std::vector<std::vector<size_t>> tensorIndices(componentTensors.size());
      std::vector<std::vector<F>> tensorValues(componentTensors.size());

      for (size_t k(0); k < elements.size(); ++k) {
        size_t component;
        size_t componentIndex;
        fromIndex(elements[k].first, component, componentIndex);
        // write to respective component tensor
        tensorIndices[component].push_back(componentIndex);
        tensorValues[component].push_back(elements[k].second);
      }

      // write data of each tensor
      for (size_t i(0); i < componentTensors.size(); ++i) {
        tensorIndices[i].reserve(tensorIndices[i].size()+1);
        tensorValues[i].reserve(tensorIndices[i].size()+1);
        get(i)->evaluate()->write(
          tensorIndices[i].size(),
          reinterpret_cast<int64_t *>(tensorIndices[i].data()),
          tensorValues[i].data()
        );
      }
    }

  protected:
    /**
     * \brief The end of the TensorUnion index range for each component.
     * This vector is used for translating component number and indices
     * into TensorUnion indicies.
     **/
    std::vector<size_t> indexEnds;

    /**
     * \Brief Builds the index ends vector needed for the
     * index translation methods getIndex and fromIndex.
     **/
    void buildIndexTranslation() {
      indexEnds.resize(componentTensors.size());
      size_t indexBase(0);
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t tensorIndexSize(1);
        auto componentTensor(get(i)->inspect());
        for (size_t d(0); d < componentTensor->getLens().size(); ++d) {
          tensorIndexSize *= componentTensor->getLen(d);
        }
        indexEnds[i] = indexBase += tensorIndexSize;
      }
    }

    /**
     * \brief Sets this TensorUnion's component tensors by copying the given
     * component tensors. Called by copy constructors and copy assignments.
     **/
    void copyComponents(
      const std::vector<Ptr<TensorExpression<F,TE>>> &components
    ) {
      componentTensors.resize(components.size());
      for (size_t i(0); i < components.size(); ++i) {
        auto componentTensor(components[i]->inspect());
        // create tensor of identical shape, NOTE: no data is copied yet
        componentTensors[i] = Tcc<TE>::template tensor<F>(
          componentTensor->getLens(), componentTensor->getName()
        );
        // copy data
        COMPILE(
          (*componentTensors[i])[componentIndices[i]] <<=
            (*components[i])[componentIndices[i]]
        )->execute();
      }
    }

    /**
     * \brief Check if two TensorUnions are transpose of each other by swapping
     * the first and the second half of the component indices.
     **/
    // TODO: Improve speed?
    void checkDualCompatibility(const TensorUnion &a) const {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size() ; i++) {
        size_t indexLens(a.get(i)->order());
        for (size_t j(0); j < indexLens; j++) {
          size_t indexPos( get(i).find(a.getIndicies(i)[j]) );
          if (indexPos == std::string::npos) {
            throw New<Exception>(
              "Indices of tensor unions do not match", SOURCE_LOCATION
            );
          }
          if (a.get(i)->lens[j] != get(i)->lens[indexPos]) {
            throw New<Exception>(
              "Shapes of component tensors does not match", SOURCE_LOCATION
            );
          }
        }
      }
    }

    void checkCompatibilityTo(const TensorUnion &a) const {
      if (
        componentTensors.size() != a.componentTensors.size() ||
        componentIndices.size() != a.componentIndices.size()
      ) {
        throw New<Exception>(
          "Number of component tensors does no match", SOURCE_LOCATION
        );
      }
      // TODO: check shapes.
    }

  public:
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    ) {
      auto pointerNode(node->toPtr<AtomicNode<Ptr<Object>>>());
      if (!pointerNode) return nullptr;
      auto t( dynamicPtrCast<TensorUnion<F,TE>>(pointerNode->value) );
      if (!t) return nullptr;
      auto componentsNode(New<MapNode>(SOURCE_LOCATION));
      auto componentsNodePath(nodePath + ".components");
      for (size_t i(0); i < t->componentTensors.size(); ++i) {
        auto componentNode(New<MapNode>(SOURCE_LOCATION));
        componentNode->setValue("indices", t->componentIndices[i]);
        auto tensorNode(
          New<PointerNode<TensorExpression<F,TE>>>(
            t->componentTensors[i], SOURCE_LOCATION
          )
        );
        componentNode->get("tensor") = TensorIo::write(
          tensorNode,
          componentsNodePath + "." + to_string(i),
          useBinary
        );
        componentsNode->get(i) = componentNode;
      }
      auto writtenNode(New<MapNode>(SOURCE_LOCATION));
      writtenNode->get("components") = componentsNode;
      return writtenNode;
    }

    static Ptr<Node> read(
      const Ptr<MapNode> &node, const std::string &nodePath
    ) {
      auto componentsNode(node->getMap("components"));
      auto componentsNodePath(nodePath + ".components");
      std::vector<Ptr<TensorExpression<F,TE>>> tensors;
      std::vector<std::string> indices;
      for (size_t i(0); componentsNode->get(i); ++i) {
        auto tensorNode(
          TensorIo::read(
            componentsNode->getMap(i)->getMap("tensor"),
            componentsNodePath + "." + to_string(i)
          )
        );
        auto pointerNode(tensorNode->toPtr<AtomicNode<Ptr<Object>>>());
        tensors.push_back(
          dynamicPtrCast<TensorExpression<F,TE>>(pointerNode->value)
        );
        indices.push_back(
          componentsNode->getMap(i)->getValue<std::string>("indices")
        );
        OUT() << i << " tensor=" << tensors.back() << endl;
        OUT() << i << " indices=" << indices.back() << endl;
      }
      auto tensorUnion(New<TensorUnion<F,TE>>(tensors,indices));
      OUT() << "tensor union=" << tensorUnion << endl;
      return New<AtomicNode<Ptr<const TensorUnion<F,TE>>>>(
        tensorUnion, node->sourceLocation
      );
    }

    class TensorUnionIo;
    friend class TensorUnionIo;
  };

  class TensorUnionIo {
  public:
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    ) {
      // multiplex different tensor types
      Ptr<Node> writtenNode;
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        writtenNode = TensorUnion<Real<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = TensorUnion<Complex<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
      } else {
        using TE = DefaultDryTensorEngine;
        writtenNode = TensorUnion<Real<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = TensorUnion<Complex<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
      }
      return nullptr;
    }

    static Ptr<Node> read(
      const Ptr<MapNode> &node, const std::string &nodePath
    ) {
      // assume at least one component
      auto scalarType(
        node->getMap("components")->getMap(0)->getMap(
          "tensor"
        )->getValue<std::string>("scalarType")
      );
      // multiplex different tensor types
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        if (scalarType == TypeTraits<Real<>>::getName()) {
          return TensorUnion<Real<>,TE>::read(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<>>::getName()) {
          return TensorUnion<Complex<>,TE>::read(node, nodePath);
        }
      } else {
        using TE = DefaultDryTensorEngine;
        if (scalarType == TypeTraits<Real<>>::getName()) {
          return TensorUnion<Real<>,TE>::read(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<>>::getName()) {
          return TensorUnion<Complex<>,TE>::read(node, nodePath);
        }
      }
      std::stringstream explanation;
      explanation << "Unsupported sclarType '" << scalarType << "'" << endl;
      throw New<Exception>(explanation.str(), SOURCE_LOCATION);
    }

    static int WRITE_REGISTERED, READ_REGISTERED;
  };

  /**
   * \brief Returns the sum of two TensorUnions a and b, where
   * neither a nor b are modified.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> operator +(
    const TensorUnion<F,TE> &a, const TensorUnion<F,TE> &b
  ) {
    TensorUnion<F,TE> result(a);
    result += b;
    return std::move(result);
  }
  /**
   * \brief Returns the sum of two TensorUnions a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator +(
    TensorUnion<F,TE> &&a, const TensorUnion<F,TE> &b
  ) {
    a += b;
    return std::move(a);
  }
  /**
   * \brief Returns the sum of two TensorUnions a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator +(
    TensorUnion<F,TE> &a, const TensorUnion<F,TE> &&b
  ) {
    b += a;
    return std::move(b);
  }

  /**
   * \brief Returns the difference between two TensorUnions a and b, where
   * neither a nor b are modified.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> operator -(
    const TensorUnion<F,TE> &a, const TensorUnion<F,TE> &b
  ) {
    TensorUnion<F,TE> result(a);
    result -= b;
    return std::move(result);
  }
  /**
   * \brief Returns the difference between two TensorUnions a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator -(
    TensorUnion<F,TE> &&a, const TensorUnion<F,TE> &b
  ) {
    a -= b;
    return std::move(a);
  }
  /**
   * \brief Returns the difference between two TensorUnions a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator -(
    const TensorUnion<F,TE> &a, TensorUnion<F,TE> &&b
  ) {
    b -= a;
    // TODO: directly invoke sum to prevent extra multiplication by -1
    b *= F(-1);
    return std::move(b);
  }

  /**
   * \brief Returns the scalar multiple of the TensorUnion a
   * right-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> operator *(const TensorUnion<F,TE> &a, const F s) {
    TensorUnion<F,TE> result(a);
    result *= s;
    return std::move(result);
  }
  /**
   * \brief Returns the scalar multiple of the TensorUnion a
   * right-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator *(TensorUnion<F,TE> &&a, const F s) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Returns the scalar multiple of the TensorUnion a
   * left-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> operator *(const F s, const TensorUnion<F,TE> &a) {
    TensorUnion<F,TE> result(a);
    result *= s;
    return result;
  }
  /**
   * \brief Returns the scalar multiple of the TensorUnion a
   * left-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F, typename TE>
  inline TensorUnion<F,TE> &&operator *(const F s, TensorUnion<F,TE> &&a) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Writes the TensorUnion a to the given stream and returns it
   * for further stream operations.
   **/
  template <typename F, typename TE>
  inline std::ostream &operator <<(
    std::ostream &stream, const TensorUnion<F,TE> &a
  ) {
    stream << "( ";
    stream << a.get(0) << "[" << a.getIndices(0) << "]";
    for (size_t i(1); i < a.componentTensors.size(); ++i) {
      stream << ", " << a.get(i) << "[" << a.getIndices(i) << "]";
    }
    return stream << " )";
  }

}

#endif

