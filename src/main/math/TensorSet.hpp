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

#ifndef TENSOR_SET_DEFINED
#define TENSOR_SET_DEFINED

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
   * \brief Represents a set of tensors and provides
   * vector space operations: addition, scalar multiplication, inner product,
   * and complex conjugation.
   **/
  template <typename F, typename TE>
  class TensorSet: public Object {
  public:
    typedef F FieldType;

    std::vector<Ptr<TensorExpression<F,TE>>> componentTensors;
    std::vector<std::string> componentIndices;

    /**
     * \brief Default constructor for an empty tensor set without elements.
     **/
    TensorSet() {
    }

    /**
     * \brief Move constructor taking possession of the tensors owned by a.
     **/
    TensorSet(
      TensorSet &&a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices)
    {
    }

    /**
     * \brief Copy constructor copying the tensors owned by a.
     **/
    TensorSet(
      const TensorSet &a
    ):
      componentTensors(a.componentTensors.size()),
      componentIndices(a.componentIndices)
    {
      copyComponents(a.componentTensors);
    }

    /**
     * \brief Move constructor taking possession of the tensors given.
     **/
    TensorSet(
      const std::vector<Ptr<TensorExpression<F,TE>>> &tensors,
      const std::vector<std::string> &indices
    ):
      componentTensors(tensors),
      componentIndices(indices)
    {
    }

    /**
     * \brief Move constructor taking possession of the tensors given
     * by the iterators.
     **/
    template <typename TensorsIterator, typename IndicesIterator>
    TensorSet(
      TensorsIterator tensorsBegin, TensorsIterator tensorsEnd,
      IndicesIterator indicesBegin, IndicesIterator indicesEnd
    ):
      componentTensors(tensorsBegin, tensorsEnd),
      componentIndices(indicesBegin, indicesEnd)
    {
    }

    /**
     * \brief Retrieves the i-th component tensor. Note that
     * the Tensor is not const since rearrangement may be
     * required also in non-modifying tensor operations.
     **/
    const Ptr<TensorExpression<F,TE>> &get(const Natural<> i) const {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component tensor.
     **/
    Ptr<TensorExpression<F,TE>> &get(const Natural<> i) {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component indices.
     **/
    const std::string &getIndices(const Natural<> i) const {
      return componentIndices[i];
    }

    /**
     * \brief Retrieves the i-th component indices as modifiable string.
     **/
    std::string &getIndices(const Natural<> i) {
      return componentIndices[i];
    }

    /**
     * \brief Move assignment operator taking possession of the tensors
     * owned by a.
     **/
    TensorSet &operator =(const TensorSet &&a) {
      componentTensors = a.componentTensors;
      componentIndices = a.componentIndices;
      return *this;
    }

    /**
     * \brief Copy assignment operator copying the tensors owned by a.
     **/
    TensorSet &operator =(const TensorSet &a) {
      componentIndices = a.componentIndices;
      copyComponents(a.componentTensors);
      return *this;
    }

    /**
     * \brief Add-to assignment operator adding each component of a
     * to the respective component of this TensorSet.
     **/
    TensorSet &operator += (const TensorSet &a) {
      ASSERT(isCompatibleTo(a), "Incompatible TensorSets");
      for (Natural<> i(0); i < componentTensors.size(); ++i) {
        COMPILE(
          (*get(i))[getIndices(i)] += (*a.get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Subtract-from assignment operator subtracting each component of a
     * from the respective component of this TensorSet.
     **/
    TensorSet &operator -= (const TensorSet &a) {
      ASSERT(isCompatibleTo(a), "Incompatible TensorSets.");
      for (Natural<> i(0); i < componentTensors.size(); ++i) {
        COMPILE(
          (*get(i))[getIndices(i)] -= (*a.get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Multiply-by assignment operator scalar multiplying each component
     * each component of this TensorSet by the given scalar.
     **/
    TensorSet &operator *= (const F s) {
      for (Natural<> i(0); i < componentTensors.size(); ++i) {
        COMPILE(
          (*get(i))[getIndices(i)] <<= s * (*get(i))[getIndices(i)]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Returns the inner product of this ket-TensorSet with the
     * given ket-TensorSet a. The elements of this TensorSet are conjugated
     * in the inner product.
     **/
    F dot(const TensorSet &a) const {
      ASSERT(isCompatibleTo(a), "Incompatible TensorSets");
      auto result( Tcc<TE>::template tensor<F>("dot") );
      for (Natural<> i(0); i < componentTensors.size(); ++i) {
        // add to result
        COMPILE(
          (*result)[""] += (*get(i))[getIndices(i)] *
            map<F>(cc4s::conj<F>, (*a.get(i))[getIndices(i)])
        )->execute();
      }
      return result->read();
    }

    /**
     * \brief Get the number of component tensors of this TensorSet.
     */
    Natural<> getComponentsCount() const {
      return componentTensors.size();
    }

    bool isCompatibleTo(const TensorSet &a) const {
      return
        componentTensors.size() == a.componentTensors.size() &&
        componentIndices.size() == a.componentIndices.size();
      // TODO: check shapes.
    }

  protected:
    /**
     * \brief Sets this TensorSet's component tensors by copying the given
     * component tensors. Called by copy constructors and copy assignments.
     **/
    void copyComponents(
      const std::vector<Ptr<TensorExpression<F,TE>>> &components
    ) {
      componentTensors.resize(components.size());
      for (Natural<> i(0); i < components.size(); ++i) {
        auto componentTensor(components[i]->inspect());
        // create tensor of identical shape, NOTE: no data is copied yet
        componentTensors[i] = Tcc<TE>::template tensor<F>(
          componentTensor->getName()
        );
        // copy data
        COMPILE(
          (*componentTensors[i])[componentIndices[i]] <<=
            (*components[i])[componentIndices[i]]
        )->execute();
      }
    }

  public:
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    ) {
      auto pointerNode(node->toPtr<AtomicNode<Ptr<Object>>>());
      if (!pointerNode) return nullptr;
      auto t( dynamicPtrCast<TensorSet<F,TE>>(pointerNode->value) );
      if (!t) return nullptr;
      auto componentsNode(New<MapNode>(SOURCE_LOCATION));
      auto componentsNodePath(nodePath + ".components");
      for (Natural<> i(0); i < t->componentTensors.size(); ++i) {
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
      for (Natural<> i(0); componentsNode->get(i); ++i) {
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
      auto tensorSet(New<TensorSet<F,TE>>(tensors,indices));
      OUT() << "tensor set=" << tensorSet << endl;
      return New<AtomicNode<Ptr<const TensorSet<F,TE>>>>(
        tensorSet, node->sourceLocation
      );
    }

    class TensorSetIo;
    friend class TensorSetIo;
  };

  class TensorSetIo {
  public:
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    ) {
      // multiplex different tensor types
      Ptr<Node> writtenNode;
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        writtenNode = TensorSet<Real<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = TensorSet<Complex<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
      } else {
        using TE = DefaultDryTensorEngine;
        writtenNode = TensorSet<Real<64>,TE>::write(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = TensorSet<Complex<64>,TE>::write(node, nodePath, useBinary);
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
          return TensorSet<Real<>,TE>::read(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<>>::getName()) {
          return TensorSet<Complex<>,TE>::read(node, nodePath);
        }
      } else {
        using TE = DefaultDryTensorEngine;
        if (scalarType == TypeTraits<Real<>>::getName()) {
          return TensorSet<Real<>,TE>::read(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<>>::getName()) {
          return TensorSet<Complex<>,TE>::read(node, nodePath);
        }
      }
      std::stringstream explanation;
      explanation << "Unsupported sclarType '" << scalarType << "'" << endl;
      throw New<Exception>(explanation.str(), SOURCE_LOCATION);
    }

    static bool WRITE_REGISTERED, READ_REGISTERED;
  };

  /**
   * \brief Returns the sum of two TensorSets a and b, where
   * neither a nor b are modified.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> operator +(
    const TensorSet<F,TE> &a, const TensorSet<F,TE> &b
  ) {
    TensorSet<F,TE> result(a);
    result += b;
    return std::move(result);
  }
  /**
   * \brief Returns the sum of two TensorSets a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator +(
    TensorSet<F,TE> &&a, const TensorSet<F,TE> &b
  ) {
    a += b;
    return std::move(a);
  }
  /**
   * \brief Returns the sum of two TensorSets a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator +(
    TensorSet<F,TE> &a, const TensorSet<F,TE> &&b
  ) {
    b += a;
    return std::move(b);
  }

  /**
   * \brief Returns the difference between two TensorSets a and b, where
   * neither a nor b are modified.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> operator -(
    const TensorSet<F,TE> &a, const TensorSet<F,TE> &b
  ) {
    TensorSet<F,TE> result(a);
    result -= b;
    return std::move(result);
  }
  /**
   * \brief Returns the difference between two TensorSets a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator -(
    TensorSet<F,TE> &&a, const TensorSet<F,TE> &b
  ) {
    a -= b;
    return std::move(a);
  }
  /**
   * \brief Returns the difference between two TensorSets a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator -(
    const TensorSet<F,TE> &a, TensorSet<F,TE> &&b
  ) {
    b -= a;
    // TODO: directly invoke sum to prevent extra multiplication by -1
    b *= F(-1);
    return std::move(b);
  }

  /**
   * \brief Returns the scalar multiple of the TensorSet a
   * right-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> operator *(const TensorSet<F,TE> &a, const F s) {
    TensorSet<F,TE> result(a);
    result *= s;
    return std::move(result);
  }
  /**
   * \brief Returns the scalar multiple of the TensorSet a
   * right-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator *(TensorSet<F,TE> &&a, const F s) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Returns the scalar multiple of the TensorSet a
   * left-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> operator *(const F s, const TensorSet<F,TE> &a) {
    TensorSet<F,TE> result(a);
    result *= s;
    return result;
  }
  /**
   * \brief Returns the scalar multiple of the TensorSet a
   * left-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F, typename TE>
  inline TensorSet<F,TE> &&operator *(const F s, TensorSet<F,TE> &&a) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Writes the TensorSet a to the given stream and returns it
   * for further stream operations.
   **/
  template <typename F, typename TE>
  inline std::ostream &operator <<(
    std::ostream &stream, const TensorSet<F,TE> &a
  ) {
    stream << "( ";
    stream << a.get(0) << "[" << a.getIndices(0) << "]";
    for (Natural<> i(1); i < a.componentTensors.size(); ++i) {
      stream << ", " << a.get(i) << "[" << a.getIndices(i) << "]";
    }
    return stream << " )";
  }

}

#endif

