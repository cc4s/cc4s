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

#include <tcc/Tcc.hpp>
#include <MathFunctions.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>
#include <TensorIo.hpp>
#include <Log.hpp>

#include <map>
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

  protected:
    std::map<std::string,Ptr<TensorExpression<F,TE>>> components;

  public:

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
      components(a.components)
    {
    }

    /**
     * \brief Copy constructor copying the tensors owned by a.
     **/
    TensorSet(
      const TensorSet &a
    ) {
      copyComponents(a.components);
    }

    /**
     * \brief Move constructor taking possession of the tensors given.
     **/
    TensorSet(
      const std::map<std::string, Ptr<TensorExpression<F,TE>>> &components_
    ):
      components(components_)
    {
    }

    Natural<> getSize() const {
      Natural<> size(0);
      for (auto pair: components) {
        if (pair.second) ++size;
      }
      return size;
    }

    std::vector<std::string> getKeys() const {
      std::vector<std::string> keys;
      keys.reserve(components.size());
      for (auto component: components) {
        if (component.second) {
          keys.push_back(component.first);
        }
      }
      return keys;
    }

    /**
     * \brief Retrieves the requested component tensor addressed by the
     * given key string.
     * Note that the Tensor is not const since rearrangement may be
     * required also in non-modifying tensor operations.
     **/
    const Ptr<TensorExpression<F,TE>> get(const std::string &key) const {
      auto iterator(components.find(key));
      if (iterator == components.end()) {
        return nullptr;
      } else {
        return iterator->second;
      }
    }

    /**
     * \brief Retrieves the requested component tensor addressed by the
     * given key string.
     **/
    Ptr<TensorExpression<F,TE>> &get(const std::string &key) {
      return components[key];
    }

    /**
     * \brief Move assignment operator taking possession of the tensors
     * owned by a.
     **/
    TensorSet &operator =(const TensorSet &&a) {
      components = a.components;
      return *this;
    }

    /**
     * \brief Copy assignment operator copying the tensors owned by a.
     **/
    TensorSet &operator =(const TensorSet &a) {
      copyComponents(a.components);
      return *this;
    }

    std::string generateIndices(const Natural<> size) const {
      char indices[size+1];
      for (Natural<> i(0); i < size; ++i) {
        indices[i] = 'a' + i;
      }
      indices[size] = 0;
      return std::string(indices);
    }

    std::string generateIndices(const std::string &key) const {
      return generateIndices(get(key)->inspect()->getLens().size());
    }

    /**
     * \brief Add-to assignment operator adding each component of a
     * to the respective component of this TensorSet.
     **/
    TensorSet &operator += (const TensorSet &a) {
      ASSERT(isCompatibleTo(a), "Incompatible TensorSets");
      for (auto component: components) {
        auto key(component.first);
        auto tensorExpression(component.second);
        auto indices(generateIndices(key));
        COMPILE(
          (*tensorExpression)[indices] += (*a.get(key))[indices]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Subtract-from assignment operator subtracting each component of a
     * from the respective component of this TensorSet.
     **/
    TensorSet &operator -= (const TensorSet &a) {
      ASSERT(isCompatibleTo(a), "Incompatible TensorSets");
      for (auto component: components) {
        auto key(component.first);
        auto tensorExpression(component.second);
        auto indices(generateIndices(key));
        COMPILE(
          (*tensorExpression)[indices] -= (*a.get(key))[indices]
        )->execute();
      }
      return *this;
    }

    /**
     * \brief Multiply-by assignment operator scalar multiplying each component
     * each component of this TensorSet by the given scalar.
     **/
    TensorSet &operator *= (const F s) {
      for (auto component: components) {
        auto key(component.first);
        auto tensorExpression(component.second);
        auto indices(generateIndices(key));
        COMPILE(
          (*tensorExpression)[indices] <<= s * (*tensorExpression)[indices]
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
      for (auto component: components) {
        auto key(component.first);
        auto tensorExpression(component.second);
        auto indices(generateIndices(key));
        // add to result
        COMPILE(
          (*result)[""] += (*tensorExpression)[indices] *
            map<F>(cc4s::conj<F>, (*a.get(key))[indices])
        )->execute();
      }
      return result->read();
    }

    /**
     * \brief Get the number of component tensors of this TensorSet.
     */
    Natural<> getComponentsCount() const {
      return components.size();
    }

    bool isCompatibleTo(const TensorSet &a) const {
      return components.size() == a.components.size();
      // TODO: check shapes.
    }

  protected:
    /**
     * \brief Sets this TensorSet's component tensors by copying the given
     * component tensors. Called by copy constructors and copy assignments.
     **/
    void copyComponents(
      const std::map<std::string, Ptr<TensorExpression<F,TE>>> &sources
    ) {
      components.clear();
      for (auto source: sources) {
        auto sourceKey(source.first);
        auto sourceTensor(source.second->inspect());
        // create tensor receiving copy of sourceTensor
        auto component(
          Tcc<TE>::template tensor<F>(sourceTensor->getName())
        );
        // copy data
        auto indices(generateIndices(sourceTensor->getLens().size()));
        COMPILE(
          (*component)[indices] <<= (*sourceTensor)[indices]
        )->execute();
        // transfer dimension info, TODO: should be done in tcc
        component->dimensions = sourceTensor->dimensions;
        // transfer meta-data, TODO: should be done in tcc
        component->getMetaData() = sourceTensor->getMetaData();
        // enter in map
        components[sourceKey] = component;
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
      for (auto component: t->components) {
        auto tensorNode(
          New<PointerNode<TensorExpression<F,TE>>>(
            component.second, SOURCE_LOCATION
          )
        );
        auto componentNode(
          TensorIo::write(
            tensorNode,
            componentsNodePath + "." + component.first,
            useBinary
          )
        );
        componentsNode->get(component.first) = componentNode;
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
      std::map<std::string, Ptr<TensorExpression<F,TE>>> components;
      for (auto key: componentsNode->getKeys()) {
        auto tensorNode(
          TensorIo::read(
            componentsNode->getMap(key),
            componentsNodePath + "." + key
          )
        );
        auto pointerNode(tensorNode->toPtr<PointerNode<Object>>());
        components[key] =
          dynamicPtrCast<Tensor<F,TE>>(pointerNode->value);
      }
      auto tensorSet(New<TensorSet<F,TE>>(components));
      return New<PointerNode<TensorSet<F,TE>>>(
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
      if (!Cc4s::dryRun) {
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
      auto componentsNode(node->getMap("components"));
      // NOTE: assumes at least one component
      auto firstKey(componentsNode->getKeys()[0]);
      auto scalarType(
        componentsNode->getMap(firstKey)->getValue<std::string>("scalarType")
      );
      // multiplex different tensor types
      if (!Cc4s::dryRun) {
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
    return result;
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
    return result;
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
    return result;
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
    std::string delimiter("");
    stream << "( ";
    for (auto component: a.components) {
      stream << delimiter << component.first;
      delimiter = ", ";
    }
    return stream << " )";
  }

}

#endif

