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

#ifndef NODE_DEFINED
#define NODE_DEFINED

#include <Object.hpp>
#include <Integer.hpp>
#include <SourceLocation.hpp>
#include <Exception.hpp>
#include <SharedPointer.hpp>
#include <TypeTraits.hpp>

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>

namespace cc4s {
  // forward declarations
  class MapNode;
  class SymbolNode;
  template <typename AtomicType> class AtomicNode;
  template <typename PointedType> class PointerNode;

  class Node: public Thisable<Node> {
  public:
    Node(
      const SourceLocation &sourceLocation_
    ): sourceLocation(sourceLocation_) {
    }

    virtual ~Node() {
    }
    virtual bool isAtomic() {
      return true;
    }
    virtual std::string toString() = 0;

    std::string comment;
    SourceLocation sourceLocation;
  };

  class SymbolNode: public Node {
  public:
    /**
     * \brief Constructor for symbol nodes.
     */
    SymbolNode(
      const std::string &value_, const SourceLocation &sourceLocation_
    ): Node(sourceLocation_), value(value_) {
    }
    std::string toString() override {
      return value;
    }
    std::string value;
  };

  template <typename AtomicType>
  class AtomicNode: public Node {
  public:
    /**
     * \brief Constructor for atomic nodes.
     */
    AtomicNode(
      const AtomicType &value_, const SourceLocation &sourceLocation_
    ): Node(sourceLocation_), value(value_) {
    }
    std::string toString() override {
      std::stringstream stream;
      // FIXME: fixed precision
      stream << std::setprecision(17) << value;
      return stream.str();
    }
    AtomicType value;
  };

  template <typename PointedType>
  class PointerNode: public AtomicNode<Ptr<Object>> {
  public:
    PointerNode(
      const Ptr<PointedType> &value_, const SourceLocation &sourceLocation_
    ): AtomicNode(value_, sourceLocation_) {
    }
  };

  class MapNode: public Node {
  public:
    MapNode(const SourceLocation &sourceLocation_): Node(sourceLocation_) {
    }
    bool isAtomic() override {
      return false;
    }
    std::string toString() override {
      // FIXME: implement
      return "{...}";
    }
    Ptr<Node> &get(const std::string &key) {
      return elements[key];
    }
    Ptr<Node> &get(const Natural<> key) {
      return elements[std::to_string(key)];
    }

    std::vector<std::string> getKeys() const {
      std::vector<std::string> keys;
      keys.reserve(elements.size());
      for (auto pair: elements) {
        if (pair.second) {
          keys.push_back(pair.first);
        }
      }
      return keys;
    }

    Natural<> getSize() const {
      Natural<> size(0);
      for (auto pair: elements) {
        if (pair.second) ++size;
      }
      return size;
    }

    // convenience member access
    std::string getSymbol(const std::string &key) {
      ASSERT_LOCATION(
        get(key), "expecting key '" + key + "'", sourceLocation
      );
      auto symbolNode(get(key)->toPtr<SymbolNode>());
      ASSERT_LOCATION(
        symbolNode, "expecting '" + key + "' to be a symbol", sourceLocation
      );
      return symbolNode->value;
    }
    void setSymbol(
      const std::string &key, const std::string &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      get(key) = New<SymbolNode>(value, sourceLocation);
    }
    template <typename Target>
    Target getValue(const std::string &key) {
      static_assert(
        TypeTraits<Target>::IS_PRIMITIVE,
        "Cannot use getValue<Target> with non-standard atomic type Target. "
        "See TypeTraits.hpp for standard atomic types. "
        "Use getPtr<Type> for pointers."
      );
      ASSERT_LOCATION(get(key), "expecting key '" + key + "'", sourceLocation);
      // first, try to convert to expected type node
      auto targetAtomNode(get(key)->toPtr<AtomicNode<Target>>());
      if (targetAtomNode) {
        return targetAtomNode->value;
      } else {
        // otherwise, try to convert to string and then to requested target
        std::stringstream stream(get(key)->toString());
        Target targetValue;
        stream >> targetValue;
        ASSERT(
          stream.str().size() - stream.tellg() != 0,
           "failed to convert '" + get(key)->toString() + "' to "
           + TypeTraits<Target>::getName() + ", unconverted: "
           + stream.str().substr(stream.tellg())
        );
        // if successful replace previous node with converted one
        get(key) = New<AtomicNode<Target>>(targetValue, sourceLocation);
        return targetValue;
      }
    }
    template <typename PointedTarget>
    Ptr<PointedTarget> getPtr(const std::string &key) {
      ASSERT_LOCATION(get(key), "expecting key '" + key + "'", sourceLocation);
      // first, try to convert to void pointer node
      Ptr<Node> node(get(key));
      Ptr<AtomicNode<Ptr<Object>>> ptrAtom(
        node->toPtr<AtomicNode<Ptr<Object>>>()
      );
      ASSERT_LOCATION(
        ptrAtom,
        "expecting object for key '" + key + "'", sourceLocation
      );
      // then, try to convert pointer to target type
      auto pointer(dynamicPtrCast<PointedTarget>(ptrAtom->value));
      // if it fails nullptr will be returned
      return pointer;
    }
    template <typename Target>
    Target getValue(const Natural<> index) {
      return getValue<Target>(std::to_string(index));
    }
    template <typename PointedTarget>
    PointedTarget getPtr(const Natural<> index) {
      return getPtr<PointedTarget>(std::to_string(index));
    }
    template <typename Target>
    Target getValue(const std::string &element, const Target &defaultValue) {
      if (get(element)) {
        // if key present, proceed
        return getValue<Target>(element);
      } else {
        // otherwise, enter default value in map
        get(element) = New<AtomicNode<Target>>(defaultValue, sourceLocation);
        return defaultValue;
      }
    }
    template <typename Target>
    Target getValue(const Natural<> index, const Target &defaultValue) {
      return getValue<Target>(std::to_string(index), defaultValue);
    }
    template <typename Target>
    void setValue(
      const std::string &key, const Target &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      static_assert(
        TypeTraits<Target>::IS_PRIMITIVE,
        "Cannot use setValue<Target> with non-primitive atomic type Target. "
        "See TypeTraits.hpp for primitive atomic types. "
        "Use setPtr<Type> for pointers."
      );
      get(key) = New<AtomicNode<Target>>(value, sourceLocation);
    }
    template <typename Target>
    void setValue(
      const Natural<> index, const Target &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      setValue(std::to_string(index), value, sourceLocation);
    }
    template <typename PointedTarget>
    void setPtr(
      const std::string &key, const Ptr<PointedTarget> &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      get(key) = New<PointerNode<PointedTarget>>(value, sourceLocation);
    }
    template <typename PointedTarget>
    void setPtr(
      const Natural<> index, const Ptr<PointedTarget> &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      get(index) = New<PointerNode<PointedTarget>>(value, sourceLocation);
    }

    Ptr<MapNode> getMap(const std::string &element) {
      ASSERT_LOCATION(
        get(element), "expecting key '" + element + "'", sourceLocation
      );
      auto mapNode(get(element)->toPtr<MapNode>());
      ASSERT_LOCATION(
        mapNode, "expecting '" + element + "' to be a map", sourceLocation
      );
      return mapNode;
    }
    Ptr<MapNode> getMap(const Natural<> element) {
      return getMap(std::to_string(element));
    }

    bool isGiven(const std::string &element) {
      return elements.find(element) != elements.end();
    }

    void push_back(const Ptr<Node> &node) {
      get(getSize()) = node;
    }

    std::vector<std::string> getKeys() {
      std::vector<std::string> keys;
      keys.reserve(elements.size());
      for (auto pair: elements) {
        if (pair.second) {
          keys.push_back(pair.first);
        }
      }
      return keys;
    }

  protected:
    std::map<std::string,Ptr<Node>> elements;
  };
}

#endif

