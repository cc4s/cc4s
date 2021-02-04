/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DATA_DEFINED
#define DATA_DEFINED

#include <util/Log.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <util/SourceLocation.hpp>
#include <util/Exception.hpp>
#include <util/SharedPointer.hpp>

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>

namespace cc4s {
  class MapNode;
  class SymbolNode;
  template <typename AtomicType> class AtomicNode;

  template <typename F>
  class TypeTraits;

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
    // provide convenience cast routines
    Ptr<MapNode> toMap() {
      return std::dynamic_pointer_cast<MapNode>(this->toPtr<Node>());
    }
    Ptr<SymbolNode> toSymbol() {
      return std::dynamic_pointer_cast<SymbolNode>(this->toPtr<Node>());
    }
    template <typename AtomicType>
    Ptr<AtomicNode<AtomicType>> toAtom() {
      return std::dynamic_pointer_cast<AtomicNode<AtomicType>>(
        this->toPtr<Node>()
      );
    }
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
    Ptr<Node> &get(const size_t key) {
      return elements[std::to_string(key)];
    }
    // TODO: proper key iterators
    std::vector<std::string> getKeys() const {
      std::vector<std::string> keys;
      keys.reserve(elements.size());
      for (auto pairs: elements) {
        keys.push_back(pairs.first);
      }
      return keys;
    }
    size_t size() const {
      return elements.size();
    }

    // convenience member access
    std::string getSymbol(const std::string &key) {
      ASSERT_LOCATION(
        get(key), "expecting key '" + key + "'", sourceLocation
      );
      auto symbolNode(get(key)->toSymbol());
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
      ASSERT_LOCATION(get(key), "expecting key '" + key + "'", sourceLocation);
      // first, try to convert to expected type node
      auto targetAtomNode(get(key)->toAtom<Target>());
      if (targetAtomNode) {
        return targetAtomNode->value;
      } else {
        // otherwise, try to convert to string and then to requested target
        std::stringstream stream(get(key)->toString());
        Target targetValue;
        stream >> targetValue;
        // FIXME: properly throw exception when nothing meaningful is done
/*
        Assert(
          stream.str().size() < get(key)->toString().size(),
           "failed to convert '" + get(key)->toString() + "' to "
           + TypeTraits<Target>::getName() + ", left: " + stream.str()
        );
        if (stream.str().size() > 0) {
          LOG()
            << "WARNING: not all data used converting '"
            << get(key)->toString() << "' to "
            << TypeTraits<Target>::getName() << std::endl;
        }
*/
        // if successful replace previous node with converted one
        get(key) = New<AtomicNode<Target>>(targetValue, sourceLocation);
        return targetValue;
      }
    }
    template <typename Target>
    Target getValue(const size_t index) {
      return getValue<Target>(std::to_string(index));
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
    Target getValue(const size_t index, const Target &defaultValue) {
      return getValue<Target>(std::to_string(index), defaultValue);
    }
    template <typename Target>
    void setValue(
      const std::string &key, const Target &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      get(key) = New<AtomicNode<Target>>(value, sourceLocation);
    }
    template <typename Target>
    void setValue(
      const size_t index, const Target &value,
      const SourceLocation &sourceLocation = SOURCE_LOCATION
    ) {
      get(index) = New<AtomicNode<Target>>(value, sourceLocation);
    }

    Ptr<MapNode> getMap(const std::string &element) {
      ASSERT_LOCATION(
        get(element), "expecting key '" + element + "'", sourceLocation
      );
      auto mapNode(get(element)->toMap());
      ASSERT_LOCATION(
        mapNode, "expecting '" + element + "' to be a map", sourceLocation
      );
      return mapNode;
    }
    Ptr<MapNode> getMap(const size_t element) {
      return getMap(std::to_string(element));
    }

    bool isGiven(const std::string &element) {
      if (get(element)) { return true; }
      else { return false; }
    }

    void push_back(const Ptr<Node> &node) {
      get(size()) = node;
    }
  protected:
    std::map<std::string,Ptr<Node>> elements;
  };

  /**
   * Traits class for tensor element types used in cc4s.
   * It provides type specific information such as type name to
   * be displayed to the user.
   */
  template <typename F>
  class TypeTraits;

  template <>
  class TypeTraits<std::string> {
  public:
    static std::string getName() { return "text"; }
  };
  template <>
  class TypeTraits<bool> {
  public:
    static std::string getName() { return "boolean"; }
  };
  template <>
  class TypeTraits<int64_t> {
  public:
    static std::string getName() { return "integer"; }
  };
  template <>
  class TypeTraits<Real<64>> {
  public:
    static std::string getName() { return "real64"; }
  };
  template <>
  class TypeTraits<Complex<64>> {
  public:
    static std::string getName() { return "complex64"; }
  };
  template <>
  class TypeTraits<Real<128>> {
  public:
    static std::string getName() { return "real128"; }
  };
  template <>
  class TypeTraits<Complex<128>> {
  public:
    static std::string getName() { return "complex128"; }
  };
/*
  template <F,TE>
  class TypeTraits<Tensor<F,TE>> {
  public:
    static std::string getName() {
      return "tensor of " + TypeTraits<F>::getName();
    }
  };
*/
}

#endif

