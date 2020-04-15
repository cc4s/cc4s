/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DATA_DEFINED
#define DATA_DEFINED

#include <util/Log.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
// TODO: find out why Exception must be included after string,map and ctf
#include <util/Exception.hpp>
#include <util/SharedPointer.hpp>

#include <string>
#include <map>
#include <vector>
#include <sstream>

namespace cc4s {
  class MapNode;
  class SymbolNode;
  template <typename AtomicType> class AtomicNode;

  template <typename F>
  class TypeTraits;

  class Node: public Thisable<Node> {
  public:
    virtual ~Node() {
    }
    virtual bool isAtomic() {
      return true;
    }
    virtual Ptr<Node> &get(const std::string &element) {
      Assert(false, "atomic node does not have element " + element);
    }
    virtual Ptr<Node> &get(const size_t element) {
      Assert(false, "atomic node does not have element " + element);
    }
    virtual std::string toString() = 0;
    // provide convenience cast routines
    Ptr<MapNode> map() {
      return std::dynamic_pointer_cast<MapNode>(this->toPtr<Node>());
    }
    Ptr<SymbolNode> symbol() {
      return std::dynamic_pointer_cast<SymbolNode>(this->toPtr<Node>());
    }
    template <typename AtomicType>
    Ptr<AtomicNode<AtomicType>> atom() {
      return std::dynamic_pointer_cast<AtomicNode<AtomicType>>(
        this->toPtr<Node>()
      );
    }
    std::string comment;
  };

  class SymbolNode: public Node {
  public:
    /**
     * \brief Constructor for symbol nodes.
     */
    SymbolNode(const std::string &value_): value(value_) {
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
    AtomicNode(const AtomicType &value_): value(value_) {
    }
    std::string toString() override {
      std::stringstream stream;
      stream << value;
      return stream.str();
    }
    AtomicType value;
  };

  class MapNode: public Node {
  public:
    MapNode() {
    }
    bool isAtomic() override {
      return false;
    }
    Ptr<Node> &get(const std::string &element) override {
      return elements[element];
    }
    Ptr<Node> &get(const size_t element) override {
      return elements[std::to_string(element)];
    }
    std::string toString() override {
      // FIXME: implement
      return "{...}";
    }
/*
    // TODO: proper key iterators
    std::vector<std::string> getKeys() const {
      std::vector<const std::string> keys();
      keys.reserve(elements.size());
      for (auto iterator: elements) {
        keys.push_back(iterator->first);
      }
      return keys;
    }
*/
    size_t size() const {
      return elements.size();
    }

    // convenience member access
    std::string getSymbol(const std::string &element) {
      Assert(get(element), "expecting key '" + element + "'");
      auto symbolNode(get(element)->symbol());
      Assert(symbolNode, "expecting '" + element + "' to be a symbol");
      return symbolNode->value;
    }
    template <typename Target>
    Target getValue(const std::string &element) {
      Assert(get(element), "expecting key '" + element + "'");
      // first, try to convert to expected type node
      auto targetAtomNode(get(element)->atom<Target>());
      if (targetAtomNode) {
        return targetAtomNode->value;
      } else {
        // otherwise, try to convert to string and then to requested target
        std::stringstream stream(get(element)->toString());
        Target targetValue;
        stream >> targetValue;
        // FIXME: properly throw exception when nothing meaningful is done
/*
        Assert(
          stream.str().size() < get(element)->toString().size(),
           "failed to convert '" + get(element)->toString() + "' to "
           + TypeTraits<Target>::getName() + ", left: " + stream.str()
        );
*/
        if (stream.str().size() > 0) {
          LOG(0, "getValue")
            << "WARNING: not all data used converting '"
            << get(element)->toString() << "' to "
            << TypeTraits<Target>::getName() << std::endl;
        }
        // if successful replace previous node with converted one
        get(element) = New<AtomicNode<Target>>(targetValue);
        return targetValue;
      }
    }
    template <typename Target>
    Target getValue(const std::string &element, const Target &defaultValue) {
      if (get(element)) {
        // if key present, proceed
        return getValue<Target>(element);
      } else {
        // otherwise, enter default value in map
        get(element) = New<AtomicNode<Target>>(defaultValue);
        return defaultValue;
      }
    }
    Ptr<MapNode> getMap(const std::string &element) {
      Assert(get(element), "expecting key '" + element + "'");
      auto mapNode(get(element)->map());
      Assert(mapNode, "expecting '" + element + "' to be a map");
      return mapNode;
    }
    Ptr<MapNode> getMap(const size_t element) {
      return getMap(std::to_string(element));
    }
  protected:
    std::map<std::string,Ptr<Node>> elements;
  };


  // root node where all data is stored
  static Ptr<MapNode> storage;

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

