/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DATA_DEFINED
#define DATA_DEFINED

#include <util/Log.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
// TODO: find out why Exception must be included after string,map and ctf
#include <util/Exception.hpp>
#include <util/SharedPointer.hpp>

// tensor engine selection
#include <engines/DryTensorEngine.hpp>
#include <engines/CtfTensorEngine.hpp>
namespace cc4s {
  typedef cc4s::CtfTensorEngine DefaultTensorEngine;
}

#include <string>
#include <map>
#include <sstream>

namespace cc4s {
  class Node: public Thisable<Node> {
  public:
    virtual ~Node() {
    }
    virtual bool isAtomic() {
      return true;
    }
    virtual Ptr<Node> get(const std::string &element) {
      Assert(false, "atomic node does not have element " + element);
    }
    virtual Ptr<Node> get(const size_t element) {
      Assert(false, "atomic node does not have element " + element);
    }
    template <typename AtomicType>
    AtomicType value() {
      return std::dynamic_pointer_cast<AtomicNode<AtomicType>>(this)->value;
    }
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
      return elements["" + element];
    }
  protected:
    std::map<std::string,Ptr<Node>> elements;
  };

  template <typename AtomicType>
  class AtomicNode: public Node {
  public:
    /**
     * \brief Constructor for atomic nodes.
     */
    AtomicNode(const AtomicType &value_): value(value_) {
    }
    AtomicType value;
  };

  // root node where all data is stored
  static Ptr<MapNode> root;


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
    static std::string getName() { return "real<64>"; }
  };
  template <>
  class TypeTraits<Complex<64>> {
  public:
    static std::string getName() { return "complex<64>"; }
  };
  template <>
  class TypeTraits<Real<128>> {
  public:
    static std::string getName() { return "real<128>"; }
  };
  template <>
  class TypeTraits<Complex<128>> {
  public:
    static std::string getName() { return "complex<128>"; }
  };
  template <F,TE>
  class TypeTraits<Tensor<F,TE>> {
  public:
    static std::string getName() {
      return "tensor of " + TypeTraits<F>::getName();
    }
  };
}

#endif

