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
  /**
   * Traits class for tensor element types used in cc4s.
   * It provides type specific information such as type name to
   * be displayed to the user.
   */
  template <typename F>
  class TypeTraits;

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

  class Data: public THISABLE(Data) {
  public:
    enum Stage {
      MENTIONED = 0, TYPED = 1, ALLOCATED = 2, 
      READY = 3,
      UNUSED = 4, LINGERING = 5
    };
    Data(const std::string &name_) {
      std::stringstream sStream;
      sStream << name_ << " of yet unknown type";
      typeName = sStream.str();
      // enter datum in global data map
      dataMap[name_] = THIS(Data);
    }
    virtual ~Data() {
    }
    std::string getName() const { return name; }
    std::string getTypeName() const { return typeName; }
    Stage getStage() const { return stage; }

    static PTR(Data) get(const std::string &name) {
      auto iterator(dataMap.find(name));
      return (iterator != dataMap.end()) ? iterator->second : nullptr;
    }
  protected:
    /**
     * \brief protected constructor for typed data.
     */
    Data(
      const std::string &name_, const std::string &typeName_
    ): name(name_), typeName(typeName_), stage(TYPED) {
      PTR(Data) mentionedData(dataMap[name_]);
      if (mentionedData) {
        if (mentionedData->getStage() != MENTIONED) {
          LOG(1,"Data") << "overwriting existing data: " << name_ << std::endl;
//          throw new EXCEPTION("Trying to overwrite existing data");
        }
      }
      dataMap[name_] = THIS(Data);
    }
    std::string name, typeName;
    Stage stage;

    static std::map<std::string, PTR(Data)> dataMap;
    static int64_t nextAnynomousDataId;
  };

  class TypedData: public Data {
  protected:
    /**
     * \brief Protected constructor for anonymous constant data.
     */
    TypedData(const std::string &typeName_): Data(nextName(), typeName_) {
    }
    /**
     * \brief Protected constructor for named data.
     */
    TypedData(
      const std::string &name_, const std::string &typeName_
    ): Data(name_, typeName_) {
    }

    static std::string nextName() {
      std::stringstream sStream;
      sStream << "Constant" << nextId++;
      return sStream.str();
    }

    /**
     * \brief next id number to be given anonymous constant data.
     * They will be named "Constant0", "Constant1", ...
     * regardless of the type.
     */
    static size_t nextId;
  };

  class TextData: public TypedData {
  public:
    TextData(const std::string &value_): TypedData("text"), value(value_) { }
    TextData(
      const std::string &name_, const std::string &value_
    ): TypedData(name_, "text"), value(value_) { }
    std::string value;
  };

  class BooleanData: public TypedData {
  public:
    BooleanData(const bool value_): TypedData("boolean"), value(value_) { }
    BooleanData(
      const std::string &name_, const bool value_
    ): TypedData(name_, "boolean"), value(value_) { }
    bool value;
  };

  class NumericData: public TypedData {
  protected:
    NumericData(const std::string &typeName_): TypedData(typeName_) { }
    NumericData(
      const std::string &name_, const std::string &typeName_
    ): TypedData(name_, typeName_) {
    }
  };

  class RealData: public NumericData {
  public:
    RealData(Real<64> value_): NumericData("real<64>"), value(value_) { }
    RealData(
      const std::string &name_, const Real<64> value_
    ): NumericData(name_, "real<64>"), value(value_) { }
    Real<64> value;
  };

  class IntegerData: public NumericData {
  public:
    IntegerData(int64_t value_): NumericData("integer"), value(value_) { }
    IntegerData(
      const std::string &name_, const int64_t value_
    ): NumericData(name_, "integer"), value(value_) { }
    int64_t value;
  };

  template <typename F=Real<>, typename TE=DefaultTensorEngine>
  class TensorData: public NumericData {
  public:
    TensorData(
      const PTR(ESC(Tensor<F,TE>)) &value_
    ): NumericData("tensor of " + TypeTraits<F>::getName()), value(value_) {
    }
    TensorData(
      const std::string &name_, const PTR(ESC(Tensor<F,TE>)) &value_
    ):
      NumericData(name_, "tensor of " + TypeTraits<F>::getName()), value(value_)
    {
    }
    PTR(ESC(Tensor<F,TE>)) value;
  };
}

#endif

