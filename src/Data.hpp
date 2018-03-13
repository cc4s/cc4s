/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DATA_DEFINED
#define DATA_DEFINED

#include <util/Log.hpp>
#include <math/Float.hpp>
#include <math/Complex.hpp>
#include <string>
#include <map>
#include <ctf.hpp>
// TODO: find out why Exception must be included after string,map and ctf
#include <util/Exception.hpp>

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
  class TypeTraits<Float64> {
  public:
    static std::string getName() { return "real"; }
  };
  template <>
  class TypeTraits<Complex64> {
  public:
    static std::string getName() { return "complex"; }
  };
  template <>
  class TypeTraits<Float128> {
  public:
    static std::string getName() { return "qreal"; }
  };
  template <>
  class TypeTraits<Complex128> {
  public:
    static std::string getName() { return "qcomplex"; }
  };

  class Data {
  public:
    enum Stage {
      MENTIONED = 0, TYPED = 1, ALLOCATED = 2, 
      READY = 3,
      UNUSED = 4, LINGERING = 5
    };
    Data(std::string const &name_);
    virtual ~Data() {
      dataMap[name] = nullptr;
    }
    std::string   getName() const { return name; }
    std::string getTypeName() const { return typeName; }
    Stage getStage() const { return stage; }

    static Data *get(std::string const &name) {
      auto iterator(dataMap.find(name));
      return (iterator != dataMap.end()) ? iterator->second : nullptr;
    }
  protected:
    /**
     * \brief protected constructor for typed data.
     */
    Data(
      std::string const &name_, std::string const &typeName_
    ): name(name_), typeName(typeName_), stage(TYPED) {
      Data *mentionedData(dataMap[name_]);
      if (mentionedData) {
        if (mentionedData->getStage() == MENTIONED) {
          delete mentionedData;
        } else {
          LOG(1,"Data") << "overwriting existing data: " << name_ << std::endl;
          delete mentionedData;
//          throw new EXCEPTION("Trying to overwrite existing data");
        }
      }
      dataMap[name_] = this;
    }
    std::string name, typeName;
    Stage stage;

    static std::map<std::string, Data *> dataMap;
    static int64_t nextAnynomousDataId;
  };

  class TypedData: public Data {
  protected:
    /**
     * \brief Protected constructor for anonymous constant data.
     */
    TypedData(std::string const &typeName_): Data(nextName(), typeName_) {
    }
    /**
     * \brief Protected constructor for named data.
     */
    TypedData(
      std::string const &name_, std::string const &typeName_
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
    static int nextId;
  };

  class TextData: public TypedData {
  public:
    TextData(std::string const &value_): TypedData("text"), value(value_) { }
    TextData(
      std::string const &name_, std::string const &value_
    ): TypedData(name_, "text"), value(value_) { }
    std::string value;
  };

  class BooleanData: public TypedData {
  public:
    BooleanData(bool const value_): TypedData("boolean"), value(value_) { }
    BooleanData(
      std::string const &name_, bool const value_
    ): TypedData(name_, "boolean"), value(value_) { }
    bool value;
  };

  class NumericData: public TypedData {
  protected:
    NumericData(std::string const &typeName_): TypedData(typeName_) { }
    NumericData(
      std::string const &name_, std::string const &typeName_
    ): TypedData(name_, typeName_) {
    }
  };

  class RealData: public NumericData {
  public:
    RealData(real value_): NumericData("real"), value(value_) { }
    RealData(
      std::string const &name_, const real value_
    ): NumericData(name_, "real"), value(value_) { }
    real value;
  };

  class IntegerData: public NumericData {
  public:
    IntegerData(int64_t value_): NumericData("integer"), value(value_) { }
    IntegerData(
      std::string const &name_, int64_t const value_
    ): NumericData(name_, "real"), value(value_) { }
    int64_t value;
  };

  template < typename F=double, typename T=CTF::Tensor<F> >
  class TensorData: public NumericData {
  public:
    TensorData(
      T *value_
    ): NumericData("tensor of " + TypeTraits<F>::getName()), value(value_) {
    }
    TensorData(
      std::string const &name_, T *value_
    ):
      NumericData(name_, "tensor of " + TypeTraits<F>::getName()), value(value_)
    {
    }
    virtual ~TensorData() {
      if (value) delete value;
    }
    T *value;
  };
}

#endif

