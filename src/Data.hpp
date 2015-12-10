/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DATA_DEFINED
#define DATA_DEFINED

#include <string>
#include <map>
#include <ctf.hpp>
#include <util/Exception.hpp>

namespace cc4s {
  class Data {
  public:
    enum Stage {
      MENTIONED = 0, TYPED = 1, ALLOCATED = 2, 
      READY = 3,
      UNUSED = 4, LINGERING = 5
    };
    Data(
      std::string const &name_
    ): name(name_), typeName("unknown"), stage(MENTIONED) {
      dataMap[name_] = this;
    }
    virtual ~Data() { }
    std::string getName() const { return name; }
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
      if (mentionedData != nullptr) {
        if (mentionedData->getStage() == MENTIONED) {
          delete mentionedData;
        } else {
          throw new Exception("Trying to overwrite existing data");
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
    TypedData(std::string const &typeName_): Data("", typeName_) {
      std::stringstream sStream;
      sStream << "Constant" << nextId++;
      name = sStream.str();
    }
    /**
     * \brief Protected constructor for named data.
     */
    TypedData(
      std::string const &name_, std::string const &typeName_
    ): Data(name_, typeName_) {
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
    RealData(double value_): NumericData("real"), value(value_) { }
    RealData(
      std::string const &name_, double const value_
    ): NumericData(name_, "real"), value(value_) { }
    double value;
  };

  class IntegerData: public NumericData {
  public:
    IntegerData(int64_t value_): NumericData("integer"), value(value_) { }
    IntegerData(
      std::string const &name_, int64_t const value_
    ): NumericData(name_, "real"), value(value_) { }
    int64_t value;
  };

  template <typename F=double>
  class TensorData: public NumericData {
  public:
    TensorData(CTF::Tensor<F> *value_): NumericData("tensor"), value(value_) {
    }
    TensorData(
      std::string const &name_, CTF::Tensor<F> *value_
    ): NumericData(name_, "tensor"), value(value_) {
    }
    CTF::Tensor<F> *value;
  };
}

#endif

