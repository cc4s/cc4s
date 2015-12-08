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
       * \brief protected constructor for derived data types.
       */
      TypedData(
        std::string const &name_, std::string const &typeName_
      ): Data(name_, typeName_) {
      }
  };

  class TextData: public TypedData {
    public:
      TextData(
        std::string const &name_, std::string const &value_
      ): TypedData(name_, "text"), value(value_) { }
      std::string value;
  };

  class NumericData: public TypedData {
  protected:
    NumericData(
      std::string const &name_, std::string const &typeName_
    ): TypedData(name_, typeName_) {
    }
  };

  class RealData: public NumericData {
    public:
      RealData(
        std::string const &name_, double value_
      ): NumericData(name_, "real"), value(value_) { }
      double value;
  };

  class IntegerData: public NumericData {
    public:
      IntegerData(
        std::string const &name_, int64_t value_
      ): NumericData(name_, "integer"), value(value_) { }
      int64_t value;
  };

  template <typename F=double>
  class TensorData: public NumericData {
    public:
      TensorData(
        std::string const &name_, CTF::Tensor<F> *value_
      ): NumericData(name_, "tensor"), value(value_) {
      }
      CTF::Tensor<F> *value;
  };
}

#endif

