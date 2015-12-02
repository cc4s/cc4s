/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <string>
#include <ctf.hpp>

namespace cc4s {
  class Data {
    public:
      Data(std::string const &name_): name(name_) { }
      virtual std::string getTypeName() const = 0;
      std::string name;
  };

  class AtomData: public Data {
    public:
      AtomData(std::string const &name): Data(name) { }
  };

  class TextData: public AtomData {
    public:
      TextData(
        std::string const &name, std::string const &value_
      ): AtomData(name), value(value_) { }
      virtual std::string getTypeName() const { return "Text"; }
      std::string value;
  };

  class NumberData: public AtomData {
    public:
      NumberData(std::string const &name): AtomData(name) { }
  };

  class RealData: public NumberData {
    public:
      RealData(
        std::string const &name, double value_
      ): NumberData(name), value(value_) { }
      virtual std::string getTypeName() const { return "Real"; }
      double value;
  };

  class IntegerData: public NumberData {
    public:
      IntegerData(
        std::string const &name, int64_t value_
      ): NumberData(name), value(value_) { }
      virtual std::string getTypeName() const { return "Integer"; }
      int64_t value;
  };

  template <typename F=double>
  class TensorData: public AtomData {
    public:
      TensorData(
        std::string const &name, CTF::Tensor<F> const &value_
      ): AtomData(name), value(new CTF::Tensor<F>(value_)) { }
      TensorData(
        std::string const &name
      ): AtomData(name), value(nullptr) { }
      virtual std::string getTypeName() const {
        std::stringstream sStream("");
        sStream << "Tensor<" << typeid(F).name() << ">";
        return sStream.str();
      }
      CTF::Tensor<F> *value;
  };

  // TODO: maybe use iterator instead
  // TODO: if not, which form: have a look at the lazy lists project HL
  // NOTE: iterators are not stateless
  class SequenceData: public Data {
    public:
      virtual bool hasNext() = 0;
      virtual Data *next() = 0;
  };

  //TODO: add symbol to data
  //TODO: add type info to data

  class Argument {
    public:
      Argument(Argument const &a): name(a.name), data(a.data) {
      }
      Argument(
        std::string const &name_, Data *data_
      ): name(name_), data(data_) { }
      std::string getName() const { return name; }
      Data *getData() const { return data; }
    protected:
      std::string name;
      Data *data;
  };

  class Algorithm {
    public:
      Algorithm(std::vector<Argument const *> const &argumentList);
      virtual ~Algorithm();

      Data *getArgument(std::string const &name);
      std::string getTextArgument(std::string const &name);
      int64_t getIntegerArgument(std::string const &name);
      double getRealArgument(std::string const &name);
      template <typename F=double>
      TensorData<F> *getTensorDataArgument(std::string const &name);
      template <typename F=double>
      CTF::Tensor<F> *getTensorArgument(std::string const &name);
    protected:
      std::map<std::string, Data *> arguments;
  };
}


#endif
