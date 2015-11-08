/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <string>
#include <ctf.hpp>

class Data {
  public:
    Data(std::string const &_name): name(_name) { }
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
      std::string const &name, std::string const &_value
    ): AtomData(name), value(_value) { }
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
      std::string const &name, double _value
    ): NumberData(name), value(_value) { }
    virtual std::string getTypeName() const { return "Real"; }
    double value;
};

class IntegerData: public NumberData {
  public:
    IntegerData(
      std::string const &name, int64_t _value
    ): NumberData(name), value(_value) { }
    virtual std::string getTypeName() const { return "Integer"; }
    int64_t value;
};

class TensorData: public AtomData {
  public:
    TensorData(
      std::string const &name, CTF::Tensor<> const &_value
    ): AtomData(name), value(_value) { }
    virtual std::string getTypeName() const { return "Tensor"; }
    CTF::Tensor<> value;
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
    Argument(std::string const &_name): name(_name) { }
    std::string getName() const { return name; }
    virtual bool isOutput() const = 0;
    virtual Data const *getData() const = 0;
  std::string name;
};

class InputArgument: public Argument {
  public:
    InputArgument(InputArgument const &i): Argument(i.name), data(i.data) {
    }
    InputArgument(
      std::string const &name, Data const *_data
    ): Argument(name), data(_data) { }
    virtual bool isOutput() const { return false; }
    virtual Data const *getData() const { return data; }
  Data const *data;
};

class OutputArgument: public Argument {
  public:
    OutputArgument(OutputArgument const &o): Argument(o.name), data(o.data) {
    }
    OutputArgument(
      std::string const &name, Data *_data
    ): Argument(name), data(_data) { }
    virtual bool isOutput() const { return true; }
    virtual Data const *getData() const { return data; }
    virtual Data *getData() { return data; }
  Data *data;
};


class Algorithm {
  public:
    Algorithm(std::vector<Argument const *> const &argumentList);
    virtual ~Algorithm();
//    virtual std::vector<std::string> getDefaultArgumentOrder() = 0;
//    virtual Data const *getDefaultInputData(std::string const &name) = 0;
//    virtual Data *getOutputData(std::string const &name) = 0;

    Data *getArgument(std::string const &name);
    std::string getTextArgument(std::string const &name);
    int64_t getIntegerArgument(std::string const &name);
    double getRealArgument(std::string const &name);
    CTF::Tensor<> const *getTensorArgument(std::string const &name);

  std::map<std::string, Data const *> inputs;
  std::map<std::string, Data *> outputs;
};



#endif
