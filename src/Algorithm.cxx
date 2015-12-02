/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Algorithm.hpp>
#include <util/Complex.hpp>
#include <util/Exception.hpp>

#include <iostream>

using namespace cc4s;

Algorithm::Algorithm(std::vector<Argument const *> const &argumentList) {
  for (auto arg(argumentList.begin()); arg != argumentList.end(); ++arg) {
    Argument const *argument = *arg;
    arguments[argument->getName()] = argument->getData();
  }
}

Algorithm::~Algorithm() {
}

Data *Algorithm::getArgument(std::string const &name) {
  auto dataIterator = arguments.find(name);
  if (dataIterator == arguments.end()) {
    std::stringstream sstream;
    sstream << "Missing argument: " << name;
    throw new Exception(sstream.str());
  }
  //FIXME: decide whether to use explicit in/out information
  return const_cast<Data *>(dataIterator->second);
}

std::string Algorithm::getTextArgument(std::string const &name) {
  Data *data = getArgument(name);
  TextData const *textData = dynamic_cast<TextData const *>(data);
  if (textData == nullptr) {
    std::stringstream sstream;
    sstream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Text, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return textData->value;
}

int64_t Algorithm::getIntegerArgument(std::string const &name) {
  Data const *data = getArgument(name);
  IntegerData const *integerData = dynamic_cast<IntegerData const *>(data);
  if (integerData == nullptr) {
    std::stringstream sstream;
    sstream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Integer, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return integerData->value;
}

double Algorithm::getRealArgument(std::string const &name) {
  Data const *data = getArgument(name);
  RealData const *realData = dynamic_cast<RealData const *>(data);
  if (realData == nullptr) {
    std::stringstream sstream;
    sstream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Real, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return realData->value;
}

template <typename F>
TensorData<F> *Algorithm::getTensorDataArgument(std::string const &name) {
  Data *data = getArgument(name);
  TensorData<F> *tensorData = dynamic_cast<TensorData<F> *>(data);
  if (tensorData == nullptr) {
    std::stringstream sStream;
    sStream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Tensor, found " << data->getTypeName() << ".";
    throw new Exception(sStream.str());
  }
  return tensorData;
}

template <typename F>
CTF::Tensor<F> *Algorithm::getTensorArgument(std::string const &name) {
  return getTensorDataArgument<F>(name)->value;
}


// instantiate
template
TensorData<double> *Algorithm::getTensorDataArgument(std::string const &name);
template
TensorData<complex> *Algorithm::getTensorDataArgument(std::string const &name);

template
CTF::Tensor<double> *Algorithm::getTensorArgument(std::string const &name);
template
CTF::Tensor<complex> *Algorithm::getTensorArgument(std::string const &name);

