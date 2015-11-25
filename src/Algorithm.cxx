/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Algorithm.hpp"
#include "Exception.hpp"

#include <iostream>

using namespace cc4s;

Algorithm::Algorithm(std::vector<Argument const *> const &arguments) {
  for (auto arg(arguments.begin()); arg != arguments.end(); ++arg) {
    Argument const *argument = *arg;
    if (argument->isOutput()) {
      outputs[argument->getName()] =
        dynamic_cast<OutputArgument *>(
          const_cast<Argument *>(argument)
        )->getData();
    } else {
      inputs[argument->getName()] = argument->getData();
    }
  }
}

Algorithm::~Algorithm() {
}

Data *Algorithm::getArgument(std::string const &name) {
  auto dataIterator = inputs.find(name);
  if (dataIterator == inputs.end()) {
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

CTF::Tensor<> const *Algorithm::getTensorArgument(std::string const &name) {
  Data const *data = getArgument(name);
  TensorData const *tensorData = dynamic_cast<TensorData const *>(data);
  if (tensorData == nullptr) {
    std::stringstream sstream;
    sstream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Tensor, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return &tensorData->value;
}
