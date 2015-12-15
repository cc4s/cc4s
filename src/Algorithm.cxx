/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Algorithm.hpp>
#include <util/Complex.hpp>
#include <util/Exception.hpp>

#include <iostream>

using namespace cc4s;

Algorithm::Algorithm(std::vector<Argument> const &argumentList) {
  for (auto arg(argumentList.begin()); arg != argumentList.end(); ++arg) {
    Argument argument = *arg;
    arguments[argument.getName()] = argument.getData();
  }
}

Algorithm::~Algorithm() {
}

bool Algorithm::isArgumentGiven(std::string const &name) {
  return arguments.find(name) != arguments.end();
}

Data *Algorithm::getArgumentData(std::string const &name) {
  auto dataIterator(arguments.find(name));
  if (dataIterator == arguments.end()) {
    std::stringstream sStream;
    sStream << "Missing argument: " << name;
//    throw new Exception(std::stringstream() << "Missing argument: " << name);
    throw new Exception(sStream.str());
  }
  Data *data = Data::get(dataIterator->second);
  if (data == nullptr) {
    std::stringstream sStream;
    sStream << "Missing data: " << dataIterator->second;
//    throw new Exception(std::stringstream() << "Missing data: " << dataIterator->second);
    throw new Exception(sStream.str());
  }
  return data;
}

std::string Algorithm::getTextArgument(std::string const &name) {
  Data *data(getArgumentData(name));
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
  Data const *data(getArgumentData(name));
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
  Data const *data(getArgumentData(name));
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
CTF::Tensor<F> *Algorithm::getTensorArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  TensorData<F> *tensorData = dynamic_cast<TensorData<F> *>(data);
  if (tensorData == nullptr) {
    std::stringstream sStream;
    sStream << "Incompatible tpye for argument: " << name << ". "
      << "Excpected Tensor, found " << data->getTypeName() << ".";
    throw new Exception(sStream.str());
  }
  return tensorData->value;
}
// instantiate
template
CTF::Tensor<double> *Algorithm::getTensorArgument(std::string const &name);
template
CTF::Tensor<complex> *Algorithm::getTensorArgument(std::string const &name);

template <typename F>
void Algorithm::allocatedTensorArgument(
  std::string const &name, CTF::Tensor<F> *tensor
) {
  Data *mentionedData(getArgumentData(name));
  new TensorData<F>(mentionedData->getName(), tensor);
  // NOTE: the constructor of TensorData enteres its location in the
  // data map and destroys the previous content, i.e. mentionedData.
}
// instantiate
template
void Algorithm::allocatedTensorArgument(
  std::string const &name, CTF::Tensor<double> *tensor
);
template
void Algorithm::allocatedTensorArgument(
  std::string const &name, CTF::Tensor<complex> *tensor
);

void Algorithm::setRealArgument(std::string const &name, double const value) {
  Data *mentionedData(getArgumentData(name));
  new RealData(mentionedData->getName(), value);  
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;

