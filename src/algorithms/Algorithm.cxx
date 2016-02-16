/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <math/Complex.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <sstream>

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
  if (!data) {
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
  if (!textData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Text, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return textData->value;
}
std::string Algorithm::getTextArgument(
  std::string const &name, std::string const &defaultValue
) {
  return isArgumentGiven(name) ? getTextArgument(name) : defaultValue;
}

int64_t Algorithm::getIntegerArgument(std::string const &name) {
  Data const *data(getArgumentData(name));
  IntegerData const *integerData = dynamic_cast<IntegerData const *>(data);
  if (!integerData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Integer, found " << data->getTypeName() << ".";
    throw new Exception(sstream.str());
  }
  return integerData->value;
}
int64_t Algorithm::getIntegerArgument(
  std::string const &name, int64_t const defaultValue
) {
  return isArgumentGiven(name) ? getIntegerArgument(name) : defaultValue;
}

double Algorithm::getRealArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return realData->value;
  IntegerData *integerData(dynamic_cast<IntegerData *>(data));
  if (integerData) return getRealArgumentFromInteger(integerData);
  return realData ? realData->value : getRealArgumentFromInteger(integerData);
  TensorData<double> *tensorData(dynamic_cast<TensorData<double> *>(data));
  if (tensorData) return getRealArgumentFromTensor(tensorData);
  std::stringstream sstream;
  sstream << "Incompatible type for argument: " << name << ". "
    << "Excpected Real, found " << data->getTypeName() << ".";
  throw new Exception(sstream.str());
}
double Algorithm::getRealArgument(
  std::string const &name, double const defaultValue
) {
  return isArgumentGiven(name) ? getRealArgument(name) : defaultValue;
}
double Algorithm::getRealArgumentFromInteger(IntegerData *integerData ) {
  double value(integerData->value);
  if (int64_t(value) != integerData->value) {
    LOG(0, "root") << "Warning: loss of precision in conversion from integer to real."
      << std::endl;
  }
  return value;
}
double Algorithm::getRealArgumentFromTensor(TensorData<double> *data) {
  Assert(
    data->value->order == 0,
    "Scalar expected in conversion from tensor to real."
  );
  int64_t index(0); double value;
  // retrieve the real value from the tensor
  data->value->read(1, &index, &value);
  return value;
}

template <typename F>
CTF::Tensor<F> *Algorithm::getTensorArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  TensorData<F> *tensorData(dynamic_cast<TensorData<F> *>(data));
  if (tensorData) return tensorData->value;
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return getTensorArgumentFromReal<F>(realData);
  // TODO: provide conversion routines from real to complex tensors
  std::stringstream sStream;
  sStream << "Incompatible type for argument: " << name << ". "
    << "Excpected tensor of " << TypeTraits<F>::getName()
    << ", found " << data->getTypeName() << ".";
  throw new Exception(sStream.str());
}
// instantiate
template
CTF::Tensor<double> *Algorithm::getTensorArgument(std::string const &);
template
CTF::Tensor<complex> *Algorithm::getTensorArgument(std::string const &);

template <typename F>
CTF::Tensor<F> *Algorithm::getTensorArgumentFromReal(RealData *realData) {
  // FIXME: left to leak memory...
  // a better solution would be to replace the RealData with the allocated
  // TensorData and support down-cast for Scalars to Real
  return new CTF::Scalar<F>(realData->value);
}
// instantiate
template
CTF::Tensor<double> *Algorithm::getTensorArgumentFromReal(RealData *);
template
CTF::Tensor<complex> *Algorithm::getTensorArgumentFromReal(RealData *);


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

