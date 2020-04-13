/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

#include <iostream>
#include <sstream>

using namespace cc4s;

Algorithm::Algorithm(const std::vector<Argument> &argumentList) {
  for (auto argument: argumentList) {
    arguments[argument.getName()] = argument.getData();
  }
}

Algorithm::~Algorithm() {
}

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
void Algorithm::dryRun() {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
}

bool Algorithm::isArgumentGiven(const std::string &name) {
  return arguments.find(name) != arguments.end();
}

Ptr<Data> Algorithm::getArgumentData(const std::string &name) {
  auto dataIterator(arguments.find(name));
  if (dataIterator == arguments.end()) {
    std::stringstream sStream;
    sStream << "Missing argument: " << name;
//    throw new EXCEPTION(std::stringstream() << "Missing argument: " << name);
    throw new EXCEPTION(sStream.str());
  }
  Ptr<Data> data = Data::get(dataIterator->second);
  if (!data) {
    std::stringstream sStream;
    sStream << "Missing data: " << dataIterator->second;
//    throw new EXCEPTION(std::stringstream() << "Missing data: " << dataIterator->second);
    throw new EXCEPTION(sStream.str());
  }
  return data;
}

const std::string &Algorithm::getTextArgument(const std::string &name) {
  Ptr<Data> data(getArgumentData(name));
  Ptr<TextData> textData( dynamic_pointer_cast<TextData>(data) );
  if (!textData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Text, found " << data->getTypeName() << ".";
    throw new EXCEPTION(sstream.str());
  }
  return textData->value;
}
const std::string &Algorithm::getTextArgument(
  const std::string &name, const std::string &defaultValue
) {
  return isArgumentGiven(name) ? getTextArgument(name) : defaultValue;
}

bool Algorithm::getBooleanArgument(const std::string &name) {
  /*
   *TODO: Do this without the getTextArgument function, because in this
   *case the parser want to have quotes on the boolean value i.e.
   *  (myflag "true")
   */
  std::string text(getTextArgument(name));
  if (
    text.compare(".TRUE.") == 0 ||
    text.compare("true") == 0 ||
    text.compare("True") == 0 ||
    text.compare("TRUE") == 0 ||
    text.compare("1") == 0 ||
    text.compare("t") == 0 ||
    text.compare("T") == 0
  ) {
    return true;
  } else {
    return false;
  }
}
bool Algorithm::getBooleanArgument(
  const std::string &name, const bool defaultValue
) {
  return isArgumentGiven(name) ? getBooleanArgument(name) : defaultValue;
}

int64_t Algorithm::getIntegerArgument(const std::string &name) {
  Ptr<Data> data(getArgumentData(name));
  auto integerData( dynamic_pointer_cast<IntegerData >(data) );
  if (!integerData) {
    std::stringstream sstream;
    sstream << "Incompatible type for argument: " << name << ". "
      << "Excpected Integer, found " << data->getTypeName() << ".";
    throw new EXCEPTION(sstream.str());
  }
  return integerData->value;
}
int64_t Algorithm::getIntegerArgument(
  const std::string &name, const int64_t defaultValue
) {
  return isArgumentGiven(name) ? getIntegerArgument(name) : defaultValue;
}

Real<> Algorithm::getRealArgument(const std::string &name) {
  Ptr<Data> data(getArgumentData(name));
  auto realData( dynamic_pointer_cast<RealData>(data) );
  if (realData) return realData->value;
  auto integerData( dynamic_pointer_cast<IntegerData>(data) );
  if (integerData) return getRealArgumentFromInteger(integerData);
  auto defaultTensorData(
    dynamic_pointer_cast<TensorData<Real<>,DefaultTensorEngine>>(data)
  );
  if (defaultTensorData) {
    return getRealArgumentFromTensor<Real<>,DefaultTensorEngine>(
      defaultTensorData
    );
  }
  auto dryTensorData(
    dynamic_pointer_cast<TensorData<Real<>,DryTensorEngine>>(data)
  );
  if (dryTensorData) {
    return getRealArgumentFromTensor<Real<>,DryTensorEngine>(
      dryTensorData
    );
  }
  std::stringstream sstream;
  sstream << "Incompatible type for argument: " << name << ". "
    << "Excpected Real, found " << data->getTypeName() << ".";
  throw new EXCEPTION(sstream.str());
}
Real<> Algorithm::getRealArgument(
  const std::string &name, const Real<> defaultValue
) {
  return isArgumentGiven(name) ? getRealArgument(name) : defaultValue;
}
Real<> Algorithm::getRealArgumentFromInteger(
  const Ptr<IntegerData> &integerData
) {
  Real<> value(integerData->value);
  if (int64_t(value) != integerData->value) {
    LOG(0, "root") << "Warning: loss of precision in conversion from integer to real."
      << std::endl;
  }
  return value;
}
template <typename F, typename TE>
F Algorithm::getRealArgumentFromTensor(
  const Ptr<const TensorData<F,TE>> &tensorData
) {
  Assert(
    tensorData->value->lens.size() == 0,
    "Scalar expected in conversion from tensor to real."
  );
  // FIXME:
  // read the real value from the tensor
/*
  CTF::Scalar<Real<>> scalar;
  scalar[""] = (*data->value)[""];
  return scalar.get_val();
*/
  return F(0);
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> Algorithm::getTensorArgument(
  const std::string &name
) {
  auto data(getArgumentData(name));
  auto tensorData( dynamic_pointer_cast<TensorData<F,TE>>(data) );
  if (tensorData) return tensorData->value;
  auto realData( dynamic_pointer_cast<RealData>(data) );
  if (realData) return getTensorArgumentFromReal<F,TE>(realData);
  std::stringstream sStream;
  sStream << "Incompatible type for argument: " << name << ". "
    << "Excpected tensor of " << TypeTraits<F>::getName()
    << ", found " << data->getTypeName() << ".";
  throw new EXCEPTION(sStream.str());
}
// instantiate
template
Ptr<Tensor<Real<64>,DefaultTensorEngine>>
Algorithm::getTensorArgument<Real<64>, DefaultTensorEngine>(
  const std::string &
);
template
Ptr<Tensor<Complex<64>,DefaultTensorEngine>>
Algorithm::getTensorArgument<Complex<64>, DefaultTensorEngine>(
  const std::string &
);
template
Ptr<Tensor<Real<64>,DryTensorEngine>>
Algorithm::getTensorArgument<Real<64>, DryTensorEngine>(
  const std::string &
);
template
Ptr<Tensor<Complex<64>,DryTensorEngine>>
Algorithm::getTensorArgument<Complex<64>, DryTensorEngine>(
  const std::string &
);

// TODO: 128 bit tensors



/**
 * \brief Converts the given real data into a scalar tensor.
 */
template <typename F, typename TE>
Ptr<Tensor<F,TE>> Algorithm::getTensorArgumentFromReal(
  const Ptr<RealData> &realData
) {
  // FIXME: write scalar value to tensor
  return Tcc<TE>::template tensor<F>(
    std::vector<size_t>(), realData->getName()
  );
}
// instantiate
template
Ptr<Tensor<Real<64>,DefaultTensorEngine>>
Algorithm::getTensorArgumentFromReal<Real<64>,DefaultTensorEngine>(
  const Ptr<RealData> &
);
template
Ptr<Tensor<Complex<64>,DefaultTensorEngine>>
Algorithm::getTensorArgumentFromReal<Complex<64>,DefaultTensorEngine>(
  const Ptr<RealData> &
);
template
Ptr<Tensor<Real<64>,DryTensorEngine>>
Algorithm::getTensorArgumentFromReal<Real<64>,DryTensorEngine>(
  const Ptr<RealData> &
);
template
Ptr<Tensor<Complex<64>,DryTensorEngine>>
Algorithm::getTensorArgumentFromReal<Complex<64>,DryTensorEngine>(
  const Ptr<RealData> &
);



template <typename F, typename TE>
void Algorithm::setTensorArgument(
  const std::string &name, const Ptr<Tensor<F,TE>> &tensor
) {
  Ptr<Data> mentionedData(getArgumentData(name));
  New<TensorData<F,TE>>(mentionedData->getName(), tensor);
  // NOTE: the constructor of TensorData enteres its location in the
  // data map and destroys the previous content, i.e. mentionedData.
}
// instantiate
template
void Algorithm::setTensorArgument<Real<64>, DefaultTensorEngine>(
  const std::string &name,
  const Ptr<Tensor<Real<64>,DefaultTensorEngine>> &tensor
);
template
void Algorithm::setTensorArgument<Complex<64>, DefaultTensorEngine>(
  const std::string &name,
  const Ptr<Tensor<Complex<64>,DefaultTensorEngine>> &tensor
);
template
void Algorithm::setTensorArgument<Real<64>, DryTensorEngine>(
  const std::string &name,
  const Ptr<Tensor<Real<64>,DryTensorEngine>> &tensor
);
template
void Algorithm::setTensorArgument<Complex<64>, DryTensorEngine>(
  const std::string &name,
  const Ptr<Tensor<Complex<64>,DryTensorEngine>> &tensor
);

// TODO: 128 bit tensors


void Algorithm::setRealArgument(const std::string &name, const Real<> value) {
  Ptr<Data> mentionedData(getArgumentData(name));
  New<RealData>(mentionedData->getName(), value);
}

void Algorithm::setIntegerArgument(
  const std::string &name, const int64_t value
) {
  Ptr<Data> mentionedData(getArgumentData(name));
  New<IntegerData>(mentionedData->getName(), value);
}

Ptr<AlgorithmFactory::AlgorithmMap> AlgorithmFactory::algorithmMap;

