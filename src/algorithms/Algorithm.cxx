/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <math/Complex.hpp>
#include <tcc/DryTensor.hpp>
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

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
void Algorithm::dryRun() {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
}

bool Algorithm::isArgumentGiven(std::string const &name) {
  return arguments.find(name) != arguments.end();
}

Data *Algorithm::getArgumentData(std::string const &name) {
  auto dataIterator(arguments.find(name));
  if (dataIterator == arguments.end()) {
    std::stringstream sStream;
    sStream << "Missing argument: " << name;
//    throw new EXCEPTION(std::stringstream() << "Missing argument: " << name);
    throw new EXCEPTION(sStream.str());
  }
  Data *data = Data::get(dataIterator->second);
  if (!data) {
    std::stringstream sStream;
    sStream << "Missing data: " << dataIterator->second;
//    throw new EXCEPTION(std::stringstream() << "Missing data: " << dataIterator->second);
    throw new EXCEPTION(sStream.str());
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
    throw new EXCEPTION(sstream.str());
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
    throw new EXCEPTION(sstream.str());
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
  throw new EXCEPTION(sstream.str());
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

template <typename F, typename T>
T *Algorithm::getTensorArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  TensorData<F, T> *tensorData(dynamic_cast<TensorData<F, T> *>(data));
  if (tensorData) return tensorData->value;
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return getTensorArgumentFromReal<F, T>(realData);
  // TODO: provide conversion routines from real to complex tensors
  std::stringstream sStream;
  sStream << "Incompatible type for argument: " << name << ". "
    << "Excpected tensor of " << TypeTraits<F>::getName()
    << ", found " << data->getTypeName() << ".";
  throw new EXCEPTION(sStream.str());
}
// instantiate
template
CTF::Tensor<double> *Algorithm::getTensorArgument<
  double, CTF::Tensor<double>
>(std::string const &);
template
CTF::Tensor<complex> *Algorithm::getTensorArgument<
  complex, CTF::Tensor<complex>
>(std::string const &);
template
DryTensor<double> *Algorithm::getTensorArgument<
  double, DryTensor<double>
>(std::string const &);
template
DryTensor<complex> *Algorithm::getTensorArgument<
  complex, DryTensor<complex>
>(std::string const &);


/**
 * \brief Traits for retrieving the Scalar, Vector and Matrix tensor type.
 */
template < typename F, typename T=CTF::Tensor<F> >
class TensorTypeTraits;

template <typename F>
class TensorTypeTraits< F, CTF::Tensor<F> > {
public:
  typedef CTF::Tensor<F> BaseType;
  typedef CTF::Scalar<F> ScalarType;
  typedef CTF::Vector<F> VectorType;
  typedef CTF::Matrix<F> MatrixType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Matrix<F> > {
public:
  typedef CTF::Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Vector<F> > {
public:
  typedef CTF::Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, CTF::Scalar<F> > {
public:
  typedef CTF::Tensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryTensor<F> > {
public:
  typedef DryTensor<F> BaseType;
  typedef DryScalar<F> ScalarType;
  typedef DryVector<F> VectorType;
  typedef DryMatrix<F> MatrixType;
};
template <typename F>
class TensorTypeTraits< F, DryMatrix<F> > {
public:
  typedef DryTensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryVector<F> > {
public:
  typedef DryTensor<F> BaseType;
};
template <typename F>
class TensorTypeTraits< F, DryScalar<F> > {
public:
  typedef DryTensor<F> BaseType;
};


/**
 * \brief Converts the given real data into a scalar tensor.
 */
template <typename F, typename T>
T *Algorithm::getTensorArgumentFromReal(RealData *realData) {
  // FIXME: left to leak memory...
  // a better solution would be to replace the RealData with the allocated
  // TensorData and support down-cast for Scalars to Real
  return new typename TensorTypeTraits<F,T>::ScalarType(realData->value);
}
// instantiate
template
CTF::Tensor<double> *Algorithm::getTensorArgumentFromReal<
  double, CTF::Tensor<double>
>(RealData *);
template
CTF::Tensor<complex> *Algorithm::getTensorArgumentFromReal<
  complex, CTF::Tensor<complex>
>(RealData *);
template
DryTensor<double> *Algorithm::getTensorArgumentFromReal<
  double, DryTensor<double>
>(RealData *);
template
DryTensor<complex> *Algorithm::getTensorArgumentFromReal<
  complex, DryTensor<complex>
>(RealData *);

template <typename F, typename T>
void Algorithm::allocatedTensorArgument(
  std::string const &name, T *tensor
) {
  Data *mentionedData(getArgumentData(name));
  new TensorData<F, typename TensorTypeTraits<F, T>::BaseType>(
    mentionedData->getName(), tensor
  );
  // NOTE: the constructor of TensorData enteres its location in the
  // data map and destroys the previous content, i.e. mentionedData.
}
// instantiate
template
void Algorithm::allocatedTensorArgument<
  double, CTF::Tensor<double>
>(std::string const &name, CTF::Tensor<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, CTF::Matrix<double>
>(std::string const &name, CTF::Matrix<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, CTF::Vector<double>
>(std::string const &name, CTF::Vector<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, CTF::Scalar<double>
>(std::string const &name, CTF::Scalar<double> *tensor);

template
void Algorithm::allocatedTensorArgument<
  complex, CTF::Tensor<complex>
>(std::string const &name, CTF::Tensor<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, CTF::Matrix<complex>
>(std::string const &name, CTF::Matrix<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, CTF::Vector<complex>
>(std::string const &name, CTF::Vector<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, CTF::Scalar<complex>
>(std::string const &name, CTF::Scalar<complex> *tensor);

template
void Algorithm::allocatedTensorArgument<
  double, DryTensor<double>
>(std::string const &name, DryTensor<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, DryMatrix<double>
>(std::string const &name, DryMatrix<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, DryVector<double>
>(std::string const &name, DryVector<double> *tensor);
template
void Algorithm::allocatedTensorArgument<
  double, DryScalar<double>
>(std::string const &name, DryScalar<double> *tensor);

template
void Algorithm::allocatedTensorArgument<
  complex, DryTensor<complex>
>(std::string const &name, DryTensor<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, DryMatrix<complex>
>(std::string const &name, DryMatrix<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, DryVector<complex>
>(std::string const &name, DryVector<complex> *tensor);
template
void Algorithm::allocatedTensorArgument<
  complex, DryScalar<complex>
>(std::string const &name, DryScalar<complex> *tensor);

void Algorithm::setRealArgument(std::string const &name, double const value) {
  Data *mentionedData(getArgumentData(name));
  new RealData(mentionedData->getName(), value);  
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;

