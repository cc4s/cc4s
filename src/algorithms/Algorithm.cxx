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

bool Algorithm::getBooleanArgument(std::string const &name) {
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
  std::string const &name, bool const &defaultValue
) {
  return isArgumentGiven(name) ? getBooleanArgument(name) : defaultValue;
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

real Algorithm::getRealArgument(std::string const &name) {
  Data *data(getArgumentData(name));
  RealData *realData(dynamic_cast<RealData *>(data));
  if (realData) return realData->value;
  IntegerData *integerData(dynamic_cast<IntegerData *>(data));
  if (integerData) return getRealArgumentFromInteger(integerData);
  TensorData<real> *tensorData(dynamic_cast<TensorData<real> *>(data));
  if (tensorData) return getRealArgumentFromTensor(tensorData);
  std::stringstream sstream;
  sstream << "Incompatible type for argument: " << name << ". "
    << "Excpected Real, found " << data->getTypeName() << ".";
  throw new EXCEPTION(sstream.str());
}
real Algorithm::getRealArgument(
  const std::string &name, const real defaultValue
) {
  return isArgumentGiven(name) ? getRealArgument(name) : defaultValue;
}
real Algorithm::getRealArgumentFromInteger(IntegerData *integerData ) {
  real value(integerData->value);
  if (int64_t(value) != integerData->value) {
    LOG(0, "root") << "Warning: loss of precision in conversion from integer to real."
      << std::endl;
  }
  return value;
}
real Algorithm::getRealArgumentFromTensor(TensorData<real> *data) {
  Assert(
    data->value->order == 0,
    "Scalar expected in conversion from tensor to real."
  );
  int64_t index(0); real value;
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
CTF::Tensor<Float64> *Algorithm::getTensorArgument<
  Float64, CTF::Tensor<Float64>
>(std::string const &);
template
CTF::Tensor<Complex64> *Algorithm::getTensorArgument<
  Complex64, CTF::Tensor<Complex64>
>(std::string const &);
template
DryTensor<Float64> *Algorithm::getTensorArgument<
  Float64, DryTensor<Float64>
>(std::string const &);
template
DryTensor<Complex64> *Algorithm::getTensorArgument<
  Complex64, DryTensor<Complex64>
>(std::string const &);

#ifndef INTEL_COMPILER
template
CTF::Tensor<Float128> *Algorithm::getTensorArgument<
  Float128, CTF::Tensor<Float128>
>(std::string const &);
template
CTF::Tensor<Complex128> *Algorithm::getTensorArgument<
  Complex128, CTF::Tensor<Complex128>
>(std::string const &);
template
DryTensor<Float128> *Algorithm::getTensorArgument<
  Float128, DryTensor<Float128>
>(std::string const &);
template
DryTensor<Complex128> *Algorithm::getTensorArgument<
  Complex128, DryTensor<Complex128>
>(std::string const &);
#endif


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
CTF::Tensor<Float64> *Algorithm::getTensorArgumentFromReal<
  Float64, CTF::Tensor<Float64>
>(RealData *);
template
CTF::Tensor<Complex64> *Algorithm::getTensorArgumentFromReal<
  Complex64, CTF::Tensor<Complex64>
>(RealData *);
template
DryTensor<Float64> *Algorithm::getTensorArgumentFromReal<
  Float64, DryTensor<Float64>
>(RealData *);
template
DryTensor<Complex64> *Algorithm::getTensorArgumentFromReal<
  Complex64, DryTensor<Complex64>
>(RealData *);

#ifndef INTEL_COMPILER
template
CTF::Tensor<Float128> *Algorithm::getTensorArgumentFromReal<
  Float128, CTF::Tensor<Float128>
>(RealData *);
template
CTF::Tensor<Complex128> *Algorithm::getTensorArgumentFromReal<
  Complex128, CTF::Tensor<Complex128>
>(RealData *);
template
DryTensor<Float128> *Algorithm::getTensorArgumentFromReal<
  Float128, DryTensor<Float128>
>(RealData *);
template
DryTensor<Complex128> *Algorithm::getTensorArgumentFromReal<
  Complex128, DryTensor<Complex128>
>(RealData *);
#endif

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
  Float64, CTF::Tensor<Float64>
>(std::string const &name, CTF::Tensor<Float64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  real, CTF::Matrix<real>
>(std::string const &name, CTF::Matrix<real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  real, CTF::Vector<real>
>(std::string const &name, CTF::Vector<real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  real, CTF::Scalar<real>
>(std::string const &name, CTF::Scalar<real> *tensor);

template
void Algorithm::allocatedTensorArgument<
  Complex64, CTF::Tensor<Complex64>
>(std::string const &name, CTF::Tensor<Complex64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
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
  Float64, DryTensor<Float64>
>(std::string const &name, DryTensor<Float64> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
template
void Algorithm::allocatedTensorArgument<
  real, DryMatrix<real>
>(std::string const &name, DryMatrix<real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  real, DryVector<real>
>(std::string const &name, DryVector<real> *tensor);
template
void Algorithm::allocatedTensorArgument<
  real, DryScalar<real>
>(std::string const &name, DryScalar<real> *tensor);

template
void Algorithm::allocatedTensorArgument<
  complex, DryTensor<Complex64>
>(std::string const &name, DryTensor<complex> *tensor);
// TODO: remove specialized tensors (matrix, vector, scalar)
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

#ifndef INTEL_COMPILER
template
void Algorithm::allocatedTensorArgument<
  Float128, CTF::Tensor<Float128>
>(std::string const &name, CTF::Tensor<Float128> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Complex128, CTF::Tensor<Complex128>
>(std::string const &name, CTF::Tensor<Complex128> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Float128, DryTensor<Float128>
>(std::string const &name, DryTensor<Float128> *tensor);
template
void Algorithm::allocatedTensorArgument<
  Complex128, DryTensor<Complex128>
>(std::string const &name, DryTensor<Complex128> *tensor);
#endif


void Algorithm::setRealArgument(std::string const &name, const real value) {
  Data *mentionedData(getArgumentData(name));
  new RealData(mentionedData->getName(), value);
}

void Algorithm::setIntegerArgument(std::string const &name, int const value) {
  Data *mentionedData(getArgumentData(name));
  new IntegerData(mentionedData->getName(), value);
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;

