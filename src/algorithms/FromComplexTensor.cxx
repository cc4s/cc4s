#include <algorithms/FromComplexTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(FromComplexTensor);

FromComplexTensor::FromComplexTensor(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

FromComplexTensor::~FromComplexTensor() {
}

/**
 * \brief Testing environement
 */
void FromComplexTensor::run() {
  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<> *RealA(new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "RealA"));
  if (isArgumentGiven("imagA")) {
    Tensor<> *ImagA(new Tensor<>(A->order, A->lens, A->sym, *A->wrld, "ImagA"));
    fromComplexTensor(*A, *RealA, *ImagA);
    allocatedTensorArgument("imagA", ImagA);
  } else {
    fromComplexTensor(*A, *RealA);
  }
  allocatedTensorArgument("RealA", RealA);
}

