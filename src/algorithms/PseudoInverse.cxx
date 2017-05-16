/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/PseudoInverse.hpp>
#include <math/PseudoInverseHermitianSvd.hpp>
#include <util/Log.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(PseudoInverse);

PseudoInverse::PseudoInverse(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PseudoInverse::~PseudoInverse() {
}

void PseudoInverse::run() {
  Tensor<complex> *A(getTensorArgument<complex>("A"));
  // FIXME: use cast operators provided by CTF as soon as supported
  if (A->order != 2) throw new EXCEPTION("Matrix expected as argument A");
  Matrix<complex> *MatrixA(static_cast<Matrix<complex> *>(A));
  PseudoInverseHermitianSvd<complex> pseudoInverse(*MatrixA);
  Tensor<complex> *inverseA(new Tensor<complex>(pseudoInverse.get()));
  allocatedTensorArgument<complex>("InverseA", inverseA);
}

