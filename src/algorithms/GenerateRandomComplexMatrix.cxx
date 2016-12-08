#include <algorithms/GenerateRandomComplexMatrix.hpp>
#include <math/RandomTensor.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(GenerateRandomComplexMatrix);

GenerateRandomComplexMatrix::GenerateRandomComplexMatrix(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

GenerateRandomComplexMatrix::~GenerateRandomComplexMatrix() {
}

/**
 * \brief Testing environement
 */
void GenerateRandomComplexMatrix::run() {
  int m(getIntegerArgument("m", 200)), n(getIntegerArgument("n", 200));
  std::string symmetry(getTextArgument("symmetric", "none"));
  int sym(NS);
  if (symmetry == "hermitian") {
    throw new EXCEPTION("Hermitian symmetry of complex tensors not yet supported.");
  }
  Matrix<complex> *C(new Matrix<complex>(m, n, sym, *Cc4s::world, "C"));
  setRandomTensor(*C);
  allocatedTensorArgument<complex>("Result", C);
}

