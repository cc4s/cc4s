#include <algorithms/GenerateRandomMatrix.hpp>
#include <math/Complex.hpp>
#include <math/RandomTensor.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(GenerateRandomMatrix);

GenerateRandomMatrix::GenerateRandomMatrix(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

GenerateRandomMatrix::~GenerateRandomMatrix() {
}

/**
 * \brief Testing environement
 */
void GenerateRandomMatrix::run() {
  int m(getIntegerArgument("m", 200)), n(getIntegerArgument("n", 200));
  std::string symmetry(getTextArgument("symmetric", "none"));
  int sym(NS);
  if (symmetry == "symmetric" || symmetry == "hermitian") {
    sym = SY;
  } else if (symmetry == "anti") {
    sym = AS;
  } else if (symmetry == "hollow") {
    sym = SH;
  }
  Matrix<> *C(new Matrix<>(m, n, sym, *Cc4s::world, "C"));
  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*C, normalDistribution, random);
  allocatedTensorArgument("Result", C);
}

