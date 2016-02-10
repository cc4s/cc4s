#include <ComplexTensorNorm.hpp>
#include <util/Complex.hpp>
#include <util/MathFunctions.hpp>
#include <util/Log.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorNorm);

ComplexTensorNorm::ComplexTensorNorm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ComplexTensorNorm::~ComplexTensorNorm() {
}

/**
 * \brief Testing environement
 */
void ComplexTensorNorm::run() {
  Tensor<complex> *A(getTensorArgument<complex>("A"));
  double norm(frobeniusNorm(*A));
  LOG(0) << "|A| = " << norm << std::endl;
  setRealArgument("Norm", norm);
}
