#include <algorithms/ComplexTensorSum.hpp>
#include <math/Complex.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorSum);

ComplexTensorSum::ComplexTensorSum(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ComplexTensorSum::~ComplexTensorSum() {
}

/**
 * \brief Testing environement
 */
void ComplexTensorSum::run() {

  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<complex> *B(getTensorArgument<complex>("B"));
  Tensor<complex> *C(getTensorArgument<complex>("Result"));
  (*C)[ getTextArgument("ResultIndex").c_str() ] =
    getRealArgument("AFactor") * (*A)[ getTextArgument("AIndex").c_str() ] +
    getRealArgument("BFactor") * (*B)[ getTextArgument("BIndex").c_str() ];

}
