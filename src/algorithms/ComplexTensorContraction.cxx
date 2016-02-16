#include <algorithms/ComplexTensorContraction.hpp>
#include <math/Complex.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorContraction);

ComplexTensorContraction::ComplexTensorContraction(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ComplexTensorContraction::~ComplexTensorContraction() {
}

/**
 * \brief Testing environement
 */
void ComplexTensorContraction::run() {

  Tensor<complex> *A(getTensorArgument<complex>("A"));
  Tensor<complex> *B(getTensorArgument<complex>("B"));
  Tensor<complex> *C(getTensorArgument<complex>("Result"));
  (*C)[ getTextArgument("ResultIndex").c_str() ] =
    (*A)[ getTextArgument("AIndex").c_str() ] *
    (*B)[ getTextArgument("BIndex").c_str() ];

}
