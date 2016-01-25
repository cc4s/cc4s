#include <TensorContraction.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(TensorContraction);

TensorContraction::TensorContraction(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorContraction::~TensorContraction() {
}

/**
 * \brief Testing environement
 */
void TensorContraction::run() {

  Tensor<> *A(getTensorArgument("A"));
  Tensor<> *B(getTensorArgument("B"));
  Tensor<> *C(getTensorArgument("Result"));
  (*C)[ getTextArgument("Resultindex").c_str() ] =
    (*A)[ getTextArgument("Aindex").c_str() ] *
    (*B)[ getTextArgument("Bindex").c_str() ];

}