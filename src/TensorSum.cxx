#include <TensorSum.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(TensorSum);

TensorSum::TensorSum(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorSum::~TensorSum() {
}

/**
 * \brief Testing environement
 */
void TensorSum::run() {

  Tensor<> *A(getTensorArgument("A"));
  Tensor<> *B(getTensorArgument("B"));
  Tensor<> *C(getTensorArgument("Result"));
  (*C)[ getTextArgument("Resultindex").c_str() ] =
    getRealArgument("Afactor") * (*A)[ getTextArgument("Aindex").c_str() ] +
    getRealArgument("Bfactor") * (*B)[ getTextArgument("Bindex").c_str() ];

}