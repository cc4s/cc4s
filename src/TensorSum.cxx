#include <TensorSum.hpp>
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
  (*C)[ getTextArgument("ResultIndex").c_str() ] =
    getRealArgument("AFactor") * (*A)[ getTextArgument("AIndex").c_str() ] +
    getRealArgument("BFactor") * (*B)[ getTextArgument("BIndex").c_str() ];

}
