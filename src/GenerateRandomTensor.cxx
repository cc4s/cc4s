#include <GenerateRandomTensor.hpp>
#include <util/RandomTensor.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(GenerateRandomTensor);

GenerateRandomTensor::GenerateRandomTensor(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

GenerateRandomTensor::~GenerateRandomTensor() {
}

/**
 * \brief Testing environement
 */
void GenerateRandomTensor::run() {

  int nv(getIntegerArgument("Nv")), no(getIntegerArgument("No")); 
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *C(new Tensor<>(4, lens, syms, *Cc4s::world, "C"));
  setRandomTensor(*C);
  allocatedTensorArgument("Result", C);
  
}