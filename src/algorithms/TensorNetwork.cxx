#include <algorithms/TensorNetwork.hpp>
#include <tcc/DryTensor.hpp>
#include <tcc/DryTensorTransaction.hpp>

#include <ctf.hpp>
#include <array>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorNetwork);

TensorNetwork::TensorNetwork(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorNetwork::~TensorNetwork() {
}

/**
 * \brief Testing environement
 */
void TensorNetwork::run() {
}


void TensorNetwork::dryRun() {
  DryTensor<> T(
    4, std::array<int,4>{{100,100,10,10}}.data(),
    std::array<int,4>{{NS,NS,NS,NS}}.data()
  );
  T.set_name("T");
  DryTensor<> Pi(
    2, std::array<int,2>{{300,100}}.data(),
    std::array<int,2>{{NS,NS}}.data()
  );
  Pi.set_name("Pi");
  DryTensor<> PiT(
    2, std::array<int,2>{{300,100}}.data(),
    std::array<int,2>{{NS,NS}}.data()
  );
  PiT.set_name("PiT");
  DryTensor<> Lambda(
    2, std::array<int,2>{{300,200}}.data(),
    std::array<int,2>{{NS,NS}}.data()
  );
  Lambda.set_name("Lambda");
  DryTensor<> LambdaT(
    2, std::array<int,2>{{300,200}}.data(),
    std::array<int,2>{{NS,NS}}.data()
  );
  LambdaT.set_name("LambdaT");

//  CompoundDryTensorExpression<> Gamma("Fac") = PiT["Ra"] * Pi["Rc"] * Lambda["RG"]

  DryTensorTransaction<> ladder(
    T["abij"] = T["cdij"] *
      PiT["Sa"] * PiT["Rb"] * LambdaT["SF"] * Lambda["RF"] * Pi["Rd"] * Pi["Sc"]
  );
}


