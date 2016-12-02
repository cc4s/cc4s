#include <algorithms/TensorNetwork.hpp>
#include <tcc/DryTensor.hpp>
#include <tcc/TensorContraction.hpp>
#include <tcc/TensorAssignment.hpp>
#include <tcc/TensorOperation.hpp>
#include <tcc/TensorAssignmentOperation.hpp>
#include <tcc/TensorContractionOperation.hpp>

#include <array>
#include <memory>
using std::shared_ptr;

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
    std::array<int,4>{{0,0,0,0}}.data()
  );
  T.set_name("T");
  DryTensor<> Pi(
    2, std::array<int,2>{{300,100}}.data(),
    std::array<int,2>{{0,0}}.data()
  );
  Pi.set_name("Pi");
  DryTensor<> PiT(
    2, std::array<int,2>{{300,100}}.data(),
    std::array<int,2>{{0,0}}.data()
  );
  PiT.set_name("PiT");
  DryTensor<> Lambda(
    2, std::array<int,2>{{300,200}}.data(),
    std::array<int,2>{{0,0}}.data()
  );
  Lambda.set_name("Lambda");
  DryTensor<> LambdaT(
    2, std::array<int,2>{{300,200}}.data(),
    std::array<int,2>{{0,0}}.data()
  );
  LambdaT.set_name("LambdaT");
  DryTensor<> Gamma(
    3, std::array<int,3>{{200, 100,100}}.data(),
    std::array<int,3>{{0,0,0}}.data()
  );
  Gamma.set_name("Gamma");

//  CompoundDryTensorExpression<> Gamma("Fac") = PiT["Ra"] * Pi["Rc"] * Lambda["RG"]

  shared_ptr<TensorOperation<>> factorOperation = compile(
    Gamma["Fqr"] = Lambda["RF"] * PiT["Rq"] * Pi["Rr"]
  );
  factorOperation->execute();

  shared_ptr<TensorOperation<>> ladderOperation = compile(
    T["abij"] =
      T["cdij"] * Pi["Rd"]  * PiT["Rb"] *
      Pi["Sc"] * PiT["Sa"] * LambdaT["SF"] * Lambda["RF"]
  );
  ladderOperation->execute();
}


