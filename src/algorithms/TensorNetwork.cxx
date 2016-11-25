#include <algorithms/TensorNetwork.hpp>
#include <util/DryTensor.hpp>
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
  DryTensor<> V(
    4, std::array<int,4>{{100,100,100,100}}.data(),
    std::array<int,4>{{NS,NS,NS,NS}}.data()
  );
  DryTensor<> T(
    4, std::array<int,4>{{100,100,10,10}}.data(),
    std::array<int,4>{{NS,NS,NS,NS}}.data()
  );
  T["abji"] = T["cdij"] * V["abcd"];
}


