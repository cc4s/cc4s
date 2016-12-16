#include <algorithms/TensorNetwork.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <Cc4s.hpp>

#include <vector>
#include <memory>
using std::shared_ptr;

using namespace cc4s;
using namespace tcc;

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
  int No(10);
  int Nv(90);
  int Np(No+Nv);
  int NF(200);
  int NR(400);
  auto machineTensorFactory(
//    CtfMachineTensorFactory<>::create(Cc4s::world)
    DryMachineTensorFactory<>::create()
  );
  auto tcc(Tcc<>::create(machineTensorFactory));

/*
  shared_ptr<Tensor<complex>> Tc(
    tcc.createTensor<complex>(std::vector<int>({Nv,Nv,No,No}), "Tc")
  );
*/
  auto T(
    tcc->createTensor(std::vector<int>({Nv,Nv,No,No}), "T")
  );
  auto D(
    tcc->createTensor(std::vector<int>({Nv,Nv,No,No}), "D")
  );
  auto Pi(
    tcc->createTensor(std::vector<int>({NR,Np}), "Pi")
  );
  auto PiT(
    tcc->createTensor(std::vector<int>({NR,Np}), "PiT")
  );
  auto Lambda(
    tcc->createTensor(std::vector<int>({NR,NF}), "Lambda")
  );
  auto LambdaT(
    tcc->createTensor(std::vector<int>({NR,NF}), "LambdaT")
  );
  auto Gamma(
    tcc->createTensor(std::vector<int>({NF,Np,Np}), "Gamma")
  );

//  CompoundDryTensorExpression<> Gamma("Fac") = PiT["Ra"] * Pi["Rc"] * Lambda["RG"]

  // compile a sequence (,) of operations. Note the required parenthesis
  auto ladderOperation = tcc->compile(
    (
/*
      (*Gamma)["Fqr"] <<= (*Pi)["Rq"] * (*Pi)["Rr"] * (*Lambda)["RF"],
      (*Pi)["Rr"] <<= (*LambdaT)["RF"] * (*PiT)["Rq"] * (*Gamma)["Fqr"],
      (*T)["abij"] -= -1/4. *(*T)["abji"],
*/
      (*D)["abij"] +=
        (*T)["cdij"] *
        (*Pi)["Rd"] * (*PiT)["Rb"] *
        (*Pi)["Sc"] * (*PiT)["Sa"] *
        (*LambdaT)["SF"] * (*Lambda)["RF"]
    )
  );
  ladderOperation->execute();

// this contraction already requires heuristics
/*
  shared_ptr<Tensor<>> Pia(
    tcc->createTensor(std::vector<int>({NR,Nv}), "Pia")
  );
  shared_ptr<Tensor<>> Pii(
    tcc->createTensor(std::vector<int>({NR,No}), "Pii")
  );
  int Nn(7);
  shared_ptr<Tensor<>> w(
    tcc->createTensor(std::vector<int>({Nn}), "w")
  );
  shared_ptr<Tensor<>> H(
    tcc->createTensor(std::vector<int>({No,Nn}), "H")
  );
  shared_ptr<Tensor<>> P(
    tcc->createTensor(std::vector<int>({Nv,Nn}), "P")
  );
  shared_ptr<Tensor<>> e(
    tcc->createTensor(std::vector<int>(), "e")
  );

  auto imaginaryTimeMp2Operation = tcc->compile(
    (*e)[""] <<=
      (*Pii)["Ri"]  * (*Pia)["Ra"] *
        (*LambdaT)["RF"] * (*Lambda)["SF"] *
      (*Pii)["Sj"] * (*Pia)["Sb"] *
        (*w)["n"] * (*P)["an"] * (*H)["in"] * (*P)["bn"] * (*H)["jn"] *
      (*Pii)["Ti"]  * (*Pia)["Ta"] *
        (*LambdaT)["TH"] * (*Lambda)["UH"] *
      (*Pii)["Uj"] * (*Pia)["Ub"]
  );
  imaginaryTimeMp2Operation->execute();
*/
}


