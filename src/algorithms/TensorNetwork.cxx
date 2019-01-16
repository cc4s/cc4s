#include <algorithms/TensorNetwork.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
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
  size_t No(10);
  size_t Nv(90);
  size_t Np(No+Nv);
  size_t NF(200);
  size_t NR(400);
  auto machineTensorFactory(
//    CtfMachineTensorFactory<>::create(Cc4s::world)
    DryMachineTensorFactory<>::create()
  );
  auto tcc(Tcc<real>::create(machineTensorFactory));

/*
  shared_ptr<Tensor<complex>> Tc(
    tcc.createTensor<complex>(std::vector<size_t>({Nv,Nv,No,No}), "Tc")
  );
*/
  auto T(
    tcc->createTensor(std::vector<size_t>({Nv,Nv,No,No}), "T")
  );
  auto D(
    tcc->createTensor(std::vector<size_t>({Nv,Nv,No,No}), "D")
  );
  auto Pi(
    tcc->createTensor(std::vector<size_t>({NR,Np}), "Pi")
  );
  auto PiT(
    tcc->createTensor(std::vector<size_t>({NR,Np}), "PiT")
  );
  auto Lambda(
    tcc->createTensor(std::vector<size_t>({NR,NF}), "Lambda")
  );
  auto LambdaT(
    tcc->createTensor(std::vector<size_t>({NR,NF}), "LambdaT")
  );
  auto Gamma(
    tcc->createTensor(std::vector<size_t>({NF,Np,Np}), "Gamma")
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
        (*Pi)["Rd"] * map<real>(conj<real>, (*Pi)["Rb"]) *
        (*Pi)["Sc"] * (*PiT)["Sa"] *
        (*LambdaT)["SF"] * (*Lambda)["RF"]
    )
  );
  ladderOperation->execute();

// this contraction already requires heuristics
/*
  shared_ptr<Tensor<>> Pia(
    tcc->createTensor(std::vector<size_t>({NR,Nv}), "Pia")
  );
  shared_ptr<Tensor<>> Pii(
    tcc->createTensor(std::vector<size_t>({NR,No}), "Pii")
  );
  size_t Nn(7);
  shared_ptr<Tensor<>> w(
    tcc->createTensor(std::vector<size_t>({Nn}), "w")
  );
  shared_ptr<Tensor<>> H(
    tcc->createTensor(std::vector<size_t>({No,Nn}), "H")
  );
  shared_ptr<Tensor<>> P(
    tcc->createTensor(std::vector<size_t>({Nv,Nn}), "P")
  );
  shared_ptr<Tensor<>> e(
    tcc->createTensor(std::vector<size_t>(), "e")
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


