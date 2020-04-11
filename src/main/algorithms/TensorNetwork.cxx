#include <algorithms/TensorNetwork.hpp>

#include <tcc/Tcc.hpp>
#include <tcc/engines/DryMachineTensor.hpp>
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
  typedef Tcc<DryEngine> TCC;
  auto T(
    TCC::tensor(std::vector<size_t>({Nv,Nv,No,No}), "T")
  );
  auto D(
    TCC::tensor(std::vector<size_t>({Nv,Nv,No,No}), "D")
  );
  auto Pi(
    TCC::tensor(std::vector<size_t>({NR,Nv}), "Pi")
  );
  auto PiT(
    TCC::tensor(std::vector<size_t>({NR,Nv}), "PiT")
  );
  auto Lambda(
    TCC::tensor(std::vector<size_t>({NR,NF}), "Lambda")
  );
  auto LambdaT(
    TCC::tensor(std::vector<size_t>({NR,NF}), "LambdaT")
  );
  auto Gamma(
    TCC::tensor(std::vector<size_t>({NF,Nv,Nv}), "Gamma")
  );
  auto Pir(
    TCC::tensor(std::vector<size_t>({NR,Np}), "Pir")
  );

  // compile a sequence (,) of operations. Note the required parenthesis
  auto ladderExpression = (
    (*LambdaT)["RF"] <<= (*Lambda)["RF"],
    (*PiT)["Rb"] <<= tcc::map(conj<Real<>>, (*Pi)["Rb"]),
    (*D)["abij"] +=
      (*T)["cdij"] *
      (*Pi)["Rd"] *
      tcc::map(conj<Real<>>, (*Pi)["Rb"]) *
      (*Pi)["Sc"] * (*PiT)["Sa"] *
      (*LambdaT)["SF"] * (*Lambda)["RF"],
    (*Pi)["Ra"] <<= (*Pi)["Ra"],
    (*Pi)["Ra"] <<= (*(*Pir)({0,No},{NR,Np}))["Ra"],
    (*(*Pir)({0,No},{NR,Np}))["Ra"] -= (*Pi)["Ra"]
  );
  LOG(1,"TensorNetwork") << "Expression = " <<
    std::string(*ladderExpression) << std::endl;
  auto ladderOperation(ladderExpression->compile());
  LOG(1,"TensorNetwork") << "Operation = " <<
    std::string(*ladderOperation) << std::endl;
  ladderOperation->execute();

// this contraction already requires heuristics
/*
  shared_ptr<Tensor<>> Pia(
    TCC::tensor(std::vector<size_t>({NR,Nv}), "Pia")
  );
  shared_ptr<Tensor<>> Pii(
    TCC::tensor(std::vector<size_t>({NR,No}), "Pii")
  );
  size_t Nn(7);
  shared_ptr<Tensor<>> w(
    TCC::tensor(std::vector<size_t>({Nn}), "w")
  );
  shared_ptr<Tensor<>> H(
    TCC::tensor(std::vector<size_t>({No,Nn}), "H")
  );
  shared_ptr<Tensor<>> P(
    TCC::tensor(std::vector<size_t>({Nv,Nn}), "P")
  );
  shared_ptr<Tensor<>> e(
    TCC::tensor(std::vector<size_t>(), "e")
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


