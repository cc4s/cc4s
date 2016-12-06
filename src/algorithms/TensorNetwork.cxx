#include <algorithms/TensorNetwork.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Assignment.hpp>
#include <tcc/Operation.hpp>
#include <tcc/AssignmentOperation.hpp>
#include <tcc/ContractionOperation.hpp>

#include <array>
#include <memory>
using std::shared_ptr;
using std::make_shared;

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
  int NR(300);
  Tcc tcc;

/*
  shared_ptr<Tensor<complex>> Tc(
    tcc.createTensor<complex>(std::vector<int>{{100,100,10,10}}, "Tc")
  );
*/
  shared_ptr<Tensor<>> T(
    tcc.createTensor<>(std::vector<int>{{100,100,10,10}}, "T")
  );
  shared_ptr<Tensor<>> Pi(
    tcc.createTensor<>(std::vector<int>{{300,100}}, "Pi")
  );
  shared_ptr<Tensor<>> PiT(
    tcc.createTensor<>(std::vector<int>{{300,100}}, "PiT")
  );
  shared_ptr<Tensor<>> Lambda(
    tcc.createTensor<>(std::vector<int>{{300,200}}, "Lambda")
  );
  shared_ptr<Tensor<>> LambdaT(
    tcc.createTensor<>(std::vector<int>{{300,200}}, "LambdaT")
  );

//  CompoundDryTensorExpression<> Gamma("Fac") = PiT["Ra"] * Pi["Rc"] * Lambda["RG"]

//  shared_ptr<Expression<>> ladderExpression((*Pi)["rq"] <<= (*Pi)["Rr"] * (*Pi)["Rq"]);

  shared_ptr<Operation<>> ladderOperation = compile(
    (*T)["abij"] <<=
      (*T)["cdij"] * (*Pi)["Rd"]  * (*PiT)["Rb"] *
      (*Pi)["Sc"] * (*PiT)["Sa"] * (*LambdaT)["SF"] * (*Lambda)["RF"]
  );
  ladderOperation->execute();
// this contraction already requires heuristics
/*
  Tensor<> Pia(std::vector<int>{{NR,Nv}}, "Pia");
  Tensor<> Pii(std::vector<int>{{NR,No}}, "Pii");
  int Nn(7);
  Tensor<> w(std::vector<int>{{Nn}}, "w");
  Tensor<> H(std::array<int>{{No,Nn}}, "H");
  Tensor<> P(std::array<int>{{Nv,Nn}}, "P");
  Tensor<> e(std::vector<int>(), "e");
  
  shared_ptr<Operation<>> imaginaryTimeMp2Operation = compile(
    e[""] <<=
      Pii["Ri"]  * Pia["Ra"] *
        LambdaT["RF"] * Lambda["SF"] *
      Pii["Sj"] * Pia["Sb"] *
        w["n"] * P["an"] * H["in"] * P["bn"] * H["jn"] *
      Pii["Ti"]  * Pia["Ta"] *
        LambdaT["TH"] * Lambda["UH"] *
      Pii["Uj"] * Pia["Ub"]
  );
  imaginaryTimeMp2Operation->execute();
*/
}


