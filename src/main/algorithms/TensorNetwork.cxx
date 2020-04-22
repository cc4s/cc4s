#include <algorithms/TensorNetwork.hpp>

#include <Cc4s.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorNetwork);

/**
 * \brief Testing environement
 */
Ptr<MapNode> TensorNetwork::run(const Ptr<MapNode> &arguments) {
  // optional argument
  auto spins(arguments->getValue<int64_t>("spins", 2));
  // tensor meta data
  auto matrix(arguments->getMap("matrix"));
  Real<> trace;
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    trace = getTrace<DryTensorEngine>(matrix);
  } else {
    trace = getTrace<DefaultTensorEngine>(matrix);
  }

  LOG(1,"TensorNetwork") << "spins=" << spins << std::endl;
  LOG(1,"TensorNetwork") << "trace=" << trace << std::endl;

  // build result
  auto result(New<MapNode>());
  result->setValue<Real<>>("trace", trace);
  return result;
}

template <typename TE>
Real<> TensorNetwork::getTrace(const Ptr<MapNode> &matrix) {
  typedef Tensor<Real<>,TE> T;
  auto matrixData(matrix->getValue<Ptr<T>>("data"));
  Assert(matrixData, "expecting matrix to be real");
  auto scalar( Tcc<TE>::template tensor<Real<>>(std::vector<size_t>({}), "s") );
  // get trace
  auto getTrace(
    ((*scalar)[""] <<= (*matrixData)["ij"])->compile()
  );
  getTrace->execute();

  // get smaller dimension:
  auto m( std::min(matrixData->getLens()[0],matrixData->getLens()[1]) );
  // set matrix:
  ((*(*matrixData)({0,0},{m,m}))["ii"] += (*scalar)[""])->compile()->execute();
/*
  auto distribute(
    ((*matrixData)["ii"] <<= (*scalar)[""])->compile()
  );
  distribute->execute();
  getTrace->execute();
  distribute->execute();
*/
  getTrace->execute();
  return scalar->read();
}

/*
  if (Cc4s::options->dryRun) {
    size_t No(10);
    size_t Nv(90);
    size_t Np(No+Nv);
    size_t NF(200);
    size_t NR(400);
    typedef Tcc<DryTensorEngine> TCC;
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
      (*PiT)["Rb"] <<= map(conj<Real<>>, (*Pi)["Rb"]),
      (*D)["abij"] +=
        (*T)["cdij"] *
        (*Pi)["Rd"] *
        map(conj<Real<>>, (*Pi)["Rb"]) *
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
  }
*/
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

