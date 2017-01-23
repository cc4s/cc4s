#include <algorithms/DoublesAmplitudesCholeskyDecomposition.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;
using std::shared_ptr;
using std::make_shared;

ALGORITHM_REGISTRAR_DEFINITION(DoublesAmplitudesCholeskyDecomposition);

DoublesAmplitudesCholeskyDecomposition::DoublesAmplitudesCholeskyDecomposition(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DoublesAmplitudesCholeskyDecomposition::
  ~DoublesAmplitudesCholeskyDecomposition()
{
}

void DoublesAmplitudesCholeskyDecomposition::run() {
  diagonlizeAmplitudes();
//  sliceLargestEigenValues();
}

void DoublesAmplitudesCholeskyDecomposition::dryRun() {
  // Read the Coulomb vertex GammaGai
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<>>("DoublesAmplitudes")
  );

  Nv = Tabij->lens[0]*Tabij->lens[1];
  // get number of field variables
  int NL(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NL == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NL = static_cast<int>(Nv * reduction + 0.5);
  }
}

void DoublesAmplitudesCholeskyDecomposition::diagonlizeAmplitudes() {
  Tensor<> *Tabij(getTensorArgument<>("DoublesAmplitudes"));

  // reorder to Taibj
  Nv = Tabij->lens[0]; No = Tabij->lens[2];
  int lens[] = { Nv,No, Nv,No };
  Taibj = make_shared<Tensor<>>(4, lens, Tabij->sym, *Tabij->wrld, "Taibj");
  (*Taibj)["aibj"] = (*Tabij)["abij"];
  // release unneeded resources early
  // Taibj.reset();

  // T(ai)(bj) = U.Lambda.U^T, seen as a matrix with compound indices
  BlacsWorld world(Taibj->wrld->rank, Taibj->wrld->np);
  NvNo = Nv*No;
  int scaTLens[2] = { NvNo, NvNo };
  auto scaTaibj(make_shared<ScaLapackMatrix<>>(*Taibj, scaTLens, &world));

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(make_shared<ScaLapackMatrix<>>(*scaTaibj));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTaibj, scaU);
  lambdas = new double[NvNo];
  eigenSystem.solve(lambdas);
  scaTaibj.reset();

  // write matrix U(ai)(F) back to CTF as tensor UaiF
  int ULens[3] = { Nv, No, NvNo };
  UaiF = make_shared<Tensor<>>(3, ULens, Tabij->sym, *Tabij->wrld, "UaiF");
  scaU->write(*UaiF);
  scaU.reset();

  // write Lambda and conj(sqrt(Lambda)) back to CTF
  lambdasCount = UaiF->wrld->rank == 0 ? NvNo : 0;
  lambdaIndices = new int64_t[lambdasCount];
  sqrtLambdas = new complex[lambdasCount];
  for (int64_t i(0); i < lambdasCount; ++i) {
    lambdaIndices[i] = i;
    sqrtLambdas[i] = conj(sqrt(complex(lambdas[i])));
  }
  sqrtLambdaF = make_shared<Tensor<complex>>(1, &NvNo, UaiF->sym, *UaiF->wrld);
  sqrtLambdaF->set_name("sqrtLambda");
  sqrtLambdaF->write(lambdasCount, lambdaIndices, sqrtLambdas);
  LambdaF = new Tensor<>(1, &NvNo, UaiF->sym, *UaiF->wrld);
  LambdaF->set_name("Lambda");
  LambdaF->write(lambdasCount, lambdaIndices, lambdas);
  delete[] lambdaIndices;
  delete[] sqrtLambdas;
  allocatedTensorArgument<>("DoublesAmplitudesEigenValues", LambdaF);


  // test with the complex tensors
  Tensor<complex> cTaibj(4, Taibj->lens, Taibj->sym, *Taibj->wrld, "cTaibj");
  toComplexTensor(*Taibj, cTaibj);

  // recompose Taibj for testing
  (*Taibj)["aibj"] += (-1.0) * (*UaiF)["aiF"] * (*LambdaF)["F"] * (*UaiF)["bjF"];
  double norm(frobeniusNorm(*Taibj));
  LOG(1, "DoublesAmplitudesCholeskyDecomposition") << "|T-U.Lambda.U*|=" << norm <<
    std::endl;

  Tensor<complex> cUaiF(3, UaiF->lens, UaiF->sym, *UaiF->wrld, "cUaiF");
  toComplexTensor(*UaiF, cUaiF);
  Tensor<complex> conjSqrtLambdaF(false, *sqrtLambdaF);
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  conjSqrtLambdaF.sum(1.0, *sqrtLambdaF,"F", 0.0,"F", fConj);
  cTaibj["aibj"] +=
    (1.0) * cUaiF["aiF"] * conjSqrtLambdaF["F"] *
    (*sqrtLambdaF)["F"] * cUaiF["bjF"];
  double cNorm(frobeniusNorm(cTaibj));
  LOG(1, "DoublesAmplitudesCholeskyDecomposition") << "|T-L*.L|=" << cNorm <<
    std::endl;
}

void DoublesAmplitudesCholeskyDecomposition::sliceLargestEigenValues() {
  // get number of field variables
  NF = getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES);
  // if fieldVariables not given use reduction
  if (NF == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NF = static_cast<int>(Nv * reduction + 0.5);
  }
  NF = std::min(std::max(0, NF), NvNo);

  // identify largest NF |eigenvalues|
  int lower(0), upper(NvNo-1);
  while (NvNo-1 - upper + lower < NF) {
    if (std::abs(lambdas[lower]) > std::abs(lambdas[upper])) {
      ++lower;
    } else {
      --upper;
    }
  }
  delete[] lambdas;
  LOG(1, "DoublesAmplitudesCholeskyDecomposition") <<
    "taken values: 0<=i<" << lower << " and " << upper+1 << "<=i<" << NvNo <<
    std::endl;

  // slice eigen vectors U and values Lambda corresponding to NF largest values
  int LaiFLens[] = { Nv, No, NF };
  auto LaiF(make_shared<Tensor<>>(3, LaiFLens, UaiF->sym, *UaiF->wrld, "LaiF"));
  auto LF(make_shared<Tensor<complex>>(1, &NF, UaiF->sym, *UaiF->wrld, "LF"));
  {
    int lowerStart[] = { 0, 0, 0 };
    int lowerEnd[] = { Nv, No, lower };
    LaiF->slice(lowerStart, lowerEnd, 0.0, *UaiF, lowerStart, lowerEnd, 1.0);
    LF->slice(
      &lowerStart[2], &lowerEnd[2], 0.0,
      *sqrtLambdaF, &lowerStart[2], &lowerEnd[2], 1.0
    );
    int upperStart[] = { 0, 0, upper+1 };
    int upperEnd[] = { Nv, No, NvNo };
    int targetStart[] = { 0, 0, lower };
    LaiF->slice(targetStart, LaiFLens, 0.0, *UaiF, upperStart, upperEnd, 1.0);
    LF->slice(
      &targetStart[2], &LaiFLens[2], 0.0,
      *sqrtLambdaF, &upperStart[2], &upperEnd[2], 1.0
    );
  }
  UaiF.reset();
  sqrtLambdaF.reset();

  auto cLaiF(
    make_shared<Tensor<complex>>(3, LaiFLens, LaiF->sym, *LaiF->wrld, "cLaiF")
  );
  toComplexTensor(*LaiF, *cLaiF);
  LaiF.reset();

  int LFaiLens[] = { NF, Nv, No };
  auto LFai(
    new Tensor<complex>(3, LFaiLens, cLaiF->sym, *cLaiF->wrld, "LFai")
  );
  (*LFai)["Fai"] = (*cLaiF)["aiF"] * (*LF)["F"];
  cLaiF.reset();
  LF.reset();

  allocatedTensorArgument<complex>("DoublesAmplitudesVertex", LFai);

  // and now from the complex tensors
  Tensor<complex> cTaibj(4, Taibj->lens, Taibj->sym, *Taibj->wrld, "cTaibj");
  toComplexTensor(*Taibj, cTaibj);
  Tensor<complex> conjLFai(false, *LFai);
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  conjLFai.sum(1.0, *LFai,"Fai", 0.0,"Fai", fConj);
  cTaibj["aibj"] += (-1.0) * conjLFai["Fai"] * (*LFai)["Fbj"];
  double cNorm(frobeniusNorm(cTaibj));
  LOG(1, "DoublesAmplitudesCholeskyDecomposition") << "|T-L*.L|=" << cNorm <<
    std::endl;
}
