#include <algorithms/DoublesAmplitudesSingularVectors.hpp>
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

ALGORITHM_REGISTRAR_DEFINITION(DoublesAmplitudesSingularVectors);

DoublesAmplitudesSingularVectors::DoublesAmplitudesSingularVectors(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DoublesAmplitudesSingularVectors::~DoublesAmplitudesSingularVectors() {
}

void DoublesAmplitudesSingularVectors::run() {
  // read the doubles amplitudes
  Tensor<> *Tabij(getTensorArgument<>("DoublesAmplitudes"));

  // reorder to Taibj
  int lens[] = {
    Tabij->lens[0], Tabij->lens[2], Tabij->lens[1], Tabij->lens[3]
  };
  Tensor<> *Taibj( new Tensor<>(4, lens, Tabij->sym, *Tabij->wrld, "Taibj") );
  (*Taibj)["aibj"] = (*Tabij)["abij"];

  // T(ai)(bj) = U Sigma U^T, seen as a matrix with compound indices
  int matrixLens[2] = { lens[0]*lens[1], lens[2]*lens[3] };
  BlacsWorld world(Taibj->wrld->rank, Taibj->wrld->np);
  ScaLapackMatrix<> *scaTaibj(new ScaLapackMatrix<>(*Taibj, matrixLens, &world));
  delete Taibj;

  // use ScaLapack routines to diagonalise the U S UT matrix, i.e. find U
  ScaLapackMatrix<> *scaU(new ScaLapackMatrix<>(*scaTaibj));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTaibj, scaU);
  int NvNo(matrixLens[0]);
  double *S(new double[NvNo]);
  eigenSystem.solve(S);
  delete scaTaibj;

  // get number of field variables
  int NL(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NL == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NL = static_cast<int>(NvNo * reduction + 0.5);
  }

  // write singular vectors back to CTF
  int Ulens[3] = { lens[0], lens[1], NvNo };
  Tensor<> *UaiF(new Tensor<>(3, Ulens, Tabij->sym, *Tabij->wrld, "UaiF"));
  scaU->write(*UaiF);
  delete scaU;

  int64_t sigmaCount( UaiF->wrld->rank == 0 ? NvNo : 0 );
  int64_t *sigmaIndices(new int64_t[sigmaCount]);
  for (int64_t i(0); i < sigmaCount; ++i) sigmaIndices[i] = i;
  Tensor<> *Sigma(new Tensor<>(1, &NvNo, UaiF->sym, *UaiF->wrld));
  Sigma->set_name("Sigma");
  Sigma->write(sigmaCount, sigmaIndices, S);
  delete[] sigmaIndices;

  allocatedTensorArgument<>("DoublesAmplitudesEigenValues", Sigma);
/*
  // slice singular vectors U corresponding to NL largest singular values S
  // TODO: what about negative eigenvalues
  int start[] = {0, 0, NvNo-NL}, end[] = {lens[0], lens[1], NvNo};
  allocatedTensorArgument<>(
    "DoublesAmplitudesSingularVectors", new Tensor<>(UaiF->slice(start, end))
  );
*/
  delete UaiF;
  // TODO: also write out the singular values
  // TODO: take sqrt and do U[aiF] * sqrt(sigma[F])
  // TODO: make it complex
  delete[] S;
}

void DoublesAmplitudesSingularVectors::dryRun() {
  // Read the Coulomb vertex GammaGai
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<>>("DoublesAmplitudes")
  );

  int NvNo(Tabij->lens[0]*Tabij->lens[1]);
  // get number of field variables
  int NL(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (NL == DEFAULT_FIELD_VARIABLES) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    NL = static_cast<int>(NvNo * reduction + 0.5);
  }

  allocatedTensorArgument<double, DryTensor<>>(
    "DoublesAmplitudesEigenValues",
    new DryVector<>(NvNo, SOURCE_LOCATION)
  );
}

