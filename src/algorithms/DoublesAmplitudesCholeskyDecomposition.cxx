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
    NL = static_cast<int>(lens[0] * reduction + 0.5);
  }

  // write singular vectors back to CTF
  int Ulens[3] = { lens[0], lens[1], NvNo };
  Tensor<> *UaiF(new Tensor<>(3, Ulens, Tabij->sym, *Tabij->wrld, "UaiF"));
  scaU->write(*UaiF);
  delete scaU;

  // identify largest NL |eigenvalues|
  int lower(0), upper(NvNo-1);
  while (NvNo-1 - upper + lower < NL) {
    if (std::abs(S[lower]) > std::abs(S[upper])) {
      ++lower;
    } else {
      --upper;
    }
  }
  LOG(1, "DoublesAmplitudesCholeskyDecomposition") <<
    "taken values: 0<=i<" << lower << " and " << upper+1 << "<=i<" << NvNo <<
    std::endl;

  // write conj(sqrt(Sigmas)) back to CTF
  int64_t sigmaCount( UaiF->wrld->rank == 0 ? NvNo : 0 );
  int64_t *sigmaIndices(new int64_t[sigmaCount]);
  complex *sigmas(new complex[sigmaCount]);
  for (int64_t i(0); i < sigmaCount; ++i) {
    sigmaIndices[i] = i;
    sigmas[i] = conj(sqrt(complex(S[i])));
  }
  Tensor<complex> *Sigma(new Tensor<complex>(1, &NvNo, UaiF->sym, *UaiF->wrld));
  Sigma->set_name("Sigma");
  Sigma->write(sigmaCount, sigmaIndices, sigmas);
  delete[] sigmaIndices;
  delete[] sigmas;
  delete[] S;

  // slice singular vectors U and values S corresponding to NL largest values
  int LaiFlens[] = { lens[0], lens[1], NL };
  Tensor<> *LaiF(new Tensor<>(3, LaiFlens, UaiF->sym, *UaiF->wrld, "LaiF"));
  Tensor<complex> *SF(new Tensor<complex>(1, &NL,UaiF->sym, *UaiF->wrld, "SF"));
  {
    int lowerStart[] = { 0, 0, 0 };
    int lowerEnd[] = { lens[0], lens[1], lower };
    LaiF->slice(lowerStart, lowerEnd, 0.0, *UaiF, lowerStart, lowerEnd, 1.0);
    SF->slice(
      &lowerStart[2],&lowerEnd[2], 0.0, *Sigma, &lowerStart[2],&lowerEnd[2], 1.0
    );
    int upperStart[] = { 0, 0, upper+1 };
    int upperEnd[] = { lens[0], lens[1], NvNo };
    int targetStart[] = { 0, 0, lower };
    LaiF->slice(targetStart, LaiFlens, 0.0, *UaiF, upperStart, upperEnd, 1.0);
    SF->slice(
      &targetStart[2],&LaiFlens[2],0.0, *Sigma, &upperStart[2],&upperEnd[2], 1.0
    );
  }
  delete UaiF;
  delete Sigma;

  Tensor<complex> *cLaiF(
    new Tensor<complex>(3, LaiFlens, LaiF->sym, *LaiF->wrld, "cLaiF")
  );
  toComplexTensor(*LaiF, *cLaiF);
  delete LaiF;

  int LFailens[] = { NL, lens[0], lens[1] };
  Tensor<complex> *LFai(
    new Tensor<complex>(3, LFailens, cLaiF->sym, *cLaiF->wrld, "LFai")
  );
  (*LFai)["Fai"] = (*cLaiF)["aiF"] * (*SF)["F"];
  delete cLaiF;
  delete SF;

  allocatedTensorArgument<complex>("DoublesAmplitudesVertex", LFai);
}

void DoublesAmplitudesCholeskyDecomposition::dryRun() {
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
    NL = static_cast<int>(Tabij->lens[0] * reduction + 0.5);
  }
}

