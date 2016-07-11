#include <algorithms/EnergyMatrixFromDoublesAmplitudes.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(EnergyMatrixFromDoublesAmplitudes);

EnergyMatrixFromDoublesAmplitudes::EnergyMatrixFromDoublesAmplitudes(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

EnergyMatrixFromDoublesAmplitudes::~EnergyMatrixFromDoublesAmplitudes() {
}

void EnergyMatrixFromDoublesAmplitudes::run() {
  Tensor<> *Tabij(getTensorArgument("DoublesAmplitudes"));
 
  // TODO: use complex conversion routines
  Tensor<> imagTabij(false, *Tabij);
  Tensor<complex> cTabij(4, Tabij->lens, Tabij->sym, *Tabij->wrld, "cTabij");
  // convert into complex tensor
  toComplexTensor(*Tabij, imagTabij, cTabij);

  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  // calculate conjugate of GammaGai
  Tensor<complex> conjGammaGai(false, *GammaGai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaGai.sum(1.0, *GammaGai,"Gai", 0.0,"Gai", fConj);

  Matrix<complex> *energyMatrix = new Matrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], *Tabij->wrld
  );
  allocatedTensorArgument<complex>("EnergyMatrix", energyMatrix);

  (*energyMatrix)["GH"] =
    conjGammaGai["Gai"] * (*GammaGai)["Hbj"] * cTabij["abij"];
  (*energyMatrix)["GH"] *= 2.0;

  Scalar<complex> energy(*energyMatrix->wrld);
  energy[""] = (*energyMatrix)["GG"];
  complex e(energy.get_val());
  LOG(1, "EMAT") << "Tr{E}=" << e << std::endl;
}

void EnergyMatrixFromDoublesAmplitudes::dryRun() {
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("DoublesAmplitudes")
  );
 
  // TODO: use complex conversion routines
  DryTensor<> imagTabij(*Tabij);
  DryTensor<complex> cTabij(4, Tabij->lens.data(), Tabij->syms.data());

  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );

  // calculate conjugate of GammaGai
  DryTensor<complex> conjGammaGai(*GammaGai);

  DryMatrix<complex> *energyMatrix = new DryMatrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "EnergyMatrix", energyMatrix
  );

  // Allocate intermediate (implicitly used)
  int syms[] = { NS, NS, NS };
  int Gai[] = { GammaGai->lens[0], Tabij->lens[0], Tabij->lens[2] };
  DryTensor<complex> GammaTGai(3, Gai, syms);
}

