#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalMp2EnergyFromCoulombIntegrals);

ThermalMp2EnergyFromCoulombIntegrals::ThermalMp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalMp2EnergyFromCoulombIntegrals::~ThermalMp2EnergyFromCoulombIntegrals() {
}

void ThermalMp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Ni(getTensorArgument("ThermalHoleOccupancies"));
  Tensor<> *Na(getTensorArgument("ThermalParticleOccupancies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));

  // first, compute \Delta^{ab}_{ij}, note that its sign is opposite than usual
  Tensor<> Tabij(false, Vabij);
  Tabij["abij"] =  (*epsa)["a"];
  Tabij["abij"] += (*epsa)["b"];
  Tabij["abij"] -= (*epsi)["i"];
  Tabij["abij"] -= (*epsi)["j"];

  // then, compute the factor from the time integration for each abij
  // from the \Delta^{ab}_{ij}
  // treat \Delta = 0 specially
  class timeIntegral {
  public:
    timeIntegral(double kT_): kT(kT_) { }
    void operator()(double &eps) {
      eps = std::abs(eps) > 1e-8 ?
        (std::exp(-eps/kT) - 1.0 + eps/kT) / (eps*eps/kT) :
        0.5/kT;
    }
  protected:
    double kT;
  };
  Transform<>(
    std::function<void(double &)>(timeIntegral(getRealArgument("Temperature")))
  ) (
    Tabij["abij"]
  );

  // finally aggregate the Fermi occupancies for a,b,i and j in the
  // same tensor
  Tabij["abij"] *= (*Na)["a"];
  Tabij["abij"] *= (*Na)["b"];
  Tabij["abij"] *= (*Ni)["i"];
  Tabij["abij"] *= (*Ni)["j"];

  // then it can be multiplied with Vabij before continuing with the contractions
  Tabij["abij"] *= (*Vabij)["abij"];

  // do the contractions pretty much the same way as in the zero T case
  // note that the tensors are larger from the overlap
  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * (*Vabij)["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * (*Vabij)["abij"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;

  LOG(0, "FT-MP2") << "FT-MP2=" << e << std::endl;
  LOG(1, "FT-MP2") << "FT-MP2d=" << dire << std::endl;
  LOG(1, "FT-MP2") << "FT-MP2x=" << exce << std::endl;

  setRealArgument("ThermalMp2Energy", e);
}

void ThermalMp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("ThermalPPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}

