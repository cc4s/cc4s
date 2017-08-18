#include <algorithms/Mp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EnergyFromCoulombIntegrals);

Mp2EnergyFromCoulombIntegrals::Mp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2EnergyFromCoulombIntegrals::~Mp2EnergyFromCoulombIntegrals() {
}

void Mp2EnergyFromCoulombIntegrals::run() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  TensorData<complex> *complexVabij(dynamic_cast<TensorData<complex> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = calculateMp2Energy<double>(*realVabij->value);
  } else {
    e = std::real(calculateMp2Energy<complex>(*complexVabij->value));
  }
  setRealArgument("Mp2Energy", e);
}

// FIXME: update dryRun to work in the complex case as well
void Mp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
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

template <typename F>
F Mp2EnergyFromCoulombIntegrals::calculateMp2Energy(CTF::Tensor<F> &Vabij) {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));

  // convert to type F (either complex or double)
  Tensor<F> Fepsi(1, &epsi->lens[0], epsi->sym, *epsi->wrld, "Fepsi");
  // NOTE: just copies if both arguments are real
  toComplexTensor(*epsi, Fepsi);
  Tensor<F> Fepsa(1, &epsa->lens[0], epsa->sym, *epsa->wrld, "Fepsa");
  toComplexTensor(*epsa, Fepsa);

  // create excitation energy
  auto Tabij(new Tensor<F>(Vabij));
  //  Tensor<F> Tabij(false, Vabij);
  Tabij->set_name("Tabij");
  (*Tabij)["abij"] =  Fepsi["i"];
  (*Tabij)["abij"] += Fepsi["j"];
  (*Tabij)["abij"] -= Fepsa["a"];
  (*Tabij)["abij"] -= Fepsa["b"];

  // use transform to divide Vabij by Tabij and store in Tabij
  CTF::Transform<F, F>(
    std::function<void(F, F &)>(
      [](F vabij, F &tabij) {
        tabij = conj(vabij / tabij);
      }
    )
  ) (
     Vabij["abij"], (*Tabij)["abij"]
  );

// TODO: below method requires less memory but does not work with current CTF
/*
  Matrix<> Dai(Vabij.lens[0], Vabij.lens[2], NS, *Vabij.wrld, "Dai");
  Dai["ai"] =  (*epsi)["i"];
  Dai["ai"] -= (*epsa)["a"];

  Tensor<F> Tabij(Vabij);
  // use transform to divide real/complex Tabij by real Dabij = Dai+Dbj
  CTF::Transform<double, double, F>(
    std::function<void(double, double, F &)>(
      [](double dai, double dbj, F &t) {
        t = conj(t / (dai + dbj));
      }
    )
  ) (
    Dai["ai"], Dai["bj"], Tabij["abij"]
  );
*/

  Scalar<F> energy(*Cc4s::world);
  F dire, exce;

  energy[""] = 2.0 * (*Tabij)["abij"] * Vabij["abij"];
  dire = energy.get_val();
  energy[""] = (*Tabij)["abji"] * Vabij["abij"];
  exce = -1.0 * energy.get_val();
  LOG(0, "MP2") << "e=" << (dire + exce) << std::endl;
  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  if (isArgumentGiven("Mp2DoublesAmplitudes")) {
    allocatedTensorArgument<F>("Mp2DoublesAmplitudes", Tabij);
  } else {
    delete Tabij;
  }

  return dire + exce;
}

