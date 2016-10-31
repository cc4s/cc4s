#include <algorithms/PerturbativeTriples.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriples);

PerturbativeTriples::PerturbativeTriples(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PerturbativeTriples::~PerturbativeTriples() {
}

void PerturbativeTriples::run() {
  Tensor<>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<>   *Tai(getTensorArgument("CcsdSinglesAmplitudes"));
  
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  int vvvooo[] = { Nv, Nv , Nv , No , No , No };
  int   syms[] = { NS, NS,  NS , NS , NS , NS };
  Tensor<> Zabcijk(6, vvvooo, syms, *Vabij->wrld, "Zabcijk");
  Zabcijk["abcijk"]  = (*Vabci)["bcek"] * (*Tabij)["aeij"];
  Zabcijk["abcijk"] -= (*Vijka)["jima"] * (*Tabij)["bcmk"];

  Tensor<> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");
  Tabcijk["abcijk"]  = (*epsi)["i"];
  Tabcijk["abcijk"] += (*epsi)["j"];
  Tabcijk["abcijk"] += (*epsi)["k"];
  Tabcijk["abcijk"] -= (*epsa)["a"];
  Tabcijk["abcijk"] -= (*epsa)["b"];
  Tabcijk["abcijk"] -= (*epsa)["c"];
  Bivar_Function<> fDivide(&divide<double>);
  Tabcijk.contract(
    1.0, Zabcijk,"abcijk", Tabcijk,"abcijk", 0.0,"abcijk", fDivide
  );

  Zabcijk["abcijk"] -= (*Vabij)["abij"] * (*Tai)["ck"];

  Scalar<> energy(*Cc4s::world);
  double e, triplese;
  double ccsde(getRealArgument("CcsdEnergy"));

  energy[""] = (1.0/36.0) * Tabcijk["abcijk"] * Zabcijk["abcijk"];
  triplese = energy.get_val();
  e = triplese + ccsde;

  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << ccsde << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << triplese << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PPPHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("CcsdSinglesAmplitudes");
  getTensorArgument<double, DryTensor<double>>("CcsdDoublesAmplitudes");

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
  int vvvooo[] = { Nv, Nv , Nv , No , No , No };
  int   syms[] = { NS, NS,  NS , NS , NS , NS };
  DryTensor<> Zabcijk(6, vvvooo, syms);
  DryTensor<> Tabcijk(6, vvvooo, syms);

  DryScalar<> energy();
}

