#include <algorithms/PerturbativeTriples.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
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

Tensor<> &PerturbativeTriples::getSinglesContribution(int i, int j, int k) {
  Tensor<> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  int TaiStart[] = { 0, i }, TaiEnd[] = { Nv, i+1 };
  int VbcjkStart[] = { 0, 0, j, k }, VbcjkEnd[] = { Nv, Nv, j+1, k+1 };
  (*SVabc)["abc"] =
    0.5 * Tai->slice(TaiStart,TaiEnd)["ai"]
    * Vabij->slice(VbcjkStart,VbcjkEnd)["bcjk"];
  return *SVabc;
}

Tensor<> &PerturbativeTriples::getDoublesContribution(int i, int j, int k) {
  Tensor<> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
  Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  int TadijStart[] = { 0, 0, i, j }, TadijEnd[] = { Nv, Nv, i+1, j+1 };
  int VbcdkStart[] = { 0, 0, 0, k }, VbcdkEnd[] = { Nv, Nv, Nv, k+1 };
  (*DVabc)["abc"] =
    Tabij->slice(TadijStart,TadijEnd)["adij"]
    * Vabci->slice(VbcdkStart,VbcdkEnd)["bcdk"];

  int TabilStart[] = { 0, 0, i, 0 }, TabilEnd[] = { Nv, Nv, i+1, No };
  int VjklcStart[] = { j, k, 0, 0 }, VjklcEnd[] = { j+1, k+1, No, Nv };
  (*DVabc)["abc"] -=
    Tabij->slice(TabilStart,TabilEnd)["abil"]
    * Vijka->slice(VjklcStart,VjklcEnd)["jklc"];

  return *DVabc;
}

Tensor<> &PerturbativeTriples::getEnergyDenominator(int i, int j, int k) {
  // reuse SVabc to hold the energy denominator
  Tensor<>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  int epsiStart[] = { i }, epsiEnd[] = { i+1 };
  (*SVabc)["abc"]  = epsi->slice(epsiStart,epsiEnd)["i"];
  int epsjStart[] = { j }, epsjEnd[] = { j+1 };
  (*SVabc)["abc"] += epsi->slice(epsjStart,epsjEnd)["j"];
  int epskStart[] = { k }, epskEnd[] = { k+1 };
  (*SVabc)["abc"] += epsi->slice(epskStart,epskEnd)["k"];
  (*SVabc)["abc"] -= (*epsa)["a"];
  (*SVabc)["abc"] -= (*epsa)["b"];
  (*SVabc)["abc"] -= (*epsa)["c"];
  return *SVabc;
}

void PerturbativeTriples::run() {
  Tensor<>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  No = epsi->lens[0];
  Nv = epsa->lens[0];
  int vvv[] = { Nv, Nv , Nv };
  int syms[] = { NS, NS,  NS };
  // doubles amplitudes contracted with V for current i,j,k
  DVabc = new Tensor<>(3, vvv, syms, *epsi->wrld, "DVabc");
  // unconnected singles amplitudes and V for current i,j,k
  SVabc = new Tensor<>(3, vvv, syms, *epsi->wrld, "SVabc");
  // triples amplitudes considering all permutations of a,b,c for given i,j,k
  Tensor<> Tabc(3, vvv, syms, *epsi->wrld, "Tabc");

  Scalar<> energy(*Cc4s::world);
  for (int i(0); i < No; ++i) {
    for (int j(0); j < No; ++j) {
      for (int k(0); k < No; ++k) {
        // get DV in respective permutation of a,b,k where i,j,k is
        // attached left to right to D.V
        getDoublesContribution(i,j,k);
        Tabc["abc"]  = (+8.0) * (*DVabc)["abc"];
        Tabc["abc"] += (-4.0) * (*DVabc)["acb"];
        Tabc["abc"] += (-4.0) * (*DVabc)["bac"];
        Tabc["abc"] += (+2.0) * (*DVabc)["bca"];
        Tabc["abc"] += (+2.0) * (*DVabc)["cab"];
        Tabc["abc"] += (-4.0) * (*DVabc)["cba"];

        // get SV in respective permutation of i,j,k where a,b,c is
        // attached left to right to S.V
        // TODO: consider symmetry of SV
/*
        getSinglesContribution(i,j,k);
        Tabc["abc"] += (+8.0) * (*SVabc)["abc"];
        Tabc["abc"] += (-4.0) * (*SVabc)["acb"];
        Tabc["abc"] += (-4.0) * (*SVabc)["bac"];
        Tabc["abc"] += (+2.0) * (*SVabc)["bca"];
        Tabc["abc"] += (+2.0) * (*SVabc)["cab"];
        Tabc["abc"] += (-4.0) * (*SVabc)["cba"];
*/
        Bivar_Function<> fDivide(&divide<double>);
        Tabc.contract(
          1.0, Tabc,"abc", getEnergyDenominator(i,j,k),"abc", 0.0,"abc", fDivide
        );

        // permute a,b,c and i,j,k together into SVabc
        (*SVabc)["abc"]  = getDoublesContribution(i,j,k)["abc"];
        (*SVabc)["abc"] += getDoublesContribution(i,k,j)["acb"];
        (*SVabc)["abc"] += getDoublesContribution(j,i,k)["bac"];
        (*SVabc)["abc"] += getDoublesContribution(j,k,i)["bca"];
        (*SVabc)["abc"] += getDoublesContribution(k,i,j)["cab"];
        (*SVabc)["abc"] += getDoublesContribution(k,j,i)["cba"];

        // contract
        energy[""] += (*SVabc)["abc"] * Tabc["abc"];
      }
    }
  }
  delete DVabc; delete SVabc;

  double eTriples(energy.get_val());
  double eCcsd(getRealArgument("CcsdEnergy"));
  double e(eCcsd + eTriples);
  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << eTriples << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::runInMemory() {
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
  Tensor<> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");

  {
    Tensor<> Zabcijk(6, vvvooo, syms, *Vabij->wrld, "Zabcijk");
    Zabcijk["abcijk"]  = (*Vabci)["bcdk"] * (*Tabij)["adij"];
    Zabcijk["abcijk"] -= (*Vijka)["jklc"] * (*Tabij)["abil"];

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

    Zabcijk["abcijk"]  = Tabcijk["abcijk"];
    Tabcijk["abcijk"] += Zabcijk["bacjik"];
    Tabcijk["abcijk"] += Zabcijk["acbikj"];
    Tabcijk["abcijk"] += Zabcijk["cbakji"];
    Tabcijk["abcijk"] += Zabcijk["cabkij"];
    Tabcijk["abcijk"] += Zabcijk["bcajki"];
  }

  Tensor<> Zai(false, *Tai);
  Zai.set_name("Zai");
  Tensor<> Zabij(false, *Tabij);
  Zabij.set_name("Zabij");

  Zai["ai"]  = (+2.0) * Tabcijk["acdikl"] * (*Vabij)["cdkl"];
  Zai["ai"] += (-1.0) * Tabcijk["acdikl"] * (*Vabij)["dckl"];
  Zai["ai"] += (-2.0) * Tabcijk["acdlki"] * (*Vabij)["cdkl"];
  Zai["ai"] += (+1.0) * Tabcijk["acdlki"] * (*Vabij)["dckl"];

  Zabij["abij"]  = (+2.0) * Tabcijk["acdijk"] * (*Vabci)["cdbk"];
  Zabij["abij"] += (-1.0) * Tabcijk["acdijk"] * (*Vabci)["dcbk"];
  Zabij["abij"] += (-1.0) * Tabcijk["acdkji"] * (*Vabci)["cdbk"];

  Zabij["abij"] += (-2.0) * Tabcijk["abcikl"] * (*Vijka)["kljc"];
  Zabij["abij"] += (+1.0) * Tabcijk["abcikl"] * (*Vijka)["lkjc"];
  Zabij["abij"] += (+1.0) * Tabcijk["abclki"] * (*Vijka)["kljc"];

  Scalar<> energy(*Cc4s::world);
  double e, triplese;
  double ccsde(getRealArgument("CcsdEnergy"));

  energy[""]  = (+2.0) *     Zai["ai"] *     (*Tai)["ai"];
  energy[""] += (+4.0) * Zabij["abij"] * (*Tabij)["abij"];
  energy[""] += (-2.0) * Zabij["abij"] * (*Tabij)["abji"];
  triplese = energy.get_val();
  e = triplese + ccsde;

  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << ccsde << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << triplese << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::runPiecuch() {
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
  Tensor<> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");
  Tensor<> Xabcijk(6, vvvooo, syms, *Vabij->wrld, "Xabcijk");

  Tabcijk["abcijk"]  = (*Vabci)["bcek"] * (*Tabij)["aeij"];
  Tabcijk["abcijk"] -= (*Vijka)["jkmc"] * (*Tabij)["abim"];

  Xabcijk["abcijk"]  = Tabcijk["abcijk"];
  Tabcijk["abcijk"] += Xabcijk["bacjik"];
  Tabcijk["abcijk"] += Xabcijk["acbikj"];
  Tabcijk["abcijk"] += Xabcijk["cbakji"];
  Tabcijk["abcijk"] += Xabcijk["cabkij"];
  Tabcijk["abcijk"] += Xabcijk["bcajki"];

  {
    Tensor<> Zabcijk(6, vvvooo, syms, *Vabij->wrld, "Zabcijk");
    Zabcijk["abcijk"]  = ( 1.0) * (*Tai)["ai"] * (*Vabij)["bcjk"];
    Zabcijk["abcijk"] += ( 1.0) * (*Tai)["bj"] * (*Vabij)["acik"];
    Zabcijk["abcijk"] += ( 1.0) * (*Tai)["ck"] * (*Vabij)["abij"];

    Xabcijk["abcijk"]  = (4.0/3.0) * Zabcijk["abcijk"];
    Xabcijk["abcijk"] +=    (-2.0) * Zabcijk["acbijk"];
    Xabcijk["abcijk"] += (2.0/3.0) * Zabcijk["bcaijk"];

    Xabcijk["abcijk"] += (4.0/3.0) * Tabcijk["abcijk"];
    Xabcijk["abcijk"] +=    (-2.0) * Tabcijk["acbijk"];
    Xabcijk["abcijk"] += (2.0/3.0) * Tabcijk["bcaijk"];

    Zabcijk["abcijk"]  = (*epsi)["i"];
    Zabcijk["abcijk"] += (*epsi)["j"];
    Zabcijk["abcijk"] += (*epsi)["k"];
    Zabcijk["abcijk"] -= (*epsa)["a"];
    Zabcijk["abcijk"] -= (*epsa)["b"];
    Zabcijk["abcijk"] -= (*epsa)["c"];
    Bivar_Function<> fDivide(&divide<double>);
    Tabcijk.contract(
      1.0, Tabcijk,"abcijk", Zabcijk,"abcijk", 0.0,"abcijk", fDivide
    );
  }

  Scalar<> energy(*Cc4s::world);
  double e, triplese;
  double ccsde(getRealArgument("CcsdEnergy"));

  energy[""]  = Xabcijk["abcijk"] * Tabcijk["abcijk"];
  triplese = energy.get_val();
  e = triplese + ccsde;

  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << ccsde << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << triplese << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::runPiecuchFactorizedInMemory() {
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
  Tensor<> DVabcijk(6, vvvooo, syms, *Vabij->wrld, "DVabcijk");
  DVabcijk["abcijk"]  = (*Vabci)["bcdk"] * (*Tabij)["adij"];
  DVabcijk["abcijk"] -= (*Vijka)["jklc"] * (*Tabij)["abil"];

  Tensor<> SVabcijk(6, vvvooo, syms, *Vabij->wrld, "SVabcijk");
  SVabcijk["abcijk"]  = 0.5 * (*Tai)["ai"] * (*Vabij)["bcjk"];

  Tensor<> Tabcijk(6, vvvooo, syms, *Vabij->wrld, "Tabcijk");
  Tabcijk["abcijk"]  = (+8.0) * DVabcijk["abcijk"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["abcikj"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["abcjik"];
  Tabcijk["abcijk"] += (+2.0) * DVabcijk["abcjki"];
  Tabcijk["abcijk"] += (+2.0) * DVabcijk["abckij"];
  Tabcijk["abcijk"] += (-4.0) * DVabcijk["abckji"];
/*
  Tabcijk["abcijk"] += (+8.0) * SVabcijk["abcijk"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["abcikj"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["abcjik"];
  Tabcijk["abcijk"] += (+2.0) * SVabcijk["abcjki"];
  Tabcijk["abcijk"] += (+2.0) * SVabcijk["abckij"];
  Tabcijk["abcijk"] += (-4.0) * SVabcijk["abckji"];
*/
  SVabcijk["abcijk"]  = (*epsi)["i"];
  SVabcijk["abcijk"] += (*epsi)["j"];
  SVabcijk["abcijk"] += (*epsi)["k"];
  SVabcijk["abcijk"] -= (*epsa)["a"];
  SVabcijk["abcijk"] -= (*epsa)["b"];
  SVabcijk["abcijk"] -= (*epsa)["c"];
  Bivar_Function<> fDivide(&divide<double>);
  Tabcijk.contract(
    1.0, Tabcijk,"abcijk", SVabcijk,"abcijk", 0.0,"abcijk", fDivide
  );

  Scalar<> energy(*Cc4s::world);
  energy[""]  = DVabcijk["abcijk"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["bacjik"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["acbikj"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["cbakji"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["cabkij"] * Tabcijk["abcijk"];
  energy[""] += DVabcijk["bcajki"] * Tabcijk["abcijk"];

  double eTriples(energy.get_val());
  double eCcsd(getRealArgument("CcsdEnergy"));
  double e(eCcsd + eTriples);
  LOG(0, "PerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "PerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  LOG(1, "PerturbativeTriples") << "triples=" << eTriples << std::endl;

  setRealArgument("PerturbativeTriplesEnergy", e);
}

void PerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PPPHCoulombIntegrals");

  DryTensor<> *Tai(
    getTensorArgument<double, DryTensor<double>>("CcsdSinglesAmplitudes")
  );
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("CcsdDoublesAmplitudes")
  );

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
  DryTensor<> Tabcijk(6, vvvooo, syms, SOURCE_LOCATION);

  {
    DryTensor<> Zabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  }

  DryTensor<> Zai(*Tai, SOURCE_LOCATION);
  DryTensor<> Zabij(*Tabij, SOURCE_LOCATION);

  DryScalar<> energy();
}


void PerturbativeTriples::dryRunPiecuch() {
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
  DryTensor<> Tabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  DryTensor<> Yabcijk(6, vvvooo, syms, SOURCE_LOCATION);

  {
    DryTensor<> Zabcijk(6, vvvooo, syms, SOURCE_LOCATION);
    DryTensor<> Xabcijk(6, vvvooo, syms, SOURCE_LOCATION);
  }

  DryScalar<> energy();
}
