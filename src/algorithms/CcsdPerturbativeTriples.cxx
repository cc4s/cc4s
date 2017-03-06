#include <algorithms/CcsdPerturbativeTriples.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdPerturbativeTriples);

CcsdPerturbativeTriples::CcsdPerturbativeTriples(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CcsdPerturbativeTriples::~CcsdPerturbativeTriples() {
}

Tensor<> &CcsdPerturbativeTriples::getSinglesContribution(int i, int j, int k) {
  Tensor<> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  int TaiStart[] = { 0, i }, TaiEnd[] = { Nv, i+1 };
  int VbcjkStart[] = { 0, 0, j, k }, VbcjkEnd[] = { Nv, Nv, j+1, k+1 };
  (*SVabc)["abc"] =
    0.5 * Tai->slice(TaiStart,TaiEnd)["ai"]
    * Vabij->slice(VbcjkStart,VbcjkEnd)["bcjk"];
  return *SVabc;
}

Tensor<> &CcsdPerturbativeTriples::getDoublesContribution(int i, int j, int k) {
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

Tensor<> &CcsdPerturbativeTriples::getEnergyDenominator(int i, int j, int k) {
  // reuse SVabc to hold the energy denominator
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  // NOTE: due to a bug we first to feed the sliced vector into a scalar
  Scalar<> eps(*epsi->wrld);
  int epsiStart[] = { i }, epsiEnd[] = { i+1 };
  eps[""] = epsi->slice(epsiStart,epsiEnd)["i"];
  (*SVabc)["abc"]  = eps.get_val();
  int epsjStart[] = { j }, epsjEnd[] = { j+1 };
  eps[""] = epsi->slice(epsjStart,epsjEnd)["j"];
  (*SVabc)["abc"] += eps.get_val();
  int epskStart[] = { k }, epskEnd[] = { k+1 };
  eps[""] = epsi->slice(epskStart,epskEnd)["k"];
  (*SVabc)["abc"] += eps.get_val();
  (*SVabc)["abc"] -= (*epsa)["a"];
  (*SVabc)["abc"] -= (*epsa)["b"];
  (*SVabc)["abc"] -= (*epsa)["c"];
  return *SVabc;
}

void CcsdPerturbativeTriples::run() {
  Tensor<>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  No = epsi->lens[0];
  Nv = epsa->lens[0];
  int vvv[] = { Nv, Nv, Nv };
  int syms[] = { NS, NS, NS };
  // doubles amplitudes contracted with V for current i,j,k
  DVabc = new Tensor<>(3, vvv, syms, *epsi->wrld, "DVabc");
  // unconnected singles amplitudes and V for current i,j,k
  SVabc = new Tensor<>(3, vvv, syms, *epsi->wrld, "SVabc");
  // triples amplitudes considering all permutations of a,b,c for given i,j,k
  Tensor<> Tabc(3, vvv, syms, *epsi->wrld, "Tabc");

  Scalar<> energy(*Cc4s::world);
  energy[""] = 0.0;
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
        getSinglesContribution(i,j,k);
        Tabc["abc"] += (+8.0) * (*SVabc)["abc"];
        Tabc["abc"] += (-4.0) * (*SVabc)["acb"];
        Tabc["abc"] += (-4.0) * (*SVabc)["bac"];
        Tabc["abc"] += (+2.0) * (*SVabc)["bca"];
        Tabc["abc"] += (+2.0) * (*SVabc)["cab"];
        Tabc["abc"] += (-4.0) * (*SVabc)["cba"];

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
  LOG(0, "CcsdPerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "CcsdPerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  LOG(1, "CcsdPerturbativeTriples") << "triples=" << eTriples << std::endl;

  setRealArgument("CcsdPerturbativeTriplesEnergy", e);
}

void CcsdPerturbativeTriples::dryRun() {
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
