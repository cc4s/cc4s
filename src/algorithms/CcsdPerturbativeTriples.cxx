#include <algorithms/CcsdPerturbativeTriples.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/Permutation.hpp>
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

namespace cc4s {
  template <int N>
  inline std::string operator *(
    const std::string &s, const cc4s::Permutation<N> &pi
  ) {
    std::string sPi(s);
    for (int i(0); i < N; ++i) sPi[i] = s[pi(i)];
    return sPi;
  }
}

Tensor<> &CcsdPerturbativeTriples::getSinglesContribution(const Map<3> &i) {
  Tensor<> *Tai(getTensorArgument("CcsdSinglesAmplitudes"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  int TaiStart[] = { 0, i(0) }, TaiEnd[] = { Nv, i(0)+1 };
  int VbcjkStart[] = { 0,0, i(1),i(2) }, VbcjkEnd[] = { Nv,Nv, i(1)+1,i(2)+1 };
  (*SVabc)["abc"] =
    0.5 * Tai->slice(TaiStart,TaiEnd)["ai"]
    * Vabij->slice(VbcjkStart,VbcjkEnd)["bcjk"];
  return *SVabc;
}

Tensor<> &CcsdPerturbativeTriples::getDoublesContribution(const Map<3> &i) {
  Tensor<> *Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));

  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  // Allocate and compute GammaGbd,GammaGck from GammaGqr
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);
  int iStart(i(2)), iEnd(i(2)+1);
  int aStart(Np-Nv),  aEnd(Np);
  int GckStart[] = {0 ,aStart,iStart};
  int   GckEnd[] = {NG,aEnd,  iEnd  };
  int GbdStart[] = {0 ,aStart,aStart};
  int   GbdEnd[] = {NG,aEnd,  aEnd  };
  Tensor<complex> GammaGck(GammaGqr->slice(GckStart,GckEnd));
  Tensor<complex> GammaGbd(GammaGqr->slice(GbdStart,GbdEnd));
  // Split GammaGbd,GammaGck into real and imaginary parts
  Tensor<> realGammaGck(
    3, GammaGck.lens, GammaGck.sym, *GammaGck.wrld, "RealGammaGck"
  );
  Tensor<> imagGammaGck(
    3, GammaGck.lens, GammaGck.sym, *GammaGck.wrld, "ImagGammaGck"
  );
  fromComplexTensor(GammaGck, realGammaGck, imagGammaGck);

  Tensor<> realGammaGbd(
    3, GammaGbd.lens, GammaGbd.sym, *GammaGbd.wrld, "RealGammaGbd"
  );
  Tensor<> imagGammaGbd(
    3, GammaGbd.lens, GammaGbd.sym, *GammaGbd.wrld, "ImagGammaGbd"
  );
  fromComplexTensor(GammaGbd, realGammaGbd, imagGammaGbd);

  int TadijStart[] = { 0,0, i(0),i(1) }, TadijEnd[] = { Nv,Nv, i(0)+1,i(1)+1 };

  (*SVabc)["abc"]  =
    Tabij->slice(TadijStart,TadijEnd)["adij"]
    * realGammaGbd["Gbd"] * realGammaGck["Gck"];
  (*SVabc)["abc"] +=
    Tabij->slice(TadijStart,TadijEnd)["adij"]
    * imagGammaGbd["Gbd"] * imagGammaGck["Gck"];;

  /*
  Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  int VbcdkStart[] = { 0, 0, 0, k }, VbcdkEnd[] = { Nv, Nv, Nv, k+1 };
  (*SVabc)["abc"] =
    Tabij->slice(TadijStart,TadijEnd)["adij"]
    * Vabci->slice(VbcdkStart,VbcdkEnd)["bcdk"];
    */

  int TabilStart[] = { 0, 0, i(0), 0 }, TabilEnd[] = { Nv, Nv, i(0)+1, No };
  int VjklcStart[] = { i(1),i(2), 0,0 }, VjklcEnd[] = { i(1)+1,i(2)+1, No,Nv };
  (*SVabc)["abc"] -=
    Tabij->slice(TabilStart,TabilEnd)["abil"]
    * Vijka->slice(VjklcStart,VjklcEnd)["jklc"];

  return *SVabc;
}

Tensor<> &CcsdPerturbativeTriples::getEnergyDenominator(const Map<3> &i) {
  // reuse SVabc to hold the energy denominator
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  // NOTE: due to a bug we first to feed the sliced vector into a scalar
  Scalar<> eps(*epsi->wrld);
  int epsiStart[] = { i(0) }, epsiEnd[] = { i(0)+1 };
  eps[""] = epsi->slice(epsiStart,epsiEnd)["i"];
  (*SVabc)["abc"]  = eps.get_val();
  int epsjStart[] = { i(1) }, epsjEnd[] = { i(1)+1 };
  eps[""] = epsi->slice(epsjStart,epsjEnd)["j"];
  (*SVabc)["abc"] += eps.get_val();
  int epskStart[] = { i(2) }, epskEnd[] = { i(2)+1 };
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

  Tensor<> *piDVabc[Permutation<3>::ORDER];
  for (int p(0); p < Permutation<3>::ORDER; ++p) {
    piDVabc[p] = new Tensor<>(3, vvv, syms, *epsi->wrld, "piDVabc");
  }
  // spin factors and Fermion sign depending on the number of
  // invariant indices between the when permuting a,b,c keeping i,j,k fixed.
  // having 2 invariant indices is impossible
  double spinAndFermiFactors[] = { +2.0, -4.0, 0.0, +8.0 };

  // true if the permutation Pi leaves the current indices i,j,k invariant
  bool ijkInvariantUnderPi[Permutation<3>::ORDER];

  Scalar<> energy(*Cc4s::world);
  energy[""] = 0.0;
  // indices i,j,k as map with 3 elements i(0),...,i(2)
  Map<3> i;
  // go through all distinct orders 0 <= i(0) <= i(1) <= i(2) < No
  for (i(0) = 0; i(0) < No; ++i(0)) {
    for (i(1) = i(0); i(1) < No; ++i(1)) {
      for (i(2) = i(1); i(2) < No; ++i(2)) {
        // get DV in all permuations of i,j,k
        // and build sum over all permutations of i,j,k together with a,b,c
        (*DVabc)["abc"] = 0.0;
        for (int p(0); p < Permutation<3>::ORDER; ++p) {
          Permutation<3> pi(p);
          int q;
          // check whether i after previsous permutation q leaves i invariant
          for (q = 0; q < p; ++q) if (i*Permutation<3>(q) == i*pi) break;
          if (q < p) {
            // permutation p equivalent to a previous q for the given i,j,k
            ijkInvariantUnderPi[p] = true;
            // use previously calculated permutation
            (*piDVabc[p])["abc"] = (*piDVabc[q])["abc"];
          } else {
            ijkInvariantUnderPi[p] = false;
            // non-equivalent: calculate for given i,j,k
            (*piDVabc[p])["abc"] = getDoublesContribution(i*pi)["abc"];
          }
          // aggregate all simultaneous permutations of i,j,k and a,b,c
          (*DVabc)["abc"] += (*piDVabc[p])[("abc"*pi).c_str()];
        }

        // energy denominator is invariant under all permutations
        Bivar_Function<> fDivide(&divide<double>);
        DVabc->contract(
          1.0, *DVabc,"abc", getEnergyDenominator(i),"abc", 0.0,"abc", fDivide
        );

        for (int p(0); p < Permutation<3>::ORDER; ++p) {
          if (!ijkInvariantUnderPi[p]) {
            // go through all permutation pi of i,j,k together with a,b,c
            Permutation<3> pi(p);
            Tabc["abc"] = 0.0;
            for (int s(0); s < Permutation<3>::ORDER; ++s) {
              // after pi, permute a,b,c with sigma leaving i,j,k fixed.
              Permutation<3> sigma(s);
              // the spin factors and the Fermion sign only depends on sigma
              double sf(spinAndFermiFactors[sigma.invariantElementsCount()]);
              // get D.V in
              // permutation sigma*pi of a,b,k and in permutation pi of i,j,k
              Tabc["abc"] += sf * (*piDVabc[p])[("abc"*sigma*pi).c_str()];

              // get D.V in
              // permutation sigma*pi of a,b,k and in permutation pi of i,j,k
              Tabc["abc"] +=
                sf * getSinglesContribution(i*pi)[("abc"*sigma*pi).c_str()];
            }
            // contract
            Scalar<> contribution(*Cc4s::world);
            contribution[""] += (*DVabc)["abc"] * Tabc["abc"];
            LOG(1, "CcsdPerturbativeTriples") <<
              "pi=" << pi << ", e" << i*pi << "=" <<
              contribution.get_val() << std::endl;
            energy[""] += contribution[""];
          }
        }
      }
    }
  }
  delete DVabc; delete SVabc;
  for (int p(0); p < Permutation<3>::ORDER; ++p) delete piDVabc[p];

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
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int vvv[] = { Nv, Nv, Nv };
  int   syms[] = { NS, NS, NS };
  DryTensor<> SVabcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DVabcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV0abcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV1abcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV2abcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV3abcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV4abcijk(3, vvv, syms, SOURCE_LOCATION);
  DryTensor<> DV5abcijk(3, vvv, syms, SOURCE_LOCATION);

  {
    DryTensor<> Tabcijk(3, vvv, syms, SOURCE_LOCATION);
  }

  DryTensor<> Zai(*Tai, SOURCE_LOCATION);
  DryTensor<> Zabij(*Tabij, SOURCE_LOCATION);

  DryScalar<> energy();
}
