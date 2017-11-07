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


void CcsdPerturbativeTriples::run() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = Calculator<double>.run();
  } else {
    e = std::real( Calculator<complex>.run() );
  }
  setRealArgument("CcsdPerturbativeTriplesEnergy", e);
}


template <>
void CcsdPerturbativeTriples::Calculator<double>::sliceTensors() {
  // slice Tai in dimension 1, i.e. in i
  Tai = new SlicedCtfTensor<F>(
    *getTensorArgument("CcsdSinglesAmplitudes"), {1}
  );
  // slice Tabij in dimension 2&3, i.e. in i,j
  Tabij = new SlicedCtfTensor<F>(
    *getTensorArgument("CcsdDoublesAmplitudes"), {2,3}
  );
  // slice Tabil in dimension 2
  Tabil = new SlicedCtfTensor<F>(
    *getTensorArgument("CcsdDoublesAmplitudes"), {2}
  );
  // slice Vabij in dimension 2&3, i.e. in i,j
  Vabij = new SlicedCtfTensor<F>(
    *getTensorArgument("PPHHCoulombIntegrals"), {2,3}
  );
  // slice Vijla in dimension 0&1, i.e. in i,j
  Vijla = new SlicedCtfTensor<F>(
    *getTensorArgument("HHHPCoulombIntegrals"), {0,1}
  );

  Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));
  // Allocate and compute GammaFab,GammaFai from GammaFqr
  int NF(GammaFqr->lens[0]);
  int Np(GammaFqr->lens[1]);
  int iStart(0), iEnd(No);
  int aStart(Np-Nv), aEnd(Np);
  int FaiStart[] = { 0, aStart, iStart};
  int FaiEnd[] = { NF, aEnd, iEnd };
  int FabStart[] = { 0, aStart, aStart };
  int FabEnd[] = { NF, aEnd, aEnd };
  Tensor<complex> GammaFai(GammaFqr->slice(FaiStart,FaiEnd));
  Tensor<complex> GammaFab(GammaFqr->slice(FabStart,FabEnd));
  // Split GammaFai,GammaFab into real and imaginary parts
  Tensor<> unslicedRealGammaFai(
    3, GammaFai.lens, GammaFai.sym, *GammaFai.wrld, "RealGammaFai"
  );
  Tensor<> unslicedImagGammaFai(unslicedRealGammaFai);
  fromComplexTensor(GammaFai, unslicedRealGammaFai, unslicedImagGammaFai);
  // slice real and imag GammaFai in dimension 2, i.e. in i
  realGammaFai = new SlicedCtfTensor<>(unslicedRealGammaFai, {2});
  imagGammaFai = new SlicedCtfTensor<>(unslicedImagGammaFai, {2});

  realGammaFab = new Tensor<>(
    3, GammaFab.lens, GammaFab.sym, *GammaFab.wrld, "RealGammaFab"
  );
  imagGammaFab = new Tensor<>(*realGammaFab);
  fromComplexTensor(GammaFab, *realGammaFab, *imagGammaFab);
}

template <>
void CcsdPerturbativeTriples::Calculator<complex>::sliceTensors() {
  // slice Tai in dimension 1, i.e. in i
  Tai = new SlicedCtfTensor<complex>(
    *getTensorArgument("CcsdSinglesAmplitudes"), {1}
  );
  // slice Tabij in dimension 2&3, i.e. in i,j
  Tabij = new SlicedCtfTensor<complex>(
    *getTensorArgument("CcsdDoublesAmplitudes"), {2,3}
  );
  // slice Tabil in dimension 2
  Tabil = new SlicedCtfTensor<complex>(
    *getTensorArgument("CcsdDoublesAmplitudes"), {2}
  );
  // slice Vabij in dimension 2&3, i.e. in i,j
  Vabij = new SlicedCtfTensor<complex>(
    *getTensorArgument("PPHHCoulombIntegrals"), {2,3}
  );
  // slice Vijla in dimension 0&1, i.e. in i,j
  Vijla = new SlicedCtfTensor<complex>(
    *getTensorArgument("HHHPCoulombIntegrals"), {0,1}
  );

  Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));
  // Allocate and compute GammaFab,GammaFai from GammaFqr
  int NF(GammaFqr->lens[0]);
  int Np(GammaFqr->lens[1]);
  int iStart(0), iEnd(No);
  int aStart(Np-Nv), aEnd(Np);
  int FaiStart[] = { 0, aStart, iStart};
  int FaiEnd[] = { NF, aEnd, iEnd };
  int FabStart[] = { 0, aStart, aStart };
  int FabEnd[] = { NF, aEnd, aEnd };
  Tensor<complex> GammaFai(GammaFqr->slice(FaiStart,FaiEnd));
  Tensor<complex> GammaFab(GammaFqr->slice(FabStart,FabEnd));
  // slice GammaFai in dimension 2, i.e. in i
  GammaFai = new SlicedCtfTensor<complex>(GammaFai, {2});
}

template <typename F>
Tensor<F> &CcsdPerturbativeTriples::getSinglesContribution(const Map<3> &i) {
  (*SVabc)["abc"] = 0.5 * (*Tai)({i(0)})["ai"] * (*Vabij)({i(1),i(2)})["bcjk"];
  return *SVabc;
}

Tensor<double> &CcsdPerturbativeTriples::getDoublesContribution(const Map<3> &i) {
  (*SVabc)["abc"] =
    (*Tabij)({i(0),i(1)})["adij"] *
    (*realGammaFab)["Fbd"] * (*realGammaFai)({i(2)})["Fck"];
  (*SVabc)["abc"] +=
    (*Tabij)({i(0),i(1)})["adij"] *
    (*imagGammaFab)["Fbd"] * (*imagGammaFai)({i(2)})["Fck"];

  (*SVabc)["abc"] -= (*Tabil)({i(0)})["abil"] * (*Vijla)({i(1),i(2)})["jklc"];
  return *SVabc;
}

Tensor<complex> &CcsdPerturbativeTriples::getDoublesContribution(const Map<3> &i) {
  // TODO: GammaFab["Fbd"] should be conjugate transposed
  (*SVabc)["abc"] =
    (*Tabij)({i(0),i(1)})["adij"] * (*GammaFab)["Fbd"] * (*GammaFai)({i(2)})["Fck"];

  // TODO: Vijla["jklc"] should be in complex Vaijk["clkj"]
  (*SVabc)["abc"] -= (*Tabil)({i(0)})["abil"] * (*Vijla)({i(1),i(2)})["jklc"];
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

template <typename F>
F CcsdPerturbativeTriples::run() {
  Tensor<>  *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<>  *epsa(getTensorArgument("ParticleEigenEnergies"));
  No = epsi->lens[0];
  Nv = epsa->lens[0];
  int vvv[] = { Nv, Nv, Nv };
  int syms[] = { NS, NS, NS };
  // doubles amplitudes contracted with V for current i,j,k
  DVabc = new Tensor<F>(3, vvv, syms, *epsi->wrld, "DVabc");
  // unconnected singles amplitudes and V for current i,j,k
  SVabc = new Tensor<F>(3, vvv, syms, *epsi->wrld, "SVabc");
  // triples amplitudes considering all permutations of a,b,c for given i,j,k
  Tensor<F> Tabc(3, vvv, syms, *epsi->wrld, "Tabc");

  sliceTensors();

  // D.V for all permutations of a,b,c together with i,j,k
  Tensor<> *piDVabc[Permutation<3>::ORDER];
  for (int p(0); p < Permutation<3>::ORDER; ++p) {
    piDVabc[p] = new Tensor<F>(3, vvv, syms, *epsi->wrld, "piDVabc");
  }
  // spin factors and Fermion sign depending on the number of
  // invariant indices when permuting a,b,c keeping i,j,k fixed.
  // having 2 invariant indices is impossible
  double spinAndFermiFactors[] = { +2.0, -4.0, 0.0, +8.0 };

  // true if the permutation Pi leaves the current indices i,j,k invariant
  bool givesDistinctIndexPermutation[Permutation<3>::ORDER];

  Scalar<F> energy(*Cc4s::world);
  energy[""] = 0.0;
  // indices i,j,k as map with 3 elements i(0),...,i(2)
  Map<3> i;
  // go through all N distinct orders 0 <= i(0) <= i(1) <= i(2) < No
  int n(0); const int N(No*(No+1)*(No+2) / 6);
  Time startTime(Time::getCurrentRealTime());
  for (i(0) = 0; i(0) < No; ++i(0)) {
    for (i(1) = i(0); i(1) < No; ++i(1)) {
      for (i(2) = i(1); i(2) < No; ++i(2)) {
        // get D.V in all permuations of i,j,k
        // and build sum over all permutations of i,j,k together with a,b,c
        (*DVabc)["abc"] = 0.0;
        for (int p(0); p < Permutation<3>::ORDER; ++p) {
          Permutation<3> pi(p);
          int q;
          // check if previsous permutation q permutes current i,j,k same as pi
          for (q = 0; q < p; ++q) if (i*Permutation<3>(q) == i*pi) break;
          if (q < p) {
            // permutation p equivalent to a previous q for the given i,j,k
            givesDistinctIndexPermutation[p] = false;
            // use previously calculated permutation
            (*piDVabc[p])["abc"] = (*piDVabc[q])["abc"];
          } else {
            givesDistinctIndexPermutation[p] = true;
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
          if (givesDistinctIndexPermutation[p]) {
            // all distinct permutation pi of i,j,k together with a,b,c
            Permutation<3> pi(p);
            Tabc["abc"] = 0.0;
            for (int s(0); s < Permutation<3>::ORDER; ++s) {
              // after pi, permute a,b,c with sigma leaving i,j,k fixed.
              Permutation<3> sigma(s);
              // the spin factors and the Fermion sign only depend on sigma
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
            energy[""] += (*DVabc)["abc"] * Tabc["abc"];
          }
        }
        ++n;
        LOG(1, "CcsdPerturbativeTriples") << n << "/" << N <<
          " distinct indices calculated, ETA=" <<
          (Time::getCurrentRealTime()-startTime)*(static_cast<double>(N)/n-1) <<
          " s" << std::endl;
      }
    }
  }

  delete DVabc; delete SVabc;
  for (int p(0); p < Permutation<3>::ORDER; ++p) delete piDVabc[p];
  delete realGammaFab; delete imagGammaFab;
  delete Tai; delete Tabij; delete Tabil;
  delete Vabij; delete Vijla; delete realGammaFai; delete imagGammaFai;

  F eTriples(energy.get_val());
  F eCcsd(getRealArgument("CcsdEnergy"));
  F e(eCcsd + eTriples);
  LOG(0, "CcsdPerturbativeTriples") << "e=" << e << std::endl;
  LOG(1, "CcsdPerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  LOG(1, "CcsdPerturbativeTriples") << "triples=" << eTriples << std::endl;

  return e;
}

void CcsdPerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("HHHPCoulombIntegrals");

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
