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
  auto epsi( getTensorArgument<double>("HoleEigenEnergies") );
  auto epsa( getTensorArgument<double>("ParticleEigenEnergies") );
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  // slice GammaFab,GammaFai from GammaFqr
  auto GammaFqr( getTensorArgument<complex>("CoulombVertex") );
  int NF(GammaFqr->lens[0]);
  int Np(GammaFqr->lens[1]);
  int iStart(0), iEnd(No);
  int aStart(Np-Nv), aEnd(Np);
  int FabStart[] = { 0, aStart, aStart };
  int FabEnd[] = { NF, aEnd, aEnd };
  int FaiStart[] = { 0, aStart, iStart};
  int FaiEnd[] = { NF, aEnd, iEnd };
  auto GammaFab( new Tensor<complex>(GammaFqr->slice(FabStart,FabEnd)) );
  auto GammaFai( new Tensor<complex>(GammaFqr->slice(FaiStart,FaiEnd)) );

  double eTriples(0.0);
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  if (realVabij) {
    eTriples = Calculator<double>(
      getTensorArgument<double>("CcsdSinglesAmplitudes"),
      getTensorArgument<double>("CcsdDoublesAmplitudes"),
      getTensorArgument<double>("PPHHCoulombIntegrals"),
      getTensorArgument<double>("PHHHCoulombIntegrals"),
      GammaFab, GammaFai,
      epsi, epsa
    ).calculate();
    LOG(1, "CcsdPerturbativeTriples") << "triples=" << eTriples << std::endl;
  } else {
    complex complexETriples(
      Calculator<complex>(
        getTensorArgument<complex>("CcsdSinglesAmplitudes"),
        getTensorArgument<complex>("CcsdDoublesAmplitudes"),
        getTensorArgument<complex>("PPHHCoulombIntegrals"),
        getTensorArgument<complex>("PHHHCoulombIntegrals"),
        GammaFab, GammaFai,
        epsi, epsa
      ).calculate()
    );
    eTriples = std::real(complexETriples);
    LOG(1, "CcsdPerturbativeTriples") << "triples=" << complexETriples << std::endl;
  }
  double eCcsd(getRealArgument("CcsdEnergy"));
  LOG(1, "CcsdPerturbativeTriples") << "ccsd=" << eCcsd << std::endl;
  double e(eCcsd + eTriples);
  LOG(0, "CcsdPerturbativeTriples") << "e=" << e << std::endl;
  setRealArgument("CcsdPerturbativeTriplesEnergy", e);
}


template <typename F>
CcsdPerturbativeTriples::Calculator<F>::Calculator(
  Tensor<F> *Tai_, Tensor<F> *Tabij_,
  Tensor<F> *Vabij_, Tensor<F> *Valij_,
  Tensor<complex> *GammaFab_,
  Tensor<complex> *GammaFai_,
  Tensor<double> *epsi_, Tensor<double> *epsa_
):
  // slice Tai in dimension 1, i.e. in i
  Tai( new SlicedCtfTensor<F>(*Tai_, {1}) ),
  // slice Tabij in dimension 2&3, i.e. in i,j
  Tabij( new SlicedCtfTensor<F>(*Tabij_, {2,3}) ),
  // slice Tabil in dimension 2
  Tabil( new SlicedCtfTensor<F>(*Tabij_, {2}) ),
  // slice Vabij in dimension 2&3, i.e. in i,j
  Vabij( new SlicedCtfTensor<F>(*Vabij_, {2,3}) ),
  // slice Valij in dimension 2&3, i.e. in i,j
  Valij( new SlicedCtfTensor<F>(*Valij_, {2,3}) ),
  // build coulomb vertex for the respective type (complex or double)
  Gamma(GammaFab_, GammaFai_)
{
  Tensor<F> unslicedEpsi(1, epsi_->lens, epsi_->sym, *epsi_->wrld, "epsi");
  // convert real tensor epsi_ into tensor of type F unslicedEpsi
  toComplexTensor(*epsi_, unslicedEpsi);
  // slice in dimension 0
  epsi = new SlicedCtfTensor<F>(unslicedEpsi, {0});

  epsa = new Tensor<F>(1, epsa_->lens, epsa_->sym, *epsa_->wrld, "epsa");
  // convert real tensor epsa_ into tensor of type F epsa
  toComplexTensor(*epsa_, *epsa);
}

template <typename F>
CcsdPerturbativeTriples::Calculator<F>::~Calculator() {
  delete Tai; delete Tabij; delete Tabil;
  delete Vabij; delete Valij;
  delete epsi; delete epsa;
}

template <typename F>
void CcsdPerturbativeTriples::Calculator<F>::addDoublesHoleContribution(
  const Map<3> &i, Tensor<F> &DVabc
) {
  DVabc["abc"] -= (*Tabil)({i(0)})["abil"] * (*Valij)({i(2),i(1)})["clkj"];
}

template <typename F>
Tensor<F> &CcsdPerturbativeTriples::Calculator<F>::getSinglesContribution(
  const Map<3> &i
) {
  (*SVabc)["abc"] = 0.5 * (*Tai)({i(0)})["ai"] * (*Vabij)({i(1),i(2)})["bcjk"];
  return *SVabc;
}

template <typename F>
Tensor<F> &CcsdPerturbativeTriples::Calculator<F>::getEnergyDenominator(
  const Map<3> &i
) {
  // reuse SVabc to hold the energy denominator
  (*SVabc)["abc"]  = (*epsi)({i(0)})["i"];
  (*SVabc)["abc"] += (*epsi)({i(1)})["j"];
  (*SVabc)["abc"] += (*epsi)({i(2)})["k"];
  (*SVabc)["abc"] -= (*epsa)["a"];
  (*SVabc)["abc"] -= (*epsa)["b"];
  (*SVabc)["abc"] -= (*epsa)["c"];
  return *SVabc;
}

template <typename F>
F CcsdPerturbativeTriples::Calculator<F>::calculate() {
  int No(epsi->slicedLens[0]);
  int Nv(epsa->lens[0]);
  int vvv[] = { Nv, Nv, Nv };
  int syms[] = { NS, NS, NS };
  // doubles amplitudes contracted with V for current i,j,k
  DVabc = new Tensor<F>(3, vvv, syms, *epsa->wrld, "DVabc");
  // unconnected singles amplitudes and V for current i,j,k
  SVabc = new Tensor<F>(3, vvv, syms, *epsa->wrld, "SVabc");
  // triples amplitudes considering all permutations of a,b,c for given i,j,k
  Tensor<F> Tabc(3, vvv, syms, *epsa->wrld, "Tabc");

  // D.V for all permutations of a,b,c together with i,j,k
  Tensor<F> *piDVabc[Permutation<3>::ORDER];
  for (int p(0); p < Permutation<3>::ORDER; ++p) {
    piDVabc[p] = new Tensor<F>(3, vvv, syms, *epsa->wrld, "piDVabc");
  }
  // spin factors and Fermion sign depending on the number of
  // invariant indices when permuting a,b,c keeping i,j,k fixed.
  // having 2 invariant indices is impossible
  double spinAndFermiFactors[] = { +2.0, -4.0, 0.0, +8.0 };

  // false if permutation Pi leaves current values of indices i,j,k invariant
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
            // permutation p equivalent to previous q for given values of i,j,k
            givesDistinctIndexPermutation[p] = false;
            // use previously calculated permutation
            (*piDVabc[p])["abc"] = (*piDVabc[q])["abc"];
          } else {
            givesDistinctIndexPermutation[p] = true;
            // non-equivalent: calculate for given values of i,j,k
            Gamma.getDoublesParticleContribution(*Tabij, i*pi, *piDVabc[p]);
            addDoublesHoleContribution(i*pi, *piDVabc[p]);
          }
          // aggregate all simultaneous permutations of i,j,k and a,b,c
          (*DVabc)["abc"] += (*piDVabc[p])[("abc"*pi).c_str()];
        }

        // energy denominator is invariant under all permutations
        CTF::Transform<F, F>(
          std::function<void(F, F &)>(
            [](F deltaabij, F &dvabij) {
              dvabij = conj(dvabij / deltaabij);
            }
          )
        ) (
          getEnergyDenominator(i)["abc"], (*DVabc)["abc"]
        );

        for (int p(0); p < Permutation<3>::ORDER; ++p) {
          if (givesDistinctIndexPermutation[p]) {
            // for all permutation pi giving distinct values of i,j,k
            Permutation<3> pi(p);
            Tabc["abc"] = 0.0;
            for (int s(0); s < Permutation<3>::ORDER; ++s) {
              // after pi, permute a,b,c with sigma leaving i,j,k fixed.
              Permutation<3> sigma(s);
              // the spin factors and the Fermion sign only depend on sigma
              double sf(spinAndFermiFactors[sigma.invariantElementsCount()]);
              // get precomputed D.V in
              // permutation sigma*pi of a,b,k and in permutation pi of i,j,k
              Tabc["abc"] += sf * (*piDVabc[p])[("abc"*sigma*pi).c_str()];

              // get S.V in
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

  return energy.get_val();
}


template <int Dummy>
CcsdPerturbativeTriples::CoulombVertex<double,Dummy>::CoulombVertex(
  Tensor<complex> *GammaFab,
  Tensor<complex> *GammaFai
) {
  // split GammaFab into real and imaginary parts
  realGammaFab = new Tensor<double>(
    3, GammaFab->lens, GammaFab->sym, *GammaFab->wrld, "RealGammaFab"
  );
  imagGammaFab = new Tensor<double>(false, *realGammaFab);
  fromComplexTensor(*GammaFab, *realGammaFab, *imagGammaFab);
  // complex GammaFab no longer needed
  delete GammaFab;

  // Split GammaFai into real and imaginary parts and slice along i
  Tensor<double> unslicedRealGammaFai(
    3, GammaFai->lens, GammaFai->sym, *GammaFai->wrld, "RealGammaFai"
  );
  Tensor<double> unslicedImagGammaFai(false, unslicedRealGammaFai);
  fromComplexTensor(*GammaFai, unslicedRealGammaFai, unslicedImagGammaFai);
  // complex GammaFai no longer needed
  delete GammaFai;
  // slice real and imag GammaFai in dimension 2, i.e. in i
  realGammaFai = new SlicedCtfTensor<double>(unslicedRealGammaFai, {2});
  imagGammaFai = new SlicedCtfTensor<double>(unslicedImagGammaFai, {2});
}

template <int Dummy>
CcsdPerturbativeTriples::CoulombVertex<complex,Dummy>::CoulombVertex(
  Tensor<complex> *GammaFab_,
  Tensor<complex> *GammaFai_
):
  conjGammaFab( GammaFab_ ),
  // slice GammaFai in dimension 2, i.e. in i
  GammaFai( new SlicedCtfTensor<complex>(*GammaFai_, {2}) )
{
  conjugate(*conjGammaFab),
  // the unsliced GammaFai_ is no longer needed
  delete GammaFai_;
}

template <int Dummy>
CcsdPerturbativeTriples::CoulombVertex<double,Dummy>::~CoulombVertex() {
  delete realGammaFab; delete imagGammaFab;
  delete realGammaFai; delete imagGammaFai;
}

template <int Dummy>
CcsdPerturbativeTriples::CoulombVertex<complex,Dummy>::~CoulombVertex() {
  delete conjGammaFab;
  delete GammaFai;
}

template <int Dummy>
void
CcsdPerturbativeTriples::CoulombVertex<double,Dummy>::
getDoublesParticleContribution(
  SlicedCtfTensor<double> &Tabij, const Map<3> &i, Tensor<double> &DVabc
) {
  DVabc["abc"] =
    Tabij({i(0),i(1)})["adij"] *
    (*realGammaFab)["Fbd"] * (*realGammaFai)({i(2)})["Fck"];
  DVabc["abc"] +=
    Tabij({i(0),i(1)})["adij"] *
    (*imagGammaFab)["Fbd"] * (*imagGammaFai)({i(2)})["Fck"];
}

template <int Dummy>
void
CcsdPerturbativeTriples::CoulombVertex<complex,Dummy>::
getDoublesParticleContribution(
  SlicedCtfTensor<complex> &Tabij, const Map<3> &i, Tensor<complex> &DVabc
) {
  DVabc["abc"] = Tabij({i(0),i(1)})["adij"] *
    (*conjGammaFab)["Fdb"] * (*GammaFai)({i(2)})["Fck"];
}


// FIXME: dryrun for complex
void CcsdPerturbativeTriples::dryRun() {
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  getTensorArgument<double, DryTensor<double>>("PHHHCoulombIntegrals");

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
