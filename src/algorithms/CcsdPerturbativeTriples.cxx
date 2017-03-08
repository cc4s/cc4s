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

template <int N>
std::string operator *(const std::string &s, const Permutation<N> pi) {
  Assert(s.length() == N, "Only strings of length N can be permuted.");
  std::string t(N);
  for (int i(0); i < N; ++i) t[i] = s[pi(i)];
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

  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  // Allocate and compute GammaGbd,GammaGck from GammaGqr
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);
  int     iStart(k), iEnd(k+1);
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

  int TadijStart[] = { 0, 0, i, j }, TadijEnd[] = { Nv, Nv, i+1, j+1 };

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

  int TabilStart[] = { 0, 0, i, 0 }, TabilEnd[] = { Nv, Nv, i+1, No };
  int VjklcStart[] = { j, k, 0, 0 }, VjklcEnd[] = { j+1, k+1, No, Nv };
  (*SVabc)["abc"] -=
    Tabij->slice(TabilStart,TabilEnd)["abil"]
    * Vijka->slice(VjklcStart,VjklcEnd)["jklc"];

  return *SVabc;
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

  const int perms[6][3] = {
    {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
  };
  const char *permIndices[6] = {
    "abc", "acb", "bac", "bca", "cab", "cba"
  };
  const double permParticleFactors[6] = {
    +8.0, -4.0, -4.0, +2.0, +2.0, -4.0
  };
  Tensor<> *permDVabc[6];
  for (int p(0); p < 6; ++p) {
    permDVabc[p] = new Tensor<>(3, vvv, syms, *epsi->wrld, "permDVabc");
  }
  bool symmetryPerm[6];

  Scalar<> energy(*Cc4s::world);
  energy[""] = 0.0;
  int i[3];
  // go through all distinct orders 0 <= i[0] <= i[1] <= i[2] < No
  for (i[0] = 0; i[0] < No; ++i[0]) {
    for (i[1] = i[0]; i[1] < No; ++i[1]) {
      for (i[2] = i[1]; i[2] < No; ++i[2]) {
        // get DV in all permuations of i,j,k
        // and build sum over all permutations of i,j,k together with a,b,c
        (*DVabc)["abc"] = 0.0;
        for (int p(0); p < 6; ++p) {
          int q;
          for (q = 0; q < p; ++q) {
            int d;
            for (d = 0; d < 3; ++d) {
              if (i[perms[p][d]] != i[perms[q][d]]) break;
            }
            if (d == 3) break;
          }
          if (q < p) {
            // permutation p equivalent to a previous q for the given i,j,k
            symmetryPerm[p] = true;
            (*permDVabc[p])["abc"] = (*permDVabc[q])["abc"];
          } else {
            symmetryPerm[p] = false;
            // non-equivalent: calculate for given i,j,k
            (*permDVabc[p])["abc"] = getDoublesContribution(
              i[perms[p][0]], i[perms[p][1]], i[perms[p][2]]
            )["abc"];
          }
          (*DVabc)["abc"] += (*permDVabc[p])[permIndices[p]];
        }

        Bivar_Function<> fDivide(&divide<double>);
        DVabc->contract(
          1.0, *DVabc,"abc",
          getEnergyDenominator(i[0],i[1],i[2]),"abc",
          0.0,"abc", fDivide
        );

        for (int p(0); p < 6; ++p) {
          if (!symmetryPerm[p]) {
            Tabc["abc"] = 0.0;
            for (int q(0); q < 6; ++q) {
              // get index names of permutation q.p
              char indices[4];
              for (int d(0); d < 3; ++d) {
                indices[d] = permIndices[p][perms[q][d]];
              }
              indices[3] = 0;
              // get D.V in permutation q.p of a,b,k and permutation p of i,j,k
              Tabc["abc"] += permParticleFactors[q] * (*permDVabc[p])[indices];

              // get S.V in permutation q.p of a,b,k and permutation p of i,j,k
              Tabc["abc"] += permParticleFactors[q] * getSinglesContribution(
                i[perms[p][0]], i[perms[p][1]], i[perms[p][2]]
              )[indices];
            }
            // contract
            Scalar<> contribution(*Cc4s::world);
            contribution[""] += (*DVabc)["abc"] * Tabc["abc"];
            LOG(1, "CcsdPerturbativeTriples") <<
              "p=" << p << ", e[" << i[perms[p][0]] << "," << i[perms[p][1]] << "," << i[perms[p][2]] << "]=" <<
              contribution.get_val() << std::endl;
            energy[""] += contribution[""];
          }
        }
      }
    }
  }
  delete DVabc; delete SVabc;
  for (int p(0); p < 6; ++p) delete permDVabc[p];

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
