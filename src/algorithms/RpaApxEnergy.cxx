#include <algorithms/RpaApxEnergy.hpp>

#include <Cc4s.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/SlicedCtfTensor.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>

#include <ctf.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(RpaApxEnergy);

RpaApxEnergy::RpaApxEnergy(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

RpaApxEnergy::~RpaApxEnergy() {
}

void RpaApxEnergy::run() {
  // create tcc infrastructure
  typedef CtfMachineTensor<complex> MT;
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(tcc::Tcc<complex>::create(machineTensorFactory));

  // read the Coulomb vertex GammaGqr
  auto GammaFqr(getTensorArgument<complex>("CoulombVertex"));

  // read the Particle/Hole Eigenenergies
  auto epsi(getTensorArgument<>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // get index ranges for particles and holes
  int NF(GammaFqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(GammaFqr->lens[1]);

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int FiaStart[] = {0, iStart,aStart};
  int FiaEnd[]   = {NF,iEnd,  aEnd};
  int FaiStart[] = {0, aStart,iStart};
  int FaiEnd[]   = {NF,aEnd,  iEnd};
  auto ctfGammaFia(GammaFqr->slice(FiaStart, FiaEnd));
  auto ctfGammaFai(GammaFqr->slice(FaiStart, FaiEnd));
  auto ctfConjGammaFia(ctfGammaFia);
  conjugate(ctfConjGammaFia); 
  auto ctfConjGammaFai(ctfGammaFai);
  conjugate(ctfConjGammaFai);
  auto GammaFia(tcc->createTensor(MT::create(ctfGammaFia)));
  auto GammaFai(tcc->createTensor(MT::create(ctfGammaFai)));
  auto conjGammaFia(tcc->createTensor(MT::create(ctfConjGammaFia)));
  auto conjGammaFai(tcc->createTensor(MT::create(ctfConjGammaFai)));

  auto realNun(getTensorArgument("ImaginaryFrequencyPoints"));
  int Nn(realNun->lens[0]);

  auto realWn(getTensorArgument("ImaginaryFrequencyWeights"));
  CTF::Tensor<complex> complexWn(1, &Nn, *realWn->wrld);
  toComplexTensor(*realWn, complexWn);
  auto Wn(tcc->createTensor(MT::create(complexWn)));

  CTF::Tensor<complex> ctfPain(3, std::vector<int>({Nv,No,Nn}).data());
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double eps, complex &d) { d = eps; }
    )
  ) (
    (*epsa)["a"], ctfPain["ain"]
  );
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double eps, complex &d) { d -= eps; }
    )
  ) (
    (*epsi)["i"], ctfPain["ain"]
  );
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      // particle/hole propagator for positive and negative nu
      [](double nu, complex &d) { d = 1.0 / (d- complex(0,1)*nu); }
    )
  ) (
    (*realNun)["n"], ctfPain["ain"]
  );
  CTF::Tensor<complex> ctfConjPain(ctfPain);
  conjugate(ctfConjPain);

  auto Pain(tcc->createTensor(MT::create(ctfPain)));
  auto conjPain(tcc->createTensor(MT::create(ctfConjPain)));

  auto mp2Direct(tcc->createTensor(std::vector<int>(), "mp2Direct"));
  auto mp2Exchange(tcc->createTensor(std::vector<int>(), "mp2Exchange"));
  chiVFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chiV");
  PxVFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "PxV");
  double spins(2.0);
  tcc->compile(
    (
      // bubble with half V on both ends:
      // particle/hole bubble propagating forwards
      (*chiVFGn)["FGn"] <<= spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      // particle/hole bubble propagating backwards, positive nu
      (*chiVFGn)["FGn"] += spins *
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"],

      // adjacent pairs exchanged
      (*PxVFGn)["FGn"] <<= spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Haj"] * (*Pain)["ain"] *
        (*GammaFia)["Hib"] * (*conjGammaFia)["Gjb"] * (*conjPain)["bjn"],

      // compute Mp2 energy for benchmark of frequency grid
      // 2 fold rotational and 2 fold mirror symmetry, 1/Pi from +nu and -nu
      (*mp2Direct)[""] <<= -0.25 / Pi() *
        (*Wn)["n"] * (*chiVFGn)["FGn"] * (*chiVFGn)["GFn"],
      // 2 fold mirror symmetry only, 1/Pi from +nu and -nu
      (*mp2Exchange)[""] <<= +0.5 / Pi() * (*Wn)["n"] * (*PxVFGn)["FFn"]
    )
  )->execute();

  CTF::Scalar<complex> ctfMp2Energy;
  ctfMp2Energy[""] = mp2Direct->getMachineTensor<MT>()->tensor[""];
  complex mp2D(ctfMp2Energy.get_val());
  ctfMp2Energy[""] = mp2Exchange->getMachineTensor<MT>()->tensor[""];
  complex mp2X(ctfMp2Energy.get_val());
  LOG(0, "RPA") << "Mp2 direct energy=" << mp2D << std::endl;
  LOG(0, "RPA") << "Mp2 exchange energy=" << mp2X << std::endl;
  setRealArgument("Mp2Energy", std::real(mp2D+mp2X));

  diagonalizeChiV();
}

void RpaApxEnergy::diagonalizeChiV() {
  typedef CtfMachineTensor<complex> MT;
  // get weights for integration
  auto wn( getTensorArgument("ImaginaryFrequencyWeights") );
  std::vector<double> weights(wn->lens[0]);
  wn->read_all(weights.data(), true);

  // slice CTF tensor of chiV and PxV along n (=3rd) dimension
  SlicedCtfTensor<complex> slicedChiVFGn(
    chiVFGn->getMachineTensor<MT>()->tensor, {2}
  );
  SlicedCtfTensor<complex> slicedPxVFGn(
    PxVFGn->getMachineTensor<MT>()->tensor, {2}
  );
  BlacsWorld world(Cc4s::world->rank, Cc4s::world->np);
  complex rpa(0), apx(0);
  for (int n(0); n < slicedChiVFGn.slicedLens[0]; ++n) {
    auto chiVFG( &slicedChiVFGn({n}) );
    auto PxVFG( &slicedPxVFGn({n}) );
    int scaLens[2] = { chiVFG->lens[0], chiVFG->lens[1] };
    auto scaChiVFG( NEW(ScaLapackMatrix<complex>, *chiVFG, scaLens, &world) );
    auto scaU( NEW( ScaLapackMatrix<complex>, *scaChiVFG) );
    ScaLapackHermitianEigenSystemDc<complex> eigenSystem(scaChiVFG, scaU);
    std::vector<double> lambdas(chiVFG->lens[0]);
    eigenSystem.solve(lambdas.data());

    // write diagonalizaing transformation to scliced ChiVFL
    scaU->write(*chiVFG);
    auto conjChiVFG(*chiVFG);
    conjugate(conjChiVFG);

    // write eigenvalues to CTF vector
    CTF::Vector<double> lambdaL(lambdas.size());
    int localNF(lambdaL.wrld->rank == 0 ? lambdas.size() : 0);
    std::vector<int64_t> lambdaIndices(localNF);
    for (int64_t i(0); i < localNF; ++i) { lambdaIndices[i] = i; }
    lambdaL.write(localNF, lambdaIndices.data(), lambdas.data());

    CTF::Vector<complex> LogChiVL(lambdas.size());
    CTF::Transform<double, complex>(
      std::function<void(double, complex &)>(
        [](double chiV, complex &logChiV) {
          logChiV = chiV < 1 ? std::log(1-chiV) + chiV : -chiV*chiV/2;
        }
      )
    ) (
      lambdaL["L"], LogChiVL["L"]
    );

    CTF::Vector<complex> InvChiVL(lambdas.size());
    CTF::Transform<double, complex>(
      std::function<void(double, complex &)>(
        [](double chiV, complex &invChiV) {
          invChiV = chiV < 1 ? 1 / (1-chiV) : 1;
        }
      )
    ) (
      lambdaL["L"], InvChiVL["L"]
    );

/*
    double en(0);
    for (int L(0); L < chiVFG->lens[0]; ++L) {
      // TODO: what is to be done with chiV eigenvalues >= 1?
      if (lambdas[L] < 1) {
        en += std::log(1 - lambdas[L]) + lambdas[L];
      } else {
        LOG(0,"RPA") << "WARNING: chiV(n=" << n << ") > 1, taking MP2 value instead." << std::endl;
        en += -lambdas[L]*lambdas[L]/2;
      }
    }
*/
    CTF::Scalar<complex> e;
    e[""] = LogChiVL["L"];
    rpa += weights[n] * e.get_val();
    e[""] = (*PxVFG)["FG"] * conjChiVFG["GL"] * InvChiVL["L"] * (*chiVFG)["FL"];
    apx += weights[n] * e.get_val();
  }
  // 2 fold mirror symmetry, 1/Pi from +nu and -nu
  rpa *= 0.5 / Pi();
  apx *= 0.5 / Pi();
  LOG(1, "RPA") << "rpa=" << rpa << std::endl;
  LOG(1, "RPA") << "apx=" << apx << std::endl;
  setRealArgument("RpaEnergy", std::real(rpa));
//  setRealArgument("ApxEnergy", std::real(apx));
}

