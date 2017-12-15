#include <algorithms/RpaApxEnergy.hpp>

#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/SlicedCtfTensor.hpp>
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
      [](double eps, complex &d){ d = eps; }
    )
  ) (
    (*epsa)["a"], ctfPain["ain"]
  );
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double eps, complex &d){ d -= eps; }
    )
  ) (
    (*epsi)["i"], ctfPain["ain"]
  );
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double nu, complex &d){ d = 2.0*d / (d*d + nu*nu); }
    )
  ) (
    (*realNun)["n"], ctfPain["ain"]
  );
  CTF::Tensor<complex> ctfConjPain(ctfPain);
  conjugate(ctfConjPain);

  auto Pain(tcc->createTensor(MT::create(ctfPain)));
  auto conjPain(tcc->createTensor(MT::create(ctfConjPain)));

  auto mp2Energy(tcc->createTensor(std::vector<int>(), "mp2Energy"));
  chiVFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chiV");
  tcc->compile(
    (
      // particle/hole bubble propagating forwards
      (*chiVFGn)["FGn"] <<=
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      // particle/hole bubble propagating backwards
      (*chiVFGn)["FGn"] +=
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"],
      // compute Mp2 direct energy for benchmark of frequency grid
      (*mp2Energy)[""] <<= -0.5 / Tau() *
        (*Wn)["n"] * (*chiVFGn)["FGn"] * (*chiVFGn)["GFn"]
    )
  )->execute();

  CTF::Scalar<complex> ctfMp2Energy;
  ctfMp2Energy[""] = mp2Energy->getMachineTensor<MT>()->tensor[""];
  complex mp2E(ctfMp2Energy.get_val());
  LOG(0, "RPA") << "Mp2 direct energy=" << mp2E << std::endl;
  setRealArgument("Mp2Energy", std::real(mp2E));

  diagonalizeChiV();
}

void RpaApxEnergy::diagonalizeChiV() {
  typedef CtfMachineTensor<complex> MT;
  // slice CTF tensor of chiV along n (=3rd) dimension
  SlicedCtfTensor<complex>(chiVFGn->getMachineTensor<MT>()->tensor, {2});
  // TODO: continue here...
}

