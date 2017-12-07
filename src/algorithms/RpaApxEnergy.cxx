#include <algorithms/RpaApxEnergy.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
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
  typedef CtfMachineTensor<complex> MachineTensor;
  auto machineTensorFactory(MachineTensor::Factory::create());
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
  auto GammaFia(tcc->createTensor(MachineTensor::create(ctfGammaFia)));
  auto GammaFai(tcc->createTensor(MachineTensor::create(ctfGammaFai)));
  auto conjGammaFia(tcc->createTensor(MachineTensor::create(ctfConjGammaFia)));
  auto conjGammaFai(tcc->createTensor(MachineTensor::create(ctfConjGammaFai)));

  auto realNun(getTensorArgument("ImaginaryFrequencyPoints"));
  int Nn(realNun->lens[0]);

  auto realWn(getTensorArgument("ImaginaryFrequencyWeights"));
  CTF::Tensor<complex> complexWn(1, &Nn, *realWn->wrld);
  toComplexTensor(*realWn, complexWn);

  CTF::Tensor<complex> Dain(3, std::vector<int>({Nv,No,Nn}).data());
  Dain["ain"]  = (*epsa)["a"];
  Dain["ain"] -= (*epsa)["i"];
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double nu, complex &d){ d = 1.0 / (d - complex(0,1)*nu); }
    )
  ) (
    Dain["ain"]
  );
  CTF::Tensor<complex> conjDain(Dain);
  conjugate(conjDain);

  auto Pain(tcc->createTensor(MachineTensor::create(Dain)));
  auto conjPain(tcc->createTensor(MachineTensor::create(conjDain)));

  auto energy(tcc->createTensor(std::vector<int>(), "energy"));
  auto chiVFGn(tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chiV"));
  tcc->compile(
    (
      (*chiVFGn)["FGn"] <<=
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      (*chiVFGn)["FGn"] +=
        (*GammaFai)["Fia"] * (*conjGammaFai)["Gia"] * (*conjPain)["ain"]
    )
  )->execute();

  CTF::Scalar<complex> ctfEnergy;
  ctfEnergy[""] = energy->getMachineTensor<MachineTensor>()->tensor[""];
  LOG(0, "RPA") << std::real(ctfEnergy.get_val());
}

