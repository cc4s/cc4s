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
#include <algorithm>

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
  auto Wn(tcc->createTensor(MachineTensor::create(complexWn)));

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
      [](double nu, complex &d){ d = 1.0 / (d - complex(0,1)*nu); }
    )
  ) (
    (*realNun)["n"], ctfPain["ain"]
  );
  CTF::Tensor<complex> ctfConjPain(ctfPain);
  conjugate(ctfConjPain);

  CTF::Tensor<complex> ctfEai(2, std::vector<int>({Nv,No}).data());
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double eps, complex &d){ d = eps; }
    )
  ) (
    (*epsa)["a"], ctfEai["ai"]
  );
  CTF::Transform<double, complex>(
    std::function<void(double, complex &)>(
      [](double eps, complex &d){ d -= eps; }
    )
  ) (
    (*epsi)["i"], ctfEai["ai"]
  );
  CTF::Transform<complex>(
    std::function<void(complex &)>(
      [](complex &d){ d = -1.0 / d; }
    )
  ) (
    ctfEai["ai"]
  );
  auto Eai(tcc->createTensor(MachineTensor::create(ctfEai)));

  auto Pain(tcc->createTensor(MachineTensor::create(ctfPain)));
  auto conjPain(tcc->createTensor(MachineTensor::create(ctfConjPain)));
  auto Eain(tcc->createTensor(MachineTensor::create(ctfPain)));

  auto energy(tcc->createTensor(std::vector<int>(), "energy"));
  auto chiVFGn(tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chiV"));
  tcc->compile(
    (
      (*chiVFGn)["FGn"] <<=
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      (*chiVFGn)["FGn"] +=
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"],
      (*energy)[""] <<= (*Wn)["n"] * (*chiVFGn)["FGn"] * (*chiVFGn)["GFn"]
/*
      (*Eain)["ain"] <<=(*Pain)["ain"],
      (*Eain)["ain"] += (*conjPain)["ain"],
      (*Eain)["ain"] <<= (*Eain)["ain"] * (*Eain)["ain"],
      (*Eai)["ai"] += 1 / Pi() * (*Wn)["n"] * (*Eain)["ain"],
      (*energy)[""] <<= 0.5 * (*Eai)["ai"] * (*Eai)["ai"]
*/
    )
  )->execute();

  CTF::Scalar<complex> ctfEnergy;
  ctfEnergy[""] = energy->getMachineTensor<MachineTensor>()->tensor[""];
  complex e(-ctfEnergy.get_val() / Pi());
//  complex e(ctfEnergy.get_val());
  LOG(0, "RPA") << e << std::endl;
  setRealArgument("Mp2Energy", std::real(e));
}

