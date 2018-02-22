#include <algorithms/RpaApxEnergy.hpp>

#include <Cc4s.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/SlicedCtfTensor.hpp>
#include <util/LapackMatrix.hpp>
#include <util/LapackInverse.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <util/MpiCommunicator.hpp>
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
  chi0VFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chi0V");
  chi1VFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chi1V");
  double spins(2.0);
  tcc->compile(
    (
      // bubble with half V on both ends:
      // sign: 1xhole, 1xinteraction, 1xloop: (-1)^3
      // particle/hole bubble propagating forwards
      (*chi0VFGn)["FGn"] <<= -spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      // particle/hole bubble propagating backwards, positive nu
      (*chi0VFGn)["FGn"] += -spins *
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"],

      // adjacent pairs exchanged
      // sign: 2xholes, 2xinteraction, 1xloop: (-1)^5
      (*chi1VFGn)["FGn"] <<= -spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Haj"] * (*Pain)["ain"] *
        (*GammaFia)["Hib"] * (*conjGammaFia)["Gjb"] * (*conjPain)["bjn"],

      // compute Mp2 energy for benchmark of frequency grid
      // 2 fold rotational and 2 fold mirror symmetry, 1/Pi from +nu and -nu
      // sign: 1xdiagram: (-1)
      (*mp2Direct)[""] <<= -0.25 / Pi() *
        (*Wn)["n"] * (*chi0VFGn)["FGn"] * (*chi0VFGn)["GFn"],
      // 2 fold mirror symmetry only, 1/Pi from +nu and -nu
      (*mp2Exchange)[""] <<= -0.5 / Pi() * (*Wn)["n"] * (*chi1VFGn)["FFn"]
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

  LOG(1, "RPA") << "slicing along imaginary frequencies..." << std::endl;
  MpiCommunicator communicator(*wn->wrld);
  // slice CTF tensor of chi0V and chi1V along n (=3rd) dimension
  auto ctfChi0VFGn( &chi0VFGn->getMachineTensor<MT>()->tensor );
  auto ctfChi1VFGn( &chi1VFGn->getMachineTensor<MT>()->tensor );
  size_t sliceElementCount(ctfChi0VFGn->lens[0] * ctfChi0VFGn->lens[1]);
  std::map<int, std::vector<complex>> localChi0VFGn;
  std::map<int, std::vector<complex>> localChi1VFGn;
  for (
    int pass(0);
    pass < (ctfChi0VFGn->lens[2]-1) / communicator.getProcesses() + 1;
    ++pass
  ) {
    // get local slice number
    int n( pass * communicator.getProcesses() + communicator.getRank() );
    std::vector<int64_t> indices;
    if (n < ctfChi0VFGn->lens[2]) {
      indices.resize(sliceElementCount);
      localChi0VFGn[n].resize(sliceElementCount);
      localChi1VFGn[n].resize(sliceElementCount);
    }
    for (size_t i(0); i < indices.size(); ++i) {
      indices[i] = n*sliceElementCount + i;
    }
    ctfChi0VFGn->read(indices.size(), indices.data(), localChi0VFGn[n].data());
    ctfChi1VFGn->read(indices.size(), indices.data(), localChi1VFGn[n].data());
  }
  // we need non-hermitian diagonlization routines, which are not parallel,
  // so rather distribute over frequencies
  complex localRpa(0), localApx(0);
  for (
    int n(communicator.getRank());
    n < ctfChi0VFGn->lens[2];
    n += communicator.getProcesses()
  ) {
    LOG(1,"RPA") << "evaluating imaginary frequency "
      << n << "/" << ctfChi0VFGn->lens[2] << std::endl;
    LapackMatrix<complex> laChi0VFG(
      ctfChi0VFGn->lens[0], ctfChi0VFGn->lens[1], localChi0VFGn[n]
    );
    LapackMatrix<complex> laChi1VFG(
      ctfChi1VFGn->lens[0], ctfChi1VFGn->lens[1], localChi1VFGn[n]
    );

    // NOTE: chi0V is not hermitian in the complex case
    LapackGeneralEigenSystem<complex> chi0VEigenSystem(laChi0VFG);

    auto RFL( chi0VEigenSystem.getRightEigenVectors() );
    // compute inverse for eigen decomposition
    LapackInverse<complex> invRFL(RFL);

    // compute matrix functions of chi0V on their eigenvalues
    // Tr{-(Log(1-X0V)+X0V)} for RPA total energy
    auto chi0VL( &chi0VEigenSystem.getEigenValues() );
    complex trLogChi0V(0);
    for (size_t l(0); l < chi0VL->size(); ++l) {
      trLogChi0V -= std::log(1.0 - (*chi0VL)[l]) + (*chi0VL)[l];
    }

    // multiply RFL with 1/(1-XV) for each eigenvalue
    for (size_t l(0); l < chi0VL->size(); ++l) {
      for (int f(0); f < RFL.getRows(); ++f) {
        RFL(f,l) /= 1.0 - (*chi0VL)[l];
      }
    }
    // multiply with invRFL to get V^-1.W = (1-chi0V)^-1
    auto invVWFG( RFL * invRFL.get() );

    // setup chi1W for APX total energy
    auto chi1WFG( laChi1VFG * invVWFG );

    // diagonalize chi1W, now we only need the eigenvalues
    // NOTE: chi1W is not hermitian in the complex case
    LapackGeneralEigenSystem<complex> chi1WEigenSystem(chi1WFG, false);

    // compute matrix functions of chi1W on their eigenvalues
    // Tr{-Log(1-X1W)} for APX total energy
    auto chi1WL( &chi1WEigenSystem.getEigenValues() );
    complex trLogChi1W(0);
    for (size_t l(0); l < chi0VL->size(); ++l) {
      trLogChi1W -= std::log(1.0 - (*chi1WL)[l]);
    }

    localRpa += weights[n] * trLogChi0V;
    localApx += weights[n] * trLogChi1W;
  }

  // wait for all processes to finish their frequencies
  communicator.barrier();

  // reduce from all nodes
  complex rpa, apx;
  communicator.allReduce(localRpa, rpa);
  communicator.allReduce(localApx, apx);

  // 2 fold mirror symmetry, 1/Pi from +nu and -nu
  // sign: 1xdiagram: (-1)
  rpa *= -0.5 / Pi();
  apx *= -0.5 / Pi();
  LOG(1, "RPA") << "rpa=" << rpa << std::endl;
  LOG(1, "RPA") << "apx=" << apx << std::endl;
  setRealArgument("RpaApxEnergy", std::real(rpa+apx));
}

