#include <algorithms/RpaApxEnergy.hpp>

#include <Cc4s.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/SlicedCtfTensor.hpp>
#include <math/IterativePseudoInverse.hpp>
#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
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
      (*chiVFGn)["FGn"] <<= -spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      // particle/hole bubble propagating backwards, positive nu
      (*chiVFGn)["FGn"] += -spins *
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
  // TODO: distribute over nu
  for (int n(0); n < slicedChiVFGn.slicedLens[0]; ++n) {
    LOG(1,"RPA") << "evaluating imaginary frequency "
      << n << "/" << slicedChiVFGn.slicedLens[0] << std::endl;
    auto chiVFG( &slicedChiVFGn({n}) );
    auto PxVFG( &slicedPxVFGn({n}) );

    LapackMatrix<complex> laChiVFG(*chiVFG);
    // NOTE: chiV is not hermitian in the complex case
    LapackGeneralEigenSystem<complex> eigenSystem(laChiVFG);

    // write diagonalizaing transformation to U
    auto UFL(*chiVFG);
    eigenSystem.getRightEigenVectors().write(UFL);
    IterativePseudoInverse<complex> invUFL(UFL);
//    LOG(1,"RPA") << "right eigen error=" << eigenSystem.rightEigenError(laChiVFG) << std::endl;
//    LOG(1,"RPA") << "left eigen error=" << eigenSystem.leftEigenError(laChiVFG) << std::endl;
//    LOG(1,"RPA") << "biorthogonal error=" << eigenSystem.biorthogonalError() << std::endl;

    // write eigenvalues to CTF vector
    CTF::Vector<complex> lambdaL(eigenSystem.getEigenValues().size());
    int localNF(lambdaL.wrld->rank == 0 ? lambdaL.lens[0] : 0);
    std::vector<int64_t> lambdaIndices(localNF);
    for (int64_t i(0); i < localNF; ++i) { lambdaIndices[i] = i; }
    lambdaL.write(
      localNF, lambdaIndices.data(), eigenSystem.getEigenValues().data()
    );

    // compute matrix functions of chiV on their eigenvalues
    // Log(1-XV)+XV for RPA total energy
    CTF::Vector<complex> LogChiVL(lambdaL.lens[0]);
    CTF::Transform<complex, complex>(
      std::function<void(complex, complex &)>(
        [](complex chiV, complex &logChiV) {
          logChiV = std::log(1.0-chiV) + chiV;
        }
      )
    ) (
      lambdaL["L"], LogChiVL["L"]
    );

    // 1/(1-XV) for APX total energy
    CTF::Vector<complex> InvChiVL(lambdaL.lens[0]);
    CTF::Transform<complex, complex>(
      std::function<void(complex, complex &)>(
        [](complex chiV, complex &invChiV) {
          invChiV = 1.0 / (1.0-chiV);
        }
      )
    ) (
      lambdaL["L"], InvChiVL["L"]
    );

/*
    // check if UFL diagonalizes chiVFG
    UFL["FGn"] = UFL["FLn"] * lambdaL["L"] * invUFL.get()["LG"];
    UFL["FGn"] += (-1.0) * (*chiVFG)["FGn"];
    double norm = frobeniusNorm(UFL);
    LOG(1,"RPA") << "|X_F^G - U*_F^L lambda_L U_L^G|=" << norm << std::endl;
*/

    CTF::Scalar<complex> e;
    // Tr{Log(1-XV)+XV}
    e[""] = LogChiVL["L"];
    rpa += weights[n] * e.get_val();
    // Tr{Px(1-XV)^-1}
    e[""] = (*PxVFG)["FGn"] * UFL["GLn"] * InvChiVL["L"] * invUFL.get()["LF"];
    apx += weights[n] * e.get_val();
  }

  // 2 fold mirror symmetry, 1/Pi from +nu and -nu
  rpa *= 0.5 / Pi();
  apx *= 0.5 / Pi();
  LOG(1, "RPA") << "rpa=" << rpa << std::endl;
  LOG(1, "RPA") << "apx=" << apx << std::endl;
  setRealArgument("RpaApxEnergy", std::real(rpa+apx));
}

