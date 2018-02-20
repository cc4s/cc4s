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
  chi0VFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chi0V");
  chi1VFGn = tcc->createTensor(std::vector<int>({NF,NF,Nn}), "chi1V");
  double spins(2.0);
  tcc->compile(
    (
      // bubble with half V on both ends:
      // particle/hole bubble propagating forwards
      (*chi0VFGn)["FGn"] <<= -spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Gai"] * (*Pain)["ain"],
      // particle/hole bubble propagating backwards, positive nu
      (*chi0VFGn)["FGn"] += -spins *
        (*GammaFia)["Fia"] * (*conjGammaFia)["Gia"] * (*conjPain)["ain"],

      // adjacent pairs exchanged
      (*chi1VFGn)["FGn"] <<= +spins *
        (*GammaFai)["Fai"] * (*conjGammaFai)["Haj"] * (*Pain)["ain"] *
        (*GammaFia)["Hib"] * (*conjGammaFia)["Gjb"] * (*conjPain)["bjn"],

      // compute Mp2 energy for benchmark of frequency grid
      // 2 fold rotational and 2 fold mirror symmetry, 1/Pi from +nu and -nu
      (*mp2Direct)[""] <<= -0.25 / Pi() *
        (*Wn)["n"] * (*chi0VFGn)["FGn"] * (*chi0VFGn)["GFn"],
      // 2 fold mirror symmetry only, 1/Pi from +nu and -nu
      (*mp2Exchange)[""] <<= +0.5 / Pi() * (*Wn)["n"] * (*chi1VFGn)["FFn"]
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

  // slice CTF tensor of chi0V and chi1V along n (=3rd) dimension
  SlicedCtfTensor<complex> slicedChi0VFGn(
    chi0VFGn->getMachineTensor<MT>()->tensor, {2}
  );
  SlicedCtfTensor<complex> slicedChi1VFGn(
    chi1VFGn->getMachineTensor<MT>()->tensor, {2}
  );
  complex rpa(0), apx(0);
  // TODO: distribute over nu
  for (int n(0); n < slicedChi0VFGn.slicedLens[0]; ++n) {
    LOG(1,"RPA") << "evaluating imaginary frequency "
      << n << "/" << slicedChi0VFGn.slicedLens[0] << std::endl;
    auto chi0VFG( &slicedChi0VFGn({n}) );
    auto chi1VFG( &slicedChi1VFGn({n}) );

    LapackMatrix<complex> laChi0VFG(*chi0VFG);
    // NOTE: chi0V is not hermitian in the complex case
    LapackGeneralEigenSystem<complex> chi0VEigenSystem(laChi0VFG);

    // write diagonalizaing transformation to U
    auto UFL(*chi0VFG);
    chi0VEigenSystem.getRightEigenVectors().write(UFL);
    // compute pseudo inverse for eigendecomposition
    IterativePseudoInverse<complex> invUFL(UFL);
//    LOG(1,"RPA") << "right eigen error=" << chi0VEigenSystem.rightEigenError(laChi0VFG) << std::endl;
//    LOG(1,"RPA") << "left eigen error=" << chi0VEigenSystem.leftEigenError(laChi0VFG) << std::endl;
//    LOG(1,"RPA") << "biorthogonal error=" << chi0VEigenSystem.biorthogonalError() << std::endl;

    // write eigenvalues to CTF vector
    CTF::Vector<complex> chi0VL(chi0VEigenSystem.getEigenValues().size());
    int localNF(chi0VL.wrld->rank == 0 ? chi0VL.lens[0] : 0);
    std::vector<int64_t> lambdaIndices(localNF);
    for (int64_t i(0); i < localNF; ++i) { lambdaIndices[i] = i; }
    chi0VL.write(
      localNF, lambdaIndices.data(), chi0VEigenSystem.getEigenValues().data()
    );

    // compute matrix functions of chi0V on their eigenvalues
    // Log(1-X0V)+X0V for RPA total energy
    CTF::Vector<complex> LogChi0VL(chi0VL.lens[0]);
    CTF::Transform<complex, complex>(
      std::function<void(complex, complex &)>(
        [](complex chi0V, complex &logChi0V) {
          logChi0V = std::log(1.0-chi0V) + chi0V;
        }
      )
    ) (
      chi0VL["L"], LogChi0VL["L"]
    );

    // 1/(1-XV) for W
    CTF::Vector<complex> InvChi0VL(chi0VL.lens[0]);
    CTF::Transform<complex, complex>(
      std::function<void(complex, complex &)>(
        [](complex chi0V, complex &invChi0V) {
          invChi0V = 1.0 / (1.0-chi0V);
        }
      )
    ) (
      chi0VL["L"], InvChi0VL["L"]
    );

    // setup chi1W for APX total energy
    auto chi1WFG(*chi1VFG);
    chi1WFG["FGn"] =
      (*chi1VFG)["FHn"] * UFL["HLn"] * InvChi0VL["L"] * invUFL.get()["LG"];

    // diagonalize chi1W
    LapackMatrix<complex> laChi1WFG(chi1WFG);
    // NOTE: chi1W is not hermitian in the complex case
    LapackGeneralEigenSystem<complex> chi1WEigenSystem(laChi1WFG);

    // write eigenvalues to CTF vector
    CTF::Vector<complex> chi1WL(chi1WEigenSystem.getEigenValues().size());
    chi1WL.write(
      localNF, lambdaIndices.data(), chi1WEigenSystem.getEigenValues().data()
    );

    // compute matrix functions of chi1W on their eigenvalues
    // -Log(1-X1W)for APX total energy
    CTF::Vector<complex> LogChi1WL(chi1WL.lens[0]);
    CTF::Transform<complex, complex>(
      std::function<void(complex, complex &)>(
        [](complex chi1W, complex &logChi1W) {
          logChi1W = std::log(1.0-chi1W);
//          logChi1W = std::log(1.0-chi1W);
//          LOG(1,"RPA") << "chi1W(L)=" << chi1W << std::endl;
        }
      )
    ) (
      chi1WL["L"], LogChi1WL["L"]
    );


    CTF::Scalar<complex> e;
    // Tr{Log(1-X0V)+XV}
    e[""] = LogChi0VL["L"];
    rpa += weights[n] * e.get_val();

    // Tr{Log(1-X1W)}
    e[""] = LogChi1WL["L"];
    apx += weights[n] * e.get_val();
  }

  // 2 fold mirror symmetry, 1/Pi from +nu and -nu
  rpa *= +0.5 / Pi();
  apx *= -0.5 / Pi();
  LOG(1, "RPA") << "rpa=" << rpa << std::endl;
  LOG(1, "RPA") << "apx=" << apx << std::endl;
  setRealArgument("RpaApxEnergy", std::real(rpa+apx));
}

