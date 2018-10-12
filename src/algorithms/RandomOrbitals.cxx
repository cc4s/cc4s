#include <algorithms/RandomOrbitals.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(RandomOrbitals);

RandomOrbitals::RandomOrbitals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

RandomOrbitals::~RandomOrbitals() {
}

void RandomOrbitals::run() {
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("FullParticleHoleCoulombVertex")
  );
  GammaGai->set_name("GammaGai");
  Tensor<complex> conjGammaGai(*GammaGai);
  conjugate(conjGammaGai);
  Tensor<complex> *cbbeta(
    getTensorArgument<complex>("RandomOrbitalsCoefficients")
  );
  int lens[] = { GammaGai->lens[0], cbbeta->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *GammaGalphai = new Tensor<complex>(
    3, lens, syms, *GammaGai->wrld, "GammaGalphai"
  );
  allocatedTensorArgument<complex>(
    "ParticleHoleCoulombVertex", GammaGalphai
  );
  LOG(1, "RandomOrbitals")
    << "Using Nr=" << cbbeta->lens[1] << " random orbitals to build " << GammaGalphai->get_name()
    << " from " << GammaGai->get_name() << " with Nv=" << GammaGai->lens[1] << " virtual orbitals" << std::endl;
   // NOTE: cbbeta are not conjugated for transforming psi*_a(x) in GammaGai
  (*GammaGalphai)["GAi"] = (*GammaGai)["Gai"] * (*cbbeta)["aA"];

  Tensor<complex> conjGammaGalphai(*GammaGalphai);
  conjugate(conjGammaGalphai);

  constexpr real eVPerHartree = 27.21138602;
  constexpr real angstromPerBohr = 0.52917721067;
  real volume(getRealArgument("Volume"));
  real gamma(getRealArgument("Gamma"));
  auto momenta(getTensorArgument<real>("Momenta"));
  // convert lengths of momenta into atomic units
  Transform<double>(
    std::function<void(double &)>(
      [volume](real &momentum) {
        momentum *= 2*Pi<>() * angstromPerBohr;
      }
    )
  )(
    (*momenta)["G"]
  );

  Tensor<complex> EG(1, lens, syms, *GammaGai->wrld, "EG");
  // compute high energy approximation with full virtual orbitals
  EG["G"] =
    conjGammaGai["Gai"] * (*GammaGai)["Gbj"] *
    conjGammaGai["Hai"] * (*GammaGai)["Hbj"];
  // multiply with approximate energy denominator = G^2/2 + G^2/2 + gamma
  // in ev
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [volume, gamma](real momentum, complex &e) {
        e = e / (gamma + momentum*momentum) * eVPerHartree;
      }
    )
  )(
    (*momenta)["G"], EG["G"]
  );
  // contract to approximate energy
  Scalar<> E;
  E[""] = EG["G"];
  real randomEnergy(E.get_val());
  LOG(0,"RandOrb") << "using virtual orbitals MP2d ~ " << randomEnergy << std::endl;

  // again with random orbitals
  EG["G"] =
    conjGammaGalphai["GAi"] * (*GammaGalphai)["GBj"] *
    conjGammaGalphai["HAi"] * (*GammaGalphai)["HBj"];
  // multiply with approximate energy denominator = G^2/2 + G^2/2 + gamma
  // in ev
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [volume, gamma](real momentum, complex &e) {
        e = e / (gamma + momentum*momentum) * eVPerHartree;
      }
    )
  )(
    (*momenta)["G"], EG["G"]
  );

  // contract to approximate energy
  E[""] = EG["G"];
  randomEnergy = E.get_val();
  LOG(0,"RandOrb") << "using random orbitals MP2d ~ " << randomEnergy << std::endl;
}

void RandomOrbitals::dryRun() {
  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("FullParticleHoleCoulombVertex")
  );
  DryTensor<complex> *cbbeta(
    getTensorArgument<complex, DryTensor<complex>>(
      "RandomOrbitalsCoefficients"
    )
  );
  int lens[] = { GammaGai->lens[0], cbbeta->lens[1], GammaGai->lens[2] }; 
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *GammaGalphai = new DryTensor<complex>(
    3, lens, syms, SOURCE_LOCATION
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleHoleCoulombVertex", GammaGalphai
  );
}

