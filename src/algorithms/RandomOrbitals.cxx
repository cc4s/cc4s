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
  Tensor<complex> *cbbeta(
    getTensorArgument<complex>("RandomOrbitalsCoefficients")
  );
  Tensor<complex> conjGammaGai(*GammaGai);
  conjugate(conjGammaGai);

  constexpr real eVPerHartree = 27.21138602;
  constexpr real angstromPerBohr = 0.52917721067;
//  constexpr real kineticEV =
//    eVPerHartree * pow(angstromPerBohr*2*Pi<>(),2);
  // volume in a.u.
  real volume(getRealArgument("Volume") / pow(angstromPerBohr,3));
  real gamma(getRealArgument("Gamma") * 2*Pi<>() * angstromPerBohr);
  auto momenta(getTensorArgument<real>("Momenta"));
  // convert lengths of momenta into atomic units
  Transform<double>(
    std::function<void(double &)>(
      [](real &momentum) {
        momentum *= 2*Pi<>() * angstromPerBohr;
      }
    )
  )(
    (*momenta)["G"]
  );

  // convert coulomb vertices into overlap densities
  Tensor<complex> SGai(*GammaGai);
  Tensor<complex> conjSGai(conjGammaGai);
  Transform<double, complex> overlapDensityFromCoulombVertex(
    std::function<void(double, complex &)>(
      [volume](real momentum, complex &gamma) {
        gamma /= sqrt(4*Pi<>() / pow(momentum,2) * eVPerHartree);
      }
    )
  );
  // (gives the unapproximated overlap densities)
  overlapDensityFromCoulombVertex((*momenta)["G"], SGai["Gai"]);
  overlapDensityFromCoulombVertex((*momenta)["G"], conjSGai["Gai"]);

  int QLens[] = {SGai.lens[0], SGai.lens[0]};
  int syms[] = { NS, NS, NS };
  Tensor<complex> QGH(2, QLens, syms, *GammaGai->wrld, "QGH");
  // compute high energy approximation with full virtual orbitals
  // compute right part of direct diagram, summing over b&j
  QGH["GH"] = SGai["Gbj"] * SGai["Hbj"];
  // copy right for computing left part
  Tensor<complex> conjQGH(QGH);
  // explicitly calculate left part from conjugate overlap densities
  // this could also be computed by: conjugate(conjQGH)
  conjQGH["GH"] = conjSGai["Gai"] * conjSGai["Hai"];
  // multiply left part together with right part to get total diagram
  QGH["GH"] *= conjQGH["GH"];
  real largeMomentumThreshold(getRealArgument("largeMomentumThreshold", 0));
  Scalar<complex> E;
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [largeMomentumThreshold](real momentum, complex &q) {
        if (largeMomentumThreshold < momentum) {
          q *= 4*Pi<>() / pow(momentum,2);
        } else {
          q = 0.0;
        }
      }
    )
  ) (
    // H is G'
    (*momenta)["H"], QGH["GH"]
  );
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [gamma, largeMomentumThreshold](real momentum, complex &q) {
        if (largeMomentumThreshold < momentum) {
          q *= 4*Pi<>() / pow(momentum,2) / (pow(gamma,2) + pow(momentum,2));
        } else {
          q = 0.0;
        }
      }
    )
  ) (
    (*momenta)["G"], QGH["GH"]
  );
  E[""] = 2.0 * eVPerHartree * QGH["GH"];
  complex randomEnergy(E.get_val());
  LOG(0,"RandOrb") << "using virtual orbitals: large momentum MP2d ~ "
    << randomEnergy << std::endl;


  ///////////////////////////////////////////////////////////////////////////
  // Random orbital stuff starts here...
  ///////////////////////////////////////////////////////////////////////////
  // again with random orbitals
  int lens[] = { GammaGai->lens[0], cbbeta->lens[1], GammaGai->lens[2] };
  Tensor<complex> SGalphai(
    3, lens, syms, *GammaGai->wrld, "SGalphai"
  );
  LOG(1, "RandomOrbitals")
    << "Using Nr=" << cbbeta->lens[1] << " random orbitals to build " << SGalphai.get_name()
    << " from " << GammaGai->get_name() << " with Nv=" << GammaGai->lens[1] << " virtual orbitals" << std::endl;
  // approximate overlap density \tilde S using random orbitals
  // NOTE: cbbeta are not conjugated for transforming psi*_a(x) in GammaGai
  SGalphai["GAi"] = SGai["Gai"] * (*cbbeta)["aA"];
  Tensor<complex> conjSGalphai(SGalphai);
  // also compute the conjugated overlap densities \tidle S
  conjugate(conjSGalphai);

  // compute high energy approximation with random orbitals
  // same equations as in the unapproximated case
  QGH["GH"] = SGalphai["Gbj"] * SGalphai["Hbj"];
  conjQGH["GH"] = conjSGalphai["Gai"] * conjSGalphai["Hai"];
  QGH["GH"] *= conjQGH["GH"];
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [largeMomentumThreshold](real momentum, complex &q) {
        if (largeMomentumThreshold < momentum) {
          q *= 4*Pi<>() / pow(momentum,2);
        } else {
          q = 0.0;
        }
      }
    )
  ) (
    (*momenta)["H"], QGH["GH"]
  );
  Transform<double, complex>(
    std::function<void(double, complex &)>(
      [gamma, largeMomentumThreshold](real momentum, complex &q) {
        if (largeMomentumThreshold < momentum) {
          q *= 4*Pi<>() / pow(momentum,2) / (pow(gamma,2) + pow(momentum,2));
        } else {
          q = 0.0;
        }
      }
    )
  ) (
    (*momenta)["G"], QGH["GH"]
  );
  E[""] = 2.0 * eVPerHartree * QGH["GH"];
  randomEnergy = E.get_val();
  LOG(0,"RandOrb") << "using random orbitals: large momentum MP2d ~ "
    << randomEnergy << std::endl;
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

