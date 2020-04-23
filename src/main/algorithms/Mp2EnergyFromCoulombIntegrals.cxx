#include <algorithms/Mp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EnergyFromCoulombIntegrals);

Ptr<MapNode> Mp2EnergyFromCoulombIntegrals::run(const Ptr<MapNode> &arguments) {
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto orbitals(coulombIntegrals->getValue<std::string>("orbitals"));
  // multiplex calls to template methods
  if (orbitals == "real") {
    if (Cc4s::options->dryRun) {
      return calculateMp2Energy<Real<>,DryTensorEngine>(arguments);
    } else {
      return calculateMp2Energy<Real<>,DefaultTensorEngine>(arguments);
    }
  } else if (orbitals == "complex") {
    if (Cc4s::options->dryRun) {
      return calculateMp2Energy<Complex<>,DryTensorEngine>(arguments);
    } else {
      return calculateMp2Energy<Complex<>,DefaultTensorEngine>(arguments);
    }
  } else {
    Assert(false, "unsupported orbitals type '" + orbitals + "'");
  }
}

template <typename F, typename TE>
Ptr<MapNode> Mp2EnergyFromCoulombIntegrals::calculateMp2Energy(
  const Ptr<MapNode> &arguments
) {
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto slices(coulombIntegrals->getMap("slices"));
  auto Vabij(slices->getValue<Ptr<TensorRecipe<F,TE>>>("abij"));
  Real<> spins(coulombIntegrals->getValue<size_t>("spins"));

  auto holeEnergies(arguments->getMap("holeEigenEnergies"));
  auto epsi(holeEnergies->getValue<Ptr<Tensor<Real<>,TE>>>("data"));
  auto particleEnergies(arguments->getMap("particleEigenEnergies"));
  auto epsa(particleEnergies->getValue<Ptr<Tensor<Real<>,TE>>>("data"));

  auto No(epsi->lens[0]);
  auto Nv(epsa->lens[0]);
  auto Dabij(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Dabij")
  );
  auto direct( Tcc<TE>::template tensor<F>("D") );
  auto exchange( Tcc<TE>::template tensor<F>("X") );
  (
    (*Dabij)["abij"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"]),
    (*Dabij)["abij"] +=  map<F>([](Real<> eps) {return F(eps);}, (*epsa)["b"]),
    (*Dabij)["abij"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"]),
    (*Dabij)["abij"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsi)["j"]),
    (*Dabij)["abij"] <<=
      map<F>(conj<F>, (*Vabij)["abij"]) *
      map<F>([](F delta) { return F(1)/delta; }, (*Dabij)["abij"]),
    (*direct)[""] <<= -0.5*spins*spins * (*Vabij)["abij"] * (*Dabij)["abij"],
    (*exchange)[""] <<= +0.5*spins * (*Vabij)["abji"] * (*Dabij)["abij"]
  )->compile()->execute();

  F D, X;
  D = direct->read();
  X = exchange->read();
  LOG(1,getName()) << "direct=" << D << std::endl;
  LOG(1,getName()) << "exchange=" << X << std::endl;
  
  auto energy(New<MapNode>());
  energy->setValue<Real<>>("direct", real<F>(D));
  energy->setValue<Real<>>("exchange", real<F>(X));
  energy->setValue<Real<>>("value", real<F>(D+X));
  energy->setValue<Real<>>("unit", holeEnergies->getValue<Real<>>("unit"));
  auto result(New<MapNode>());
  result->get("energy") = energy;
  return result;
}

/*
template <typename F>
F Mp2EnergyFromCoulombIntegrals::calculateMp2Energy(CTF::Tensor<F> &Vabij) {
  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));

  // convert to type F (either complex or double)
  Tensor<F> Fepsi(1, &epsi->lens[0], epsi->sym, *epsi->wrld, "Fepsi");
  // NOTE: just copies if both arguments are real
  toComplexTensor(*epsi, Fepsi);
  Tensor<F> Fepsa(1, &epsa->lens[0], epsa->sym, *epsa->wrld, "Fepsa");
  toComplexTensor(*epsa, Fepsa);

  // create excitation energy
  auto Tabij(new Tensor<F>(Vabij));
  //  Tensor<F> Tabij(false, Vabij);
  Tabij->set_name("Tabij");
  (*Tabij)["abij"] =  Fepsi["i"];
  (*Tabij)["abij"] += Fepsi["j"];
  (*Tabij)["abij"] -= Fepsa["a"];
  (*Tabij)["abij"] -= Fepsa["b"];

  // use transform to divide Vabij by Tabij and store in Tabij
  CTF::Transform<F, F>(
    std::function<void(F, F &)>(
      [](F vabij, F &tabij) {
        tabij = conj(vabij / tabij);
      }
    )
  ) (
    Vabij["abij"], (*Tabij)["abij"]
  );
*/
// TODO: below method requires less memory but does not work with current CTF
/*
  Matrix<> Dai(Vabij.lens[0], Vabij.lens[2], NS, *Vabij.wrld, "Dai");
  Dai["ai"] =  (*epsi)["i"];
  Dai["ai"] -= (*epsa)["a"];

  Tensor<F> Tabij(Vabij);
  // use transform to divide real/complex Tabij by real Dabij = Dai+Dbj
  CTF::Transform<double, double, F>(
    std::function<void(double, double, F &)>(
      [](double dai, double dbj, F &t) {
        t = conj(t / (dai + dbj));
      }
    )
  ) (
    Dai["ai"], Dai["bj"], Tabij["abij"]
  );
*/
/*
  Scalar<F> energy(*Cc4s::world);
  F dire, exce;

  energy[""] = 0.5 * spins * spins * (*Tabij)["abij"] * Vabij["abij"];
  dire = energy.get_val();
  energy[""] = 0.5 * spins * (*Tabij)["abji"] * Vabij["abij"];
  exce = -1.0 * energy.get_val();
  LOG(0, "MP2") << "e=" << (dire + exce) << std::endl;
  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  // use transform to divide Vabij by Tabij and store in Tabij
  CTF::Transform<F, F>(
    std::function<void(F, F &)>(
      [](F vabij, F &tabij) {
        tabij = conj(tabij);
      }
    )
  ) (
    Vabij["abij"], (*Tabij)["abij"]
  );


  if (isArgumentGiven("Mp2DoublesAmplitudes")) {
    allocatedTensorArgument<F>("Mp2DoublesAmplitudes", Tabij);
  } else {
    delete Tabij;
  }

  return dire + exce;
}
*/

