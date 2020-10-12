#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombIntegrals)

Ptr<FockVector<Real<>, DryTensorEngine>>
DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DryTensorEngine>(iteration, amplitudes);
}
Ptr<FockVector<Complex<>, DryTensorEngine>>
DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DryTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DryTensorEngine>(iteration, amplitudes);
}

Ptr<FockVector<Real<>, DefaultTensorEngine>>
DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DefaultTensorEngine>(iteration, amplitudes);
}
Ptr<FockVector<Complex<>, DefaultTensorEngine>>
DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DefaultTensorEngine>(iteration, amplitudes);
}


template <typename F, typename TE>
Ptr<FockVector<F,TE>> DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration, const Ptr<const FockVector<F,TE>> &amplitudes
) {
  // read all required integrals
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));
  auto Vphhp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phhp"));
  auto Vhhpp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhpp"));

  // get spins
  auto orbitalType(
    coulombIntegrals->getMap(
      "indices"
    )->getMap("orbital")->getValue<std::string>("type")
  );
  Real<> spins;
  if (orbitalType == "spatial") {
    spins = 2;
  } else if (orbitalType == "spin") {
    spins = 1;
  } else {
    ASSERT_LOCATION(
      false, "unsupported orbital type '" + orbitalType + "'",
      coulombIntegrals->getMap(
        "indices"
      )->getMap("orbital")->get("type")->sourceLocation
    );
  }

  // get amplitude parts
  auto Tpphh( amplitudes->get(1) );

  // construct residuum
  auto residuum( New<FockVector<F,TE>>(*amplitudes) );
  *residuum *= F(0);
  auto Rpphh( residuum->get(1) );

  bool linearized(arguments->getValue<bool>("linearized", false));
  bool adjacentPairsExchange(
    arguments->getValue<bool>("adjacentPairsExchange", false)
  );
  if (linearized) {
//    OUT() << "Solving linearized T2 Amplitude Equations" << std::endl;
  } else {
//    OUT() << "Solving T2 Amplitude Equations" << std::endl;
  }

  // TODO: deal with starting amplitudes
  if (iteration > 0 || restart) { // || isArgumentGiven("startingDoublesAmplitudes")) {
    auto Whhpp( Tcc<TE>::template tensor<F>("Whhpp") );
    // for the remaining iterations compute the drCCD residuum
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"],
      (*Rpphh)["abij"] += spins * (*Vphhp)["akic"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += spins * (*Vphhp)["bkjc"] * (*Tpphh)["acik"],
      (linearized) ? (
        // linearized: nothing more to do
        Tcc<TE>::sequence()
      ) : (
        // otherwise: do quadratic contribution
        (*Whhpp)["ijab"] <<= spins * (*Vhhpp)["ijab"],
        (adjacentPairsExchange) ? (
          // adjacent pairs correction: also exchange holes in Whhpp
          (*Whhpp)["ijab"] -= (*Vhhpp)["jiab"],
          Tcc<TE>::sequence()
        ) : (
          // otherwise: do nothing else with Whhpp
          Tcc<TE>::sequence()
        ),
        // compute quadratic contribution
        (*Rpphh)["abij"] +=
          spins * (*Whhpp)["klcd"] * (*Tpphh)["acik"] * (*Tpphh)["dblj"],
        Tcc<TE>::sequence()
      )
    )->execute();
// TODO: adjacent pairs exchange
/*
      Tensor<F> Wijab(false, *Vijab);
      Wijab["ijab"] = ;
      if (getIntegerArgument("adjacentPairsExchange", 0)) {
        Wijab["ijab"] -= (*Vijab)["jiab"];
      }
*/
  } else {
    // no amplitudes given: start with MP2 amplitudes
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  }

  return residuum;
}

