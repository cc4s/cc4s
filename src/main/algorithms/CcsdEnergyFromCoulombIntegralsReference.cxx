#include <algorithms/CcsdEnergyFromCoulombIntegralsReference.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEnergyFromCoulombIntegralsReference);

Ptr<FockVector<Real<>, DryTensorEngine>>
CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DryTensorEngine>(iteration, amplitudes);
}

Ptr<FockVector<Real<>, DefaultTensorEngine>>
CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DefaultTensorEngine>(iteration, amplitudes);
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

template <typename F, typename TE>
Ptr<FockVector<F,TE>> CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration, const Ptr<const FockVector<F,TE>> &amplitudes
) {
  // get singles and doubles part of the amplitudes
  auto Tai( amplitudes->get(0) );
  auto Tabij( amplitudes->get(1) );

  // get amplitude parts
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  // construct residuum
  auto residuum( New<FockVector<F,TE>>(*amplitudes) );
  *residuum *= F(0);
  auto Rph( residuum->get(0) );
  auto Rpphh( residuum->get(1) );

  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));

  if (iteration == 0) {
   //TODO
   // && !isArgumentGiven("initialDoublesAmplitudes"))  {
    // For first iteration compute only the MP2 amplitudes
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  }
  return residuum;
}
