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

Ptr<FockVector<Complex<>, DryTensorEngine>>
CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DryTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DryTensorEngine>(iteration, amplitudes);
}

Ptr<FockVector<Real<>, DefaultTensorEngine>>
CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DefaultTensorEngine>(iteration, amplitudes);
}


Ptr<FockVector<Complex<>, DefaultTensorEngine>>
CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DefaultTensorEngine>(iteration, amplitudes);
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

  auto onlyPpl(arguments->getValue<size_t>("onlyPpl", 0) );


  if ((iteration == 0) && (onlyPpl == 0)) {

    //TODO
    // && !isArgumentGiven("initialDoublesAmplitudes"))  {
    // For first iteration compute only the MP2 amplitudes
    // FIXME: if does nothing currently
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  } else if (onlyPpl == 1) {

    LOG(1, getCapitalizedAbbreviation())
        << "Calculate only PPL diagrams" << std::endl;

    auto Vpppp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pppp"));
    auto Vppph(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("ppph"));
    auto Xabcd( Tcc<TE>::template tensor<F>("Xabcd") );

    COMPILE(
      // Build Xabcd intermediate
      (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vppph)["cdak"] * (*Tai)["bk"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vppph)["dcbk"] * (*Tai)["ak"],

      (*Rpphh)["abij"] <<= (*Xabcd)["abcd"] * (*Tabij)["cdij"],
      (*Rpphh)["abij"]  += (*Xabcd)["abcd"] * (*Tai)["ci"] * (*Tai)["dj"]
    )->execute();

  } else {
    auto Vpppp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pppp"));
    auto Vphph(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phph"));
    auto Vhhhh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhhh"));
    auto Vhhhp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhhp"));
    auto Vppph(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("ppph"));
    auto Vhhpp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhpp"));
    auto Vpphp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphp"));
    auto Vphhh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phhh"));
    auto Vhphp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hphp"));
    auto Vphhp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phhp"));
    auto Vphpp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phpp"));
    auto Vhhph(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhph"));
    auto Vhppp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hppp"));
    // Hirata intermediates
    auto Lac( Tcc<TE>::template tensor<F>("Lac") );
    auto Kac( Tcc<TE>::template tensor<F>("Kac") );
    auto Lki( Tcc<TE>::template tensor<F>("Lki") );
    auto Kki( Tcc<TE>::template tensor<F>("Kki") );
    auto Kck( Tcc<TE>::template tensor<F>("Kck") );
    auto Xklij( Tcc<TE>::template tensor<F>("Xklij") );
    auto Xakci( Tcc<TE>::template tensor<F>("Xakci") );
    auto Xakic( Tcc<TE>::template tensor<F>("Xakic") );
    auto Xabcd( Tcc<TE>::template tensor<F>("Xabcd") );

    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;

    COMPILE(
      // Build Kac
      (*Kac)["ac"] <<= (-2.0) * (*Vhhpp)["klcd"] * (*Tabij)["adkl"],
      (*Kac)["ac"] += ( 1.0) * (*Vhhpp)["kldc"] * (*Tabij)["adkl"],
      (*Kac)["ac"] += (-2.0) * (*Vhhpp)["klcd"] * (*Tai)["ak"] * (*Tai)["dl"],
      (*Kac)["ac"] += ( 1.0) * (*Vhhpp)["kldc"] * (*Tai)["ak"] * (*Tai)["dl"],

      // Build Lac
      (*Lac)["ac"] <<= (*Kac)["ac"],
      (*Lac)["ac"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tai)["dk"],
      (*Lac)["ac"] += (-1.0) * (*Vphpp)["akdc"] * (*Tai)["dk"],

      // Build Kki
      (*Kki)["ki"] <<= (2.0) * (*Vhhpp)["klcd"] * (*Tabij)["cdil"],
      (*Kki)["ki"] += (-1.0) * (*Vhhpp)["kldc"] * (*Tabij)["cdil"],
      (*Kki)["ki"] += ( 2.0) * (*Vhhpp)["klcd"] * (*Tai)["ci"] * (*Tai)["dl"],
      (*Kki)["ki"] += (-1.0) * (*Vhhpp)["kldc"] * (*Tai)["ci"] * (*Tai)["dl"],

      // Build Lki
      (*Lki)["ki"] <<= (*Kki)["ki"],
      (*Lki)["ki"] += ( 2.0) * (*Vhhhp)["klic"] * (*Tai)["cl"],
      (*Lki)["ki"] += (-1.0) * (*Vhhph)["klci"] * (*Tai)["cl"],

      // Contract Lac with T2 Amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Lac)["ac"] * (*Tabij)["cbij"],

      // Contract Lki with T2 Amplitudes
      (*Rpphh)["abij"] += (-1.0) * (*Lki)["ki"] * (*Tabij)["abkj"],

      // Contract Coulomb integrals with T2 amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Vpphp)["abic"] * (*Tai)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vhphp)["kbic"] * (*Tai)["ak"] * (*Tai)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhh)["akij"] * (*Tai)["bk"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhp)["akic"] * (*Tai)["cj"] * (*Tai)["bk"],

      // Build Xakic
      (*Xakic)["akic"] <<= (*Vphhp)["akic"],
      (*Xakic)["akic"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tai)["al"],
      (*Xakic)["akic"] += ( 1.0) * (*Vphpp)["akdc"] * (*Tai)["di"],
      (*Xakic)["akic"] += (-0.5) * (*Vhhpp)["lkdc"] * (*Tabij)["dail"],
      (*Xakic)["akic"] += (-1.0) * (*Vhhpp)["lkdc"] * (*Tai)["di"] * (*Tai)["al"],
      (*Xakic)["akic"] += ( 1.0) * (*Vhhpp)["lkdc"] * (*Tabij)["adil"],
      (*Xakic)["akic"] += (-0.5) * (*Vhhpp)["lkcd"] * (*Tabij)["adil"],
      (*Rpphh)["abij"] += ( 2.0) * (*Xakic)["akic"] * (*Tabij)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakic)["akic"] * (*Tabij)["bckj"],

      // Build Xakci
      (*Xakci)["akci"] <<= (*Vphph)["akci"],
      (*Xakci)["akci"] += (-1.0) * (*Vhhph)["lkci"] * (*Tai)["al"],
      (*Xakci)["akci"] += ( 1.0) * (*Vphpp)["akcd"] * (*Tai)["di"],
      (*Xakci)["akci"] += (-0.5) * (*Vhhpp)["lkcd"] * (*Tabij)["dail"],
      (*Xakci)["akci"] += (-1.0) * (*Vhhpp)["lkcd"] * (*Tai)["di"] * (*Tai)["al"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tabij)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tabij)["ackj"],

      // Symmetrize Rpphh by applying permutation operator
      (*Rpphh)["abij"] += (*Rpphh)["baji"],

      //////////////////////////////////////////////////////////////////////
      // Now add all terms to Rpphh that do not need to be symmetrized with
      // the permutation operator
      //////////////////////////////////////////////////////////////////////

      // Rpphh are the Tabij amplitudes for the next iteration and need to be build
      (*Rpphh)["abij"] += (*Vpphh)["abij"],

      // Build Xklij intermediate
      (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
      (*Xklij)["klij"] += (*Vhhhp)["klic"] * (*Tai)["cj"],
      (*Xklij)["klij"] += (*Vhhph)["klcj"] * (*Tai)["ci"],
      (*Xklij)["klij"] += (*Vhhpp)["klcd"] * (*Tabij)["cdij"],
      (*Xklij)["klij"] += (*Vhhpp)["klcd"] * (*Tai)["ci"] * (*Tai)["dj"],

      // Contract Xklij with T2 Amplitudes
      (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tabij)["abkl"],

      // Contract Xklij with T1 Amplitudes
      (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tai)["ak"] * (*Tai)["bl"],

      // Build Xabcd intermediate
      (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vphpp)["akcd"] * (*Tai)["bk"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vhppp)["kbcd"] * (*Tai)["ak"],

      // Contract Xabcd with T2 and T1 Amplitudes
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tabij)["cdij"],
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tai)["ci"] * (*Tai)["dj"]
    )->execute();

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T1 Amplitude Equations" << std::endl;
    COMPILE(
      // Contract Kac and Kki with T1 amplitudes
      (*Rph)["ai"] += ( 1.0) * (*Kac)["ac"] * (*Tai)["ci"],
      (*Rph)["ai"] += (-1.0) * (*Kki)["ki"] * (*Tai)["ak"],

      //Build Kck
      (*Kck)["ck"] <<= ( 2.0) * (*Vhhpp)["klcd"] * (*Tai)["dl"],
      (*Kck)["ck"]  += (-1.0) * (*Vhhpp)["kldc"] * (*Tai)["dl"],

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rph)["ai"] += ( 2.0) * (*Kck)["ck"] * (*Tabij)["caki"],
      (*Rph)["ai"] += (-1.0) * (*Kck)["ck"] * (*Tabij)["caik"],
      (*Rph)["ai"] += ( 1.0) * (*Kck)["ck"] * (*Tai)["ci"] * (*Tai)["ak"],
      (*Rph)["ai"] += ( 2.0) * (*Vphhp)["akic"] * (*Tai)["ck"],
      (*Rph)["ai"] += (-1.0) * (*Vphph)["akci"] * (*Tai)["ck"],
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tabij)["cdik"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tabij)["cdik"],
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tai)["ci"] * (*Tai)["dk"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tai)["ci"] * (*Tai)["dk"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tabij)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tabij)["ackl"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tai)["ak"] * (*Tai)["cl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tai)["ak"] * (*Tai)["cl"]
    )->execute();

  }
  return residuum;
}
