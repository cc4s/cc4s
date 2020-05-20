#include <algorithms/CcsdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEnergyFromCoulombIntegrals);

Ptr<FockVector<Real<>, DryTensorEngine>>
CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
) {
  return getResiduum<DryTensorEngine>(iteration, amplitudes);
}
// TODO: overrides for Complex<> types. Currently the Reference methods are used
/*
Ptr<FockVector<Complex<>, DryTensorEngine>>
CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DryTensorEngine>> &amplitudes
) {
  return CcsdEnergyFromCoulombIntegralsReference().getResiduum<
    Complex<>, DryTensorEngine
  >(
    iteration, amplitudes
  );
}
*/

Ptr<FockVector<Real<>, DefaultTensorEngine>>
CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<DefaultTensorEngine>(iteration, amplitudes);
}

/*
Ptr<FockVector<Complex<>, DefaultTensorEngine>>
CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration,
  const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
) {
  return CcsdEnergyFromCoulombIntegralsReference().getResiduum<
    Complex<>, DefaultTensorEngine
  >(
    iteration, amplitudes
  );
}
*/

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

template <typename TE>
Ptr<FockVector<Real<>,TE>> CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration, const Ptr<const FockVector<Real<>,TE>> &amplitudes
) {

  // get amplitude parts
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  // construct residuum
  auto residuum( New<FockVector<Real<>,TE>>(*amplitudes) );
  *residuum *= 0.0;
  auto Rph( residuum->get(0) );
  auto Rpphh( residuum->get(1) );

  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("pphh"));

  if (iteration == 0) {
    //TODO
    // && !isArgumentGiven("initialDoublesAmplitudes"))  {
    // For first iteration compute only the MP2 amplitudes
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  } else {
    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;
    auto slicedCoulombVertex(arguments->getMap("slicedCoulombVertex"));
    auto slices(slicedCoulombVertex->getMap("slices"));
    auto orbitals(slicedCoulombVertex->getValue<std::string>("orbitals"));
    auto GammaGpp(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("pp"));
    auto GammaGph(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("ph"));
    auto GammaGhh(slices->getValue<Ptr<TensorRecipe<Complex<>,TE>>>("hh"));

    auto Vphph(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("phph"));
    auto Vhhhh(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("hhhh"));
    auto Vhhhp(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("hhhp"));
    //intermediate integrals
    auto Vphpp(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("phpp"));
    auto Vhppp(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("hppp"));
    auto Vpppp(coulombSlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("pppp"));

    //Gamma -> Real/Imag
    auto realGammaGpp( Tcc<TE>::template tensor<Real<>>("realGammaGpp") );
    auto imagGammaGpp( Tcc<TE>::template tensor<Real<>>("imagGammaGpp") );
    auto realGammaGph( Tcc<TE>::template tensor<Real<>>("realGammaGph") );
    auto imagGammaGph( Tcc<TE>::template tensor<Real<>>("imagGammaGph") );
    auto realGammaGhh( Tcc<TE>::template tensor<Real<>>("realGammaGhh") );
    auto imagGammaGhh( Tcc<TE>::template tensor<Real<>>("imagGammaGhh") );

    COMPILE(
      (*realGammaGpp)["Gab"] <<= map<Real<>>(real<Complex<>>, (*GammaGpp)["Gab"] ),
      (*imagGammaGpp)["Gab"] <<= map<Real<>>(imag<Complex<>>, (*GammaGpp)["Gab"] ),
      (*realGammaGph)["Gai"] <<= map<Real<>>(real<Complex<>>, (*GammaGph)["Gai"] ),
      (*imagGammaGph)["Gai"] <<= map<Real<>>(imag<Complex<>>, (*GammaGph)["Gai"] ),
      (*realGammaGhh)["Gij"] <<= map<Real<>>(real<Complex<>>, (*GammaGhh)["Gij"] ),
      (*imagGammaGhh)["Gij"] <<= map<Real<>>(imag<Complex<>>, (*GammaGhh)["Gij"] )
    )->execute();
    auto realDressedGammaGpp( Tcc<TE>::template tensor<Real<>>("realDressedGammaGpp") );
    auto imagDressedGammaGpp( Tcc<TE>::template tensor<Real<>>("imagDressedGammaGpp") );
    auto realDressedGammaGph( Tcc<TE>::template tensor<Real<>>("realDressedGammaGph") );
    auto imagDressedGammaGph( Tcc<TE>::template tensor<Real<>>("imagDressedGammaGph") );
    auto realDressedGammaGhh( Tcc<TE>::template tensor<Real<>>("realDressedGammaGhh") );
    auto imagDressedGammaGhh( Tcc<TE>::template tensor<Real<>>("imagDressedGammaGhh") );
    auto realDressedGammaGphTwo(Tcc<TE>::template tensor<Real<>>("realDressedGammaGhhTwo") );
    auto imagDressedGammaGphTwo(Tcc<TE>::template tensor<Real<>>("imagDressedGammaGhhTwo") );
    // define intermediates
    auto Kac( Tcc<TE>::template tensor<Real<>>("Kac") ); //kappa_ac
    auto Kki( Tcc<TE>::template tensor<Real<>>("Kki") ); //kappa_ki
    auto Lac( Tcc<TE>::template tensor<Real<>>("Lac") ); //lambda_ac
    auto Lki( Tcc<TE>::template tensor<Real<>>("Lki") ); //lambda_ki
    auto Xabij( Tcc<TE>::template tensor<Real<>>("Xabij") ); // T2+T1*T1
    auto Yabij( Tcc<TE>::template tensor<Real<>>("Yabij") ); // T2+2*T1*T1
    auto Zabij( Tcc<TE>::template tensor<Real<>>("Zabij") );
    auto Xakic( Tcc<TE>::template tensor<Real<>>("Xakic") );
    auto Xakci( Tcc<TE>::template tensor<Real<>>("Xakci") );
    auto Xklij( Tcc<TE>::template tensor<Real<>>("Xklij") );
    auto Kck( Tcc<TE>::template tensor<Real<>>("Kck") );  // T1 intermediate
    //intermediate intermediate
    auto Xabcd( Tcc<TE>::template tensor<Real<>>("Xabcd") );
    COMPILE(
      (*Xabij)["abij"] <<= (*Tpphh)["abij"],
      (*Xabij)["abij"] += (*Tph)["ai"] * (*Tph)["bj"],
      (*Yabij)["abij"] <<= (*Tpphh)["abij"],
      (*Yabij)["abij"] += ( 2.0) * (*Tph)["ai"] * (*Tph)["bj"],
      /////////////////////////////////////
      // Lac and Kac for doubles amplitudes
      /////////////////////////////////////
      // Build Kac
      (*Kac)["ac"] <<= (-2.0) * (*Vpphh)["cdkl"] * (*Xabij)["adkl"],
      (*Kac)["ac"]  += ( 1.0) * (*Vpphh)["dckl"] * (*Xabij)["adkl"],
      // Lac
      (*Lac)["ac"] <<= (*Kac)["ac"],
      // TODO distiguish
      (*Lac)["ac"] +=
        ( 2.0) * (*realGammaGpp)["Gca"] * (*realGammaGph)["Gdk"] * (*Tph)["dk"],
      (*Lac)["ac"] +=
        ( 2.0) * (*imagGammaGpp)["Gca"] * (*imagGammaGph)["Gdk"] * (*Tph)["dk"],
      //TODO use an already introduced intermediate
      (*Lac)["ac"] +=
        (-1.0) * (*realGammaGph)["Gck"] * (*realGammaGpp)["Gda"] * (*Tph)["dk"],
      (*Lac)["ac"] +=
        (-1.0) * (*imagGammaGph)["Gck"] * (*imagGammaGpp)["Gda"] * (*Tph)["dk"],
      // Build Kki
      (*Kki)["ki"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Xabij)["cdil"],
      (*Kki)["ki"] +=  (-1.0) * (*Vpphh)["dckl"] * (*Xabij)["cdil"],
      // Build Lki
      (*Lki)["ki"] <<= (*Kki)["ki"],
      //TODO distinguish
      (*Lki)["ki"] += ( 2.0) * (*Vhhhp)["klic"] * (*Tph)["cl"],
      (*Lki)["ki"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["cl"],
      // Contract Lac with T2 Amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Lac)["ac"] * (*Tpphh)["cbij"],
      // Contract Lki with T2 Amplitudes
      (*Rpphh)["abij"] += (-1.0) * (*Lki)["ki"] * (*Tpphh)["abkj"],
      ///////////////////////////////////////
      // T2 Terms without Hirata intermediate
      ///////////////////////////////////////
      //TODO use an already introduced intermediate
      (*realDressedGammaGph)["Gai"] <<= (*realGammaGph)["Gai"],
      (*realDressedGammaGph)["Gai"] += (-1.0) * (*realGammaGhh)["Gki"] * (*Tph)["ak"],
      (*imagDressedGammaGph)["Gai"] <<= (*imagGammaGph)["Gai"],
      (*imagDressedGammaGph)["Gai"] += (-1.0) * (*imagGammaGhh)["Gki"] * (*Tph)["ak"],

      (*Rpphh)["abij"] +=
        ( 1.0) * (*realDressedGammaGph)["Gai"] * (*realGammaGpp)["Gbc"] * (*Tph)["cj"],
      (*Rpphh)["abij"] +=
        ( 1.0) * (*imagDressedGammaGph)["Gai"] * (*imagGammaGpp)["Gbc"] * (*Tph)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vhhhp)["jika"] * (*Tph)["bk"],
      (*Rpphh)["abij"] += (-1.0) * (*Tph)["bk"] * (*Vpphh)["acik"] * (*Tph)["cj"],
      ////////
      // Xakic
      ////////
      (*realDressedGammaGphTwo)["Gai"] <<= (*realGammaGph)["Gai"],
      (*imagDressedGammaGphTwo)["Gai"] <<= (*imagGammaGph)["Gai"],
      (*realDressedGammaGphTwo)["Gai"] += (-1.0) * (*realGammaGhh)["Gil"] * (*Tph)["al"],
      (*imagDressedGammaGphTwo)["Gai"] += (-1.0) * (*imagGammaGhh)["Gil"] * (*Tph)["al"],
      (*realDressedGammaGphTwo)["Gai"] += ( 1.0) * (*realGammaGpp)["Gad"] * (*Tph)["di"],
      (*imagDressedGammaGphTwo)["Gai"] += ( 1.0) * (*imagGammaGpp)["Gad"] * (*Tph)["di"],
      // FIXME: there is a better way for the contractions (see complex code)
      (*Xakic)["akic"] <<=
        ( 1.0) * (*realDressedGammaGphTwo)["Gai"] * (*realGammaGph)["Gck"],
      (*Xakic)["akic"] +=
        ( 1.0) * (*imagDressedGammaGphTwo)["Gai"] * (*imagGammaGph)["Gck"],


      (*Xakic)["akic"] += (-0.5) * (*Vpphh)["dclk"] * (*Yabij)["dail"],
      (*Xakic)["akic"] += ( 1.0) * (*Vpphh)["dclk"] * (*Tpphh)["adil"],
      //TODO if (!distinguishable) {
      (*Xakic)["akic"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["adil"],
      // Contract and Xakic intermediates with T2 amplitudes Tabij
      (*Zabij)["cbkj"] <<= ( 2.0) * (*Tpphh)["cbkj"],
      (*Zabij)["cbkj"] += (-1.0) * (*Tpphh)["bckj"],
      (*Rpphh)["abij"] += ( 1.0) * (*Xakic)["akic"] * (*Zabij)["cbkj"],


      ////////
      // Xakci
      ////////
      (*realDressedGammaGpp)["Gab"] <<= (*realGammaGpp)["Gab"],
      (*imagDressedGammaGpp)["Gab"] <<= (*imagGammaGpp)["Gab"],
      (*realDressedGammaGpp)["Gac"]  += (-1.0) * (*realGammaGph)["Gcl"] * (*Tph)["al"],
      (*imagDressedGammaGpp)["Gac"]  += (-1.0) * (*imagGammaGph)["Gcl"] * (*Tph)["al"],
      (*realDressedGammaGhh)["Gij"] <<= (*realGammaGhh)["Gij"],
      (*imagDressedGammaGhh)["Gij"] <<= (*imagGammaGhh)["Gij"],
      (*realDressedGammaGhh)["Gki"]  += ( 1.0) * (*realGammaGph)["Gdk"] * (*Tph)["di"],
      (*imagDressedGammaGhh)["Gki"]  += ( 1.0) * (*imagGammaGph)["Gdk"] * (*Tph)["di"],
      (*Xakci)["akci"] <<= (*realDressedGammaGpp)["Gac"] * (*realDressedGammaGhh)["Gki"],
      (*Xakci)["akci"] +=  (*imagDressedGammaGpp)["Gac"] * (*imagDressedGammaGhh)["Gki"],
      //TODO  if (!distinguishable) {
      (*Xakci)["akci"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["dail"],

      // Contract and Xakci intermediates with T2 amplitudes Tabij
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tpphh)["ackj"],

      // Symmetrize Rabij by applying permutation operator
      (*Rpphh)["abij"] += (*Rpphh)["baji"],

      // Add Vabij to Rabij (MP2 term)
      (*Rpphh)["abij"] += (*Vpphh)["abij"],
      ///////
      //Xklij
      ///////
      (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
      (*Xklij)["klij"]  += (*Vhhhp)["klic"] * (*Tph)["cj"],
      (*Xklij)["klij"]  += (*Vhhhp)["lkjc"] * (*Tph)["ci"],
      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
      (*Rpphh)["abij"]  += (*Xklij)["klij"] * (*Xabij)["abkl"],
      // Construct last term
      //TODO if (!distinguishable) {
      (*Xklij)["klij"] <<= (*Vpphh)["cdkl"] * (*Xabij)["cdij"],
      // TODO ppl
      (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Xabij)["abkl"],

      // Build Xabcd intermediate
      (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vphpp)["akcd"] * (*Tph)["bk"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vhppp)["kbcd"] * (*Tph)["ak"],

      // Contract Xabcd with T2 and T1 Amplitudes
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tpphh)["cdij"],
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tph)["ci"] * (*Tph)["dj"],
      //*******************************************************************
      //***********************  T1 amplitude equations  ******************
      //*******************************************************************

      // Contract Kac and Kki with T1 amplitudes
      (*Rph)["ai"] <<= ( 1.0) * (*Kac)["ac"] * (*Tph)["ci"],
      (*Rph)["ai"] += (-1.0) * (*Kki)["ki"] * (*Tph)["ak"],

      // Build Kck
      (*Kck)["ck"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Tph)["dl"],
      (*Kck)["ck"] +=  (-1.0) * (*Vpphh)["dckl"] * (*Tph)["dl"],

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rph)["ai"] += ( 2.0) * (*Kck)["ck"] * (*Tpphh)["caki"],
      (*Rph)["ai"] += (-1.0) * (*Kck)["ck"] * (*Tpphh)["caik"],
      (*Rph)["ai"] += ( 1.0) * (*Tph)["ak"] * (*Kck)["ck"] * (*Tph)["ci"],
      (*Rph)["ai"] += ( 2.0) * (*Vpphh)["acik"] * (*Tph)["ck"],
      //TODO replace Vphph
      (*Rph)["ai"] += (-1.0) * (*Vphph)["ciak"] * (*Tph)["ck"],

      (*Rph)["ai"] +=
        ( 2.0) * (*realGammaGpp)["Gca"] * (*realGammaGph)["Gdk"] * (*Xabij)["cdik"],
      (*Rph)["ai"] +=
        ( 2.0) * (*imagGammaGpp)["Gca"] * (*imagGammaGph)["Gdk"] * (*Xabij)["cdik"],
      (*Rph)["ai"] +=
        (-1.0) * (*realGammaGpp)["Gda"] * (*realGammaGph)["Gck"] * (*Xabij)["cdik"],
      (*Rph)["ai"] +=
        (-1.0) * (*imagGammaGpp)["Gda"] * (*imagGammaGph)["Gck"] * (*Xabij)["cdik"],

      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Xabij)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhhp)["lkic"] * (*Xabij)["ackl"]

    )->execute();
  }
  return residuum;
}


//template <typename Complex, typename TE>
//Ptr<FockVector<Complex,TE>> CcsdEnergyFromCoulombIntegrals::getResiduum(
//  const int iteration, const Ptr<const FockVector<Complex,TE>> &amplitudes
//) {
//
//  // get amplitude parts
//  auto Tph( amplitudes->get(0) );
//  auto Tpphh( amplitudes->get(1) );
//
//  // construct residuum
//  auto residuum( New<FockVector<Complex,TE>>(*amplitudes) );
//  *residuum *= 0.0;
//  auto Rph( residuum->get(0) );
//  auto Rpphh( residuum->get(1) );
//
//  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
//  auto coulombSlices(coulombIntegrals->getMap("slices"));
//  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<Complex,TE>>>("pphh"));
//
//  if (iteration == 0) {
//    //TODO
//    // && !isArgumentGiven("initialDoublesAmplitudes"))  {
//    // For first iteration compute only the MP2 amplitudes
//    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
//    COMPILE(
//      (*Rpphh)["abij"] += (*Vpphh)["abij"]
//    )->execute();
//  }
//
//
//  return residuum;
//}
