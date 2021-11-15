/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithms/CoupledClusterSinglesDoublesReference.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CoupledClusterSinglesDoublesReference)

Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>>
CoupledClusterSinglesDoublesReference::getResiduum(
  const int iteration,
  const Ptr<const TensorUnion<Real<>, DefaultDryTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DefaultDryTensorEngine>(iteration, amplitudes);
}

Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>>
CoupledClusterSinglesDoublesReference::getResiduum(
  const int iteration,
  const Ptr<const TensorUnion<Complex<>, DefaultDryTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DefaultDryTensorEngine>(iteration, amplitudes);
}

Ptr<TensorUnion<Real<>, DefaultTensorEngine>>
CoupledClusterSinglesDoublesReference::getResiduum(
  const int iteration,
  const Ptr<const TensorUnion<Real<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Real<>, DefaultTensorEngine>(iteration, amplitudes);
}


Ptr<TensorUnion<Complex<>, DefaultTensorEngine>>
CoupledClusterSinglesDoublesReference::getResiduum(
  const int iteration,
  const Ptr<const TensorUnion<Complex<>, DefaultTensorEngine>> &amplitudes
) {
  return getResiduum<Complex<>, DefaultTensorEngine>(iteration, amplitudes);
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

template <typename F, typename TE>
Ptr<TensorUnion<F,TE>> CoupledClusterSinglesDoublesReference::getResiduum(
  const int iteration, const Ptr<const TensorUnion<F,TE>> &amplitudes
) {
  // get amplitude parts
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  // construct residuum
  auto residuum( New<TensorUnion<F,TE>>(*amplitudes) );
  *residuum *= F(0);
  auto Rph( residuum->get(0) );
  auto Rpphh( residuum->get(1) );

  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));

  auto onlyPpl(arguments->getValue<size_t>("onlyPpl", 0) );


  if ( (iteration == 0) && !restart  && (onlyPpl == 0)) {

    //TODO
    // && !isArgumentGiven("initialDoublesAmplitudes"))  {
    // For first iteration compute only the MP2 amplitudes
//    OUT() << "MP2 T2 Amplitudes" << std::endl;
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  } else if (onlyPpl == 1) {

//    OUT() << "Calculate only PPL diagrams" << std::endl;

    auto Vpppp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pppp"));
    auto Vphpp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("phpp"));
    auto Vhppp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hppp"));
    auto Xabcd( Tcc<TE>::template tensor<F>("Xabcd") );

    COMPILE(
      // Build Xabcd intermediate
      (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vphpp)["akcd"] * (*Tph)["bk"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vhppp)["kbcd"] * (*Tph)["ak"],
      (*Rpphh)["abij"] <<= (*Xabcd)["abcd"] * (*Tpphh)["cdij"],
      (*Rpphh)["abij"]  += (*Xabcd)["abcd"] * (*Tph)["ci"] * (*Tph)["dj"]
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

//    OUT() << "Solving T2 Amplitude Equations" << std::endl;

    COMPILE(
      // Build Kac
      (*Kac)["ac"] <<= (-2.0) * (*Vhhpp)["klcd"] * (*Tpphh)["adkl"],
      (*Kac)["ac"] += ( 1.0) * (*Vhhpp)["kldc"] * (*Tpphh)["adkl"],
      (*Kac)["ac"] += (-2.0) * (*Vhhpp)["klcd"] * (*Tph)["ak"] * (*Tph)["dl"],
      (*Kac)["ac"] += ( 1.0) * (*Vhhpp)["kldc"] * (*Tph)["ak"] * (*Tph)["dl"],

      // Build Lac
      (*Lac)["ac"] <<= (*Kac)["ac"],
      (*Lac)["ac"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tph)["dk"],
      (*Lac)["ac"] += (-1.0) * (*Vphpp)["akdc"] * (*Tph)["dk"],

      // Build Kki
      (*Kki)["ki"] <<= (2.0) * (*Vhhpp)["klcd"] * (*Tpphh)["cdil"],
      (*Kki)["ki"] += (-1.0) * (*Vhhpp)["kldc"] * (*Tpphh)["cdil"],
      (*Kki)["ki"] += ( 2.0) * (*Vhhpp)["klcd"] * (*Tph)["ci"] * (*Tph)["dl"],
      (*Kki)["ki"] += (-1.0) * (*Vhhpp)["kldc"] * (*Tph)["ci"] * (*Tph)["dl"],

      // Build Lki
      (*Lki)["ki"] <<= (*Kki)["ki"],
      (*Lki)["ki"] += ( 2.0) * (*Vhhhp)["klic"] * (*Tph)["cl"],
      (*Lki)["ki"] += (-1.0) * (*Vhhph)["klci"] * (*Tph)["cl"],

      // Contract Lac with T2 Amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Lac)["ac"] * (*Tpphh)["cbij"],

      // Contract Lki with T2 Amplitudes
      (*Rpphh)["abij"] += (-1.0) * (*Lki)["ki"] * (*Tpphh)["abkj"],

      // Contract Coulomb integrals with T2 amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Vpphp)["abic"] * (*Tph)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vhphp)["kbic"] * (*Tph)["ak"] * (*Tph)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhh)["akij"] * (*Tph)["bk"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhp)["akic"] * (*Tph)["cj"] * (*Tph)["bk"],

      // Build Xakic
      (*Xakic)["akic"] <<= (*Vphhp)["akic"],
      (*Xakic)["akic"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["al"],
      (*Xakic)["akic"] += ( 1.0) * (*Vphpp)["akdc"] * (*Tph)["di"],
      (*Xakic)["akic"] += (-0.5) * (*Vhhpp)["lkdc"] * (*Tpphh)["dail"],
      (*Xakic)["akic"] += (-1.0) * (*Vhhpp)["lkdc"] * (*Tph)["di"] * (*Tph)["al"],
      (*Xakic)["akic"] += ( 1.0) * (*Vhhpp)["lkdc"] * (*Tpphh)["adil"],
      (*Xakic)["akic"] += (-0.5) * (*Vhhpp)["lkcd"] * (*Tpphh)["adil"],
      (*Rpphh)["abij"] += ( 2.0) * (*Xakic)["akic"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakic)["akic"] * (*Tpphh)["bckj"],

      // Build Xakci
      (*Xakci)["akci"] <<= (*Vphph)["akci"],
      (*Xakci)["akci"] += (-1.0) * (*Vhhph)["lkci"] * (*Tph)["al"],
      (*Xakci)["akci"] += ( 1.0) * (*Vphpp)["akcd"] * (*Tph)["di"],
      (*Xakci)["akci"] += (-0.5) * (*Vhhpp)["lkcd"] * (*Tpphh)["dail"],
      (*Xakci)["akci"] += (-1.0) * (*Vhhpp)["lkcd"] * (*Tph)["di"] * (*Tph)["al"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tpphh)["ackj"],

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
      (*Xklij)["klij"] += (*Vhhhp)["klic"] * (*Tph)["cj"],
      (*Xklij)["klij"] += (*Vhhph)["klcj"] * (*Tph)["ci"],
      (*Xklij)["klij"] += (*Vhhpp)["klcd"] * (*Tpphh)["cdij"],
      (*Xklij)["klij"] += (*Vhhpp)["klcd"] * (*Tph)["ci"] * (*Tph)["dj"],

      // Contract Xklij with T2 Amplitudes
      (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tpphh)["abkl"],

      // Contract Xklij with T1 Amplitudes
      (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tph)["ak"] * (*Tph)["bl"],

      // Build Xabcd intermediate
      (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vphpp)["akcd"] * (*Tph)["bk"],
      (*Xabcd)["abcd"] += (-1.0) * (*Vhppp)["kbcd"] * (*Tph)["ak"],

      // Contract Xabcd with T2 and T1 Amplitudes
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tpphh)["cdij"],
      (*Rpphh)["abij"] += (*Xabcd)["abcd"] * (*Tph)["ci"] * (*Tph)["dj"]
    )->execute();

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
//    OUT() << "Solving T1 Amplitude Equations" << std::endl;
    COMPILE(
      // Contract Kac and Kki with T1 amplitudes
      (*Rph)["ai"] += ( 1.0) * (*Kac)["ac"] * (*Tph)["ci"],
      (*Rph)["ai"] += (-1.0) * (*Kki)["ki"] * (*Tph)["ak"],

      //Build Kck
      (*Kck)["ck"] <<= ( 2.0) * (*Vhhpp)["klcd"] * (*Tph)["dl"],
      (*Kck)["ck"]  += (-1.0) * (*Vhhpp)["kldc"] * (*Tph)["dl"],

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rph)["ai"] += ( 2.0) * (*Kck)["ck"] * (*Tpphh)["caki"],
      (*Rph)["ai"] += (-1.0) * (*Kck)["ck"] * (*Tpphh)["caik"],
      (*Rph)["ai"] += ( 1.0) * (*Kck)["ck"] * (*Tph)["ci"] * (*Tph)["ak"],
      (*Rph)["ai"] += ( 2.0) * (*Vphhp)["akic"] * (*Tph)["ck"],
      (*Rph)["ai"] += (-1.0) * (*Vphph)["akci"] * (*Tph)["ck"],
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tpphh)["cdik"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tpphh)["cdik"],
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tph)["ci"] * (*Tph)["dk"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tph)["ci"] * (*Tph)["dk"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tpphh)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tpphh)["ackl"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tph)["ak"] * (*Tph)["cl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tph)["ak"] * (*Tph)["cl"]
    )->execute();

  }
  return residuum;
}
