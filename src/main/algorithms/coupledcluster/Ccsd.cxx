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

// TODO: discuss 2cc and dcsd options
// TODO: cleanup TODOs

#include <algorithms/coupledcluster/Ccsd.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

template <typename TE>
CoupledClusterMethodRegistrar<
  Real<>,TE,Ccsd<Real<>,TE>
> Ccsd<Real<>,TE>::registrar_("Ccsd");
template <typename TE>
CoupledClusterMethodRegistrar<
  Complex<>,TE,Ccsd<Complex<>,TE>
> Ccsd<Complex<>,TE>::registrar_("Ccsd");

template <typename TE>
std::string Ccsd<Real<>,TE>::describeOptions() {
  std::stringstream stream;
  auto eigenEnergies(this->arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->template getPtr<TensorExpression<Real<>,TE>>("h"));
  auto No(epsi->inspect()->getLen(0));
  auto integralsSliceSize(
    this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
  );
  stream
    << "integralsSliceSize: " << integralsSliceSize;
  return stream.str();
}

template <typename TE>
std::string Ccsd<Complex<>,TE>::describeOptions() {
  std::stringstream stream;
  auto eigenEnergies(this->arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->template getPtr<TensorExpression<Real<>,TE>>("h"));
  auto No(epsi->inspect()->getLen(0));
  auto integralsSliceSize(
    this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
  );
  stream
    << "integralsSliceSize: " << integralsSliceSize;
  return stream.str();
}


//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

template <typename TE>
Ptr<TensorUnion<Real<>,TE>> Ccsd<Real<>,TE>::getResiduum(
  const Ptr<TensorUnion<Real<>,TE>> &amplitudes
) {
  // construct residuum. Shape will be assumed upon first use.
  auto Rph( Tcc<TE>::template tensor<Real<>>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<Real<>>("Rpphh") );
  auto residuum(
    New<TensorUnion<Real<>,TE>>(
      std::vector<Ptr<TensorExpression<Real<>,TE>>>({Rph, Rpphh}),
      std::vector<std::string>({"ai", "abij"})
    )
  );

  auto coulombIntegrals(this->arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->template getPtr<TensorExpression<Real<>,TE>>("pphh"));
  bool ppl(this->arguments->template getValue<bool>("ppl", true));
  bool twoCc(this->arguments->template getValue<bool>("2cc", false));
  bool dcsd(this->arguments->template getValue<bool>("dcsd", false));

  ASSERT_LOCATION(
    ((ppl && !(twoCc && dcsd)) || (!ppl && !twoCc && !dcsd)),
    "!ppl, 2cc and dcsd are all mutually exclusive.",
    this->arguments->sourceLocation
  );

  if (!amplitudes) {
    // no previous amplitudes given
    COMPILE(
      (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= (*Vpphh)["abij"]
    )->execute();
  } else {
    // TODO: check if given amplitudes contain expected parts
    // get amplitude parts
    auto Tph( amplitudes->get(0) );
    auto Tpphh( amplitudes->get(1) );
    Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");

//    OUT() << "\tSolving T2 Amplitude Equations" << std::endl;
    auto slicedCoulombVertex(this->arguments->getMap("slicedCoulombVertex"));
    auto slices(slicedCoulombVertex->getMap("slices"));
    auto orbitals(coulombIntegrals->template getValue<std::string>("scalarType"));
    auto GammaGpp(slices->template getPtr<TensorExpression<Complex<>,TE>>("pp"));
    auto GammaGph(slices->template getPtr<TensorExpression<Complex<>,TE>>("ph"));
    auto GammaGhh(slices->template getPtr<TensorExpression<Complex<>,TE>>("hh"));

    auto Vphph(coulombSlices->template getPtr<TensorExpression<Real<>,TE>>("phph"));
    auto Vhhhh(coulombSlices->template getPtr<TensorExpression<Real<>,TE>>("hhhh"));
    auto Vhhhp(coulombSlices->template getPtr<TensorExpression<Real<>,TE>>("hhhp"));

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
    auto realDressedGammaGpp(
      Tcc<TE>::template tensor<Real<>>("realDressedGammaGpp")
    );
    auto imagDressedGammaGpp(
      Tcc<TE>::template tensor<Real<>>("imagDressedGammaGpp")
    );
    auto realDressedGammaGph(
      Tcc<TE>::template tensor<Real<>>("realDressedGammaGph")
    );
    auto imagDressedGammaGph(
      Tcc<TE>::template tensor<Real<>>("imagDressedGammaGph")
    );
    auto realDressedGammaGhh(
      Tcc<TE>::template tensor<Real<>>("realDressedGammaGhh")
    );
    auto imagDressedGammaGhh( Tcc<TE>::template tensor<Real<>>("imagDressedGammaGhh") );
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
      (*Yabij)["abij"] += ( 2.0) * (*Tph)["ai"] * (*Tph)["bj"]
    )->execute();


    /////////////////////////////////////
    // Lac and Kac for doubles amplitudes
    /////////////////////////////////////
    if(twoCc) {
 //     OUT() << "\tUsing the 2CC approximation"  << std::endl;
      COMPILE(
        // Build Kki, leave A-term unchanged
        (*Kki)["ki"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Xabij)["cdil"],
        (*Kki)["ki"] +=  (-1.0) * (*Vpphh)["dckl"] * (*Xabij)["cdil"],
        // Build Lki
        (*Lki)["ki"] <<= (*Kki)["ki"],
        // Build Kac with C-term removed (https://doi.org/10.1063/1.4979078, figure 1)
        (*Kac)["ac"] <<= (-2.0) * (*Vpphh)["cdkl"] * (*Tph)["ak"] * (*Tph)["dl"],
        (*Kac)["ac"]  += ( 1.0) * (*Vpphh)["dckl"] * (*Tph)["ak"] * (*Tph)["dl"],
        // Lac
        (*Lac)["ac"] <<= (*Kac)["ac"],
        //Add left-out term of Kac for T1-equation later
        (*Kac)["ac"]  += (-2.0) * (*Vpphh)["cdkl"] * (*Tpphh)["adkl"],
        (*Kac)["ac"]  += ( 1.0) * (*Vpphh)["dckl"] * (*Tpphh)["adkl"],
        //Xakic with D-term removed
        (*Xakic)["akic"] <<= (-1.0) * (*Vpphh)["dclk"] * (*Tph)["di"] * (*Tph)["al"]
        //Leave out D^{ex}-term from Xakci (i.e don't do anything at this point)
      )->execute();
    }
    else if(dcsd){
//      OUT() << "\tUsing the DCSD approximation"  << std::endl;
      COMPILE(
        // Build Kki, prefactor A-term by 0.5
        (*Kki)["ki"] <<= ( 1.0) * (*Vpphh)["cdkl"] * (*Yabij)["cdil"],
        (*Kki)["ki"] +=  (-0.5) * (*Vpphh)["dckl"] * (*Yabij)["cdil"],
        // Build Lki
        (*Lki)["ki"] <<= (*Kki)["ki"],
        //Add left out term of Kki for T1-equation later
        (*Kki)["ki"] += ( 1.0) * (*Vpphh)["cdkl"] * (*Tpphh)["cdil"],
        (*Kki)["ki"] += (-0.5) * (*Vpphh)["dckl"] * (*Tpphh)["cdil"],
        // Build Kac, prefactor C-term by 0.5
        (*Kac)["ac"] <<= (-1.0) * (*Vpphh)["cdkl"] * (*Yabij)["adkl"],
        (*Kac)["ac"]  += ( 0.5) * (*Vpphh)["dckl"] * (*Yabij)["adkl"],
        // Lac
        (*Lac)["ac"] <<= (*Kac)["ac"],
        //Add left-out term of Kac for T1-equation later
        (*Kac)["ac"] += (-1.0) * (*Vpphh)["cdkl"] * (*Tpphh)["adkl"],
        (*Kac)["ac"] += (0.5) * (*Vpphh)["dckl"] * (*Tpphh)["adkl"],
        //Xakic with D^{ex}-term removed
        (*Xakic)["akic"] <<= (-0.5) * (*Vpphh)["dclk"] * (*Yabij)["dail"],
        (*Xakic)["akic"] += ( 1.0) * (*Vpphh)["dclk"] * (*Tpphh)["adil"]
        //Leave out D^{ex}-term from Xakci (i.e dont't do anything at this point)
      )->execute();
    }
    else{
//      OUT() << "\tNOT using the 2CC approximation"  << std::endl;
      COMPILE(
        // Build Kki
        (*Kki)["ki"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Xabij)["cdil"],
        (*Kki)["ki"] +=  (-1.0) * (*Vpphh)["dckl"] * (*Xabij)["cdil"],
        // Build Lki
        (*Lki)["ki"] <<= (*Kki)["ki"],
        //Build Kac
        (*Kac)["ac"] <<= (-2.0) * (*Vpphh)["cdkl"] * (*Xabij)["adkl"],
        (*Kac)["ac"]  += ( 1.0) * (*Vpphh)["dckl"] * (*Xabij)["adkl"],
        // Lac
        (*Lac)["ac"] <<= (*Kac)["ac"],
        //Xakic
        (*Xakic)["akic"] <<= (-0.5) * (*Vpphh)["dclk"] * (*Yabij)["dail"],
        (*Xakic)["akic"] += ( 1.0) * (*Vpphh)["dclk"] * (*Tpphh)["adil"], //if 2CC: throw out (D^c term)
        //TODO if (!distinguishable) {
        (*Xakic)["akic"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["adil"], //if 2CC: throw out (D^{ex} term)
        //Xakci
        //TODO  if (!distinguishable) {
        (*Xakci)["akci"] <<= (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["dail"]
      )->execute();
    }

    COMPILE(

      // TODO distiguish
      (*Lac)["ac"] +=
        ( 2.0) * (*realGammaGpp)["Gca"] * (*realGammaGph)["Gdk"] * (*Tph)["dk"],
      (*Lac)["ac"] +=
        ( 2.0) * (*imagGammaGpp)["Gca"] * (*imagGammaGph)["Gdk"] * (*Tph)["dk"],
      (*Lac)["ac"] +=
        (-1.0) * (*realGammaGph)["Gck"] * (*realGammaGpp)["Gda"] * (*Tph)["dk"],
      (*Lac)["ac"] +=
        (-1.0) * (*imagGammaGph)["Gck"] * (*imagGammaGpp)["Gda"] * (*Tph)["dk"],


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
      // Add further dressing to dressed Vertex
      (*realDressedGammaGph)["Gai"] += ( 1.0) * (*realGammaGpp)["Gad"] * (*Tph)["di"],
      (*imagDressedGammaGph)["Gai"] += ( 1.0) * (*imagGammaGpp)["Gad"] * (*Tph)["di"],
      // FIXME: there is a better way for the contractions (see complex code)
      (*Xakic)["akic"] +=
        ( 1.0) * (*realDressedGammaGph)["Gai"] * (*realGammaGph)["Gck"],
      (*Xakic)["akic"] +=
        ( 1.0) * (*imagDressedGammaGph)["Gai"] * (*imagGammaGph)["Gck"],

      // Contract and Xakic intermediates with T2 amplitudes Tabij
      (*Zabij)["cbkj"] <<= ( 2.0) * (*Tpphh)["cbkj"],
      (*Zabij)["cbkj"] += (-1.0) * (*Tpphh)["bckj"],
      (*Rpphh)["abij"] += ( 1.0) * (*Xakic)["akic"] * (*Zabij)["cbkj"], //affected by 2CC
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
      (*Xakci)["akci"] += (*realDressedGammaGpp)["Gac"] * (*realDressedGammaGhh)["Gki"],
      (*Xakci)["akci"] +=  (*imagDressedGammaGpp)["Gac"] * (*imagDressedGammaGhh)["Gki"],


      // Contract and Xakci intermediates with T2 amplitudes Tabij
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tpphh)["cbkj"], //affected by 2CC
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tpphh)["ackj"], //affected by 2CC

      // Symmetrize Rabij by applying permutation operator
      (*Rpphh)["abij"] += (*Rpphh)["baji"]
    )->execute();

      // Add Vabij to Rabij (MP2 term)
      if (ppl) {
        COMPILE(
          (*Rpphh)["abij"] += (*Vpphh)["abij"]
        )->execute();
      }
    if(dcsd && ppl){
//      OUT() << "\tdcsd AND ppl are true"  << std::endl;
      COMPILE(
        ///////
        //Xklij
        ///////
        (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
        (*Xklij)["klij"]  += (*Vhhhp)["klic"] * (*Tph)["cj"],
        (*Xklij)["klij"]  += (*Vhhhp)["lkjc"] * (*Tph)["ci"],
        // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
        (*Rpphh)["abij"]  += (*Xklij)["klij"] * (*Xabij)["abkl"],
        //Construct last term without B-term
        (*Xklij)["klij"] <<= (*Vpphh)["cdkl"] * (*Tph)["ci"] * (*Tph)["dj"],
        (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tpphh)["abkl"]//,
//        //Add left-out term for second Xklij-contribution
//        (*Xklij)["klij"] += (*Vpphh)["cdkl"] * (*Tpphh)["cdij"],
//        (*Rpphh)["abij"] += (*Xklij)["klij"] * (*Tph)["ak"] * (*Tph)["bl"]
      )->execute();
    }
    else{
      COMPILE(
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
        (ppl) ? (
          (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Tpphh)["abkl"]
        ) : (
          (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Xabij)["abkl"]
        )
      )->execute();
    }

//    COMPILE(
//      ///////
//      //Xklij
//      ///////
//      (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
//      (*Xklij)["klij"]  += (*Vhhhp)["klic"] * (*Tph)["cj"],
//      (*Xklij)["klij"]  += (*Vhhhp)["lkjc"] * (*Tph)["ci"],
//      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
//      (*Rpphh)["abij"]  += (*Xklij)["klij"] * (*Xabij)["abkl"],
//      // Construct last term
//      //TODO if (!distinguishable) {
//      (*Xklij)["klij"] <<= (*Vpphh)["cdkl"] * (*Xabij)["cdij"],
//      (ppl) ? (
//        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Tpphh)["abkl"]
//      ) : (
//        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Xabij)["abkl"]
//      )
//    )->execute();

    if (ppl) {
//      OUT() << "\tAdding Particle-particle contraction"  << std::endl;
      Natural<> Nv(realDressedGammaGpp->lens[1]);
      Natural<> NG(realDressedGammaGpp->lens[0]);
      Natural<> No(Rpphh->inspect()->getLen(2));
      Natural<> sliceSize(
        this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
      );
      Natural<> numberSlices(Natural<>(ceil(1.0*Nv/sliceSize)));
      std::vector<Ptr<Tensor<Real<>, TE>>> realSlicedGammaGpp;
      std::vector<Ptr<Tensor<Real<>, TE>>> imagSlicedGammaGpp;
      //Slice GammaGab and store it in a vector
      for (Natural<> v(0); v < numberSlices; v++){
        Natural<> xStart = v*sliceSize;
        Natural<> xEnd = std::min((v+1)*sliceSize,Nv);
        auto dummyr( Tcc<TE>::template tensor<Real<>>("dummyr") );
        auto dummyi( Tcc<TE>::template tensor<Real<>>("dummyi") );
        COMPILE(
          (*dummyr)["Gxb"] <<=
            (*(*realDressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"],
          (*dummyi)["Gxb"] <<=
            (*(*imagDressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"]
        )->execute();
        realSlicedGammaGpp.push_back(dummyr);
        imagSlicedGammaGpp.push_back(dummyi);
      }
      // loop over slices
      for (Natural<> m(0); m < numberSlices; m++)
      for (Natural<> n(m); n < numberSlices; n++){
        auto Vxycd( Tcc<TE>::template tensor<Real<>>("Vxycd") );
        auto Rxyij( Tcc<TE>::template tensor<Real<>>("Rxyij") );
        auto Ryxji( Tcc<TE>::template tensor<Real<>>("Ryxji") );
        Natural<> a(n*sliceSize); Natural<> b(m*sliceSize);
        Natural<> Nx(realSlicedGammaGpp[n]->lens[1]);
        Natural<> Ny(realSlicedGammaGpp[m]->lens[1]);
        COMPILE(
          (*Vxycd)["xycd"] <<=
            (*realSlicedGammaGpp[n])["Gxc"] * (*realSlicedGammaGpp[m])["Gyd"],
          (*Vxycd)["xycd"]  +=
            (*imagSlicedGammaGpp[n])["Gxc"] * (*imagSlicedGammaGpp[m])["Gyd"],
          (*Rxyij)["xyij"] <<= (*Vxycd)["xycd"] * (*Xabij)["cdij"],
          (*(*Rpphh)({a, b, 0, 0},{a+Nx, b+Ny, No, No}))["xyij"] += (*Rxyij)["xyij"],
          // if a>b: add the same slice at (b,a,j,i):
          (a>b) ? (
            (*Ryxji)["yxji"] <<= (*Rxyij)["xyij"],
            (*(*Rpphh)({b, a, 0, 0},{b+Ny, a+Nx, No, No}))["xyij"] += (*Ryxji)["xyij"]
          ) : (
            Tcc<TE>::sequence()
          )
        )->execute();
      }
    }

    //*******************************************************************
    //***********************  T1 amplitude equations  ******************
    //*******************************************************************
//    OUT() << "\tSolving T1 Amplitude Equations" << std::endl;
    COMPILE(
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
      //TODO maybe replace Vphph
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


template <typename TE>
Ptr<TensorUnion<Complex<>,TE>> Ccsd<Complex<>,TE>::getResiduum(
  const Ptr<TensorUnion<Complex<>,TE>> &amplitudes
) {
  // construct residuum. Shape will be assumed upon first use.
  auto Rph( Tcc<TE>::template tensor<Complex<>>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<Complex<>>("Rpphh") );
  auto residuum(
    New<TensorUnion<Complex<>,TE>>(
      std::vector<Ptr<TensorExpression<Complex<>,TE>>>({Rph, Rpphh}),
      std::vector<std::string>({"ai", "abij"})
    )
  );

  auto coulombIntegrals(this->arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("pphh"));
  bool ppl(this->arguments->template getValue<bool>("ppl", true));

  if (!amplitudes) {
    // no previous amplitudes given
    COMPILE(
      (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= (*Vpphh)["abij"]
    )->execute();
  } else {
    // TODO: check if given amplitudes contain expected parts
    // get amplitude parts
    auto Tph( amplitudes->get(0) );
    auto Tpphh( amplitudes->get(1) );
    Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");

//    OUT() << "\tSolving T2 Amplitude Equations" << std::endl;
    auto slicedCoulombVertex(this->arguments->getMap("slicedCoulombVertex"));
    auto slices(slicedCoulombVertex->getMap("slices"));
    auto orbitals(coulombIntegrals->template getValue<std::string>("scalarType"));
    auto GammaGpp(slices->template getPtr<TensorExpression<Complex<>,TE>>("pp"));
    auto GammaGph(slices->template getPtr<TensorExpression<Complex<>,TE>>("ph"));
    auto GammaGhp(slices->template getPtr<TensorExpression<Complex<>,TE>>("hp"));
    auto GammaGhh(slices->template getPtr<TensorExpression<Complex<>,TE>>("hh"));
//we need all conjugate transposed.
//strange labeling, though
    auto cTGammaGph( Tcc<TE>::template tensor<Complex<>>("cTGammaGph"));
    auto cTGammaGhp( Tcc<TE>::template tensor<Complex<>>("cTGammaGhp"));
    auto cTGammaGpp( Tcc<TE>::template tensor<Complex<>>("cTGammaGpp"));
    auto cTGammaGhh( Tcc<TE>::template tensor<Complex<>>("cTGammaGhh"));
    COMPILE(
      (*cTGammaGpp)["Gab"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGpp)["Gba"]),
      (*cTGammaGhp)["Gia"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGph)["Gai"]),
      (*cTGammaGph)["Gai"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGhp)["Gia"]),
      (*cTGammaGhh)["Gij"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGhh)["Gji"])
    )->execute();
    auto cTDressedGammaGph( Tcc<TE>::template tensor<Complex<>>("cTDressedGammaGph"));
    auto cTDressedGammaGpp( Tcc<TE>::template tensor<Complex<>>("cTDressedGammaGpp"));
    auto dressedGammaGhh( Tcc<TE>::template tensor<Complex<>>("dressedGammaGhh"));
    auto dressedGammaGpp( Tcc<TE>::template tensor<Complex<>>("dressedGammaGpp"));

    auto Vphhp(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("phhp"));
    auto Vhhpp(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("hhpp"));
    auto Vphph(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("phph"));
    auto Vhhhh(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("hhhh"));
    auto Vhhhp(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("hhhp"));
    auto Vphhh(coulombSlices->template getPtr<TensorExpression<Complex<>,TE>>("phhh"));
    // Hirata intermediates
    auto Lac( Tcc<TE>::template tensor<Complex<>>("Lac") );
    auto Kac( Tcc<TE>::template tensor<Complex<>>("Kac") );
    auto Lki( Tcc<TE>::template tensor<Complex<>>("Lki") );
    auto Kki( Tcc<TE>::template tensor<Complex<>>("Kki") );
    auto Kck( Tcc<TE>::template tensor<Complex<>>("Kck") );
    auto Xabij( Tcc<TE>::template tensor<Complex<>>("Xabij") ); // T2+T1*T1
    auto Yabij( Tcc<TE>::template tensor<Complex<>>("Yabij") ); // T2+2*T1*T1
    auto Zabij( Tcc<TE>::template tensor<Complex<>>("Zabij") );
    auto Xklij( Tcc<TE>::template tensor<Complex<>>("Xklij") );
    auto Xakci( Tcc<TE>::template tensor<Complex<>>("Xakci") );
    auto Xakic( Tcc<TE>::template tensor<Complex<>>("Xakic") );
    auto Xabcd( Tcc<TE>::template tensor<Complex<>>("Xabcd") );

    COMPILE(
      (*Xabij)["abij"] <<= (*Tpphh)["abij"],
      (*Xabij)["abij"] += (*Tph)["ai"] * (*Tph)["bj"],
      (*Yabij)["abij"] <<= (*Tpphh)["abij"],
      (*Yabij)["abij"] += ( 2.0) * (*Tph)["ai"] * (*Tph)["bj"],
      // Build Kac
      (*Kac)["ac"] <<= (-2.0) * (*Vhhpp)["klcd"] * (*Xabij)["adkl"],
      (*Kac)["ac"] += ( 1.0) * (*Vhhpp)["kldc"] * (*Xabij)["adkl"],
      // Build Lac
      (*Lac)["ac"] <<= (*Kac)["ac"],
      (*Lac)["ac"] +=  (2.0) * (*cTGammaGpp)["Gac"] * (*GammaGhp)["Gkd"] * (*Tph)["dk"],
      (*Lac)["ac"] += (-1.0) * (*cTGammaGpp)["Gad"] * (*GammaGhp)["Gkc"] * (*Tph)["dk"],
      // Build Kki
      (*Kki)["ki"] <<= (2.0) * (*Vhhpp)["klcd"] * (*Xabij)["cdil"],
      (*Kki)["ki"] += (-1.0) * (*Vhhpp)["kldc"] * (*Xabij)["cdil"],

      // Build Lki
      (*Lki)["ki"] <<= (*Kki)["ki"],
      (*Lki)["ki"] += ( 2.0) * (*Vhhhp)["klic"] * (*Tph)["cl"],
      (*Lki)["ki"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["cl"],

      // Contract Lac with T2 Amplitudes
      (*Rpphh)["abij"] += ( 1.0) * (*Lac)["ac"] * (*Tpphh)["cbij"],

      // Contract Lki with T2 Amplitudes
      (*Rpphh)["abij"] += (-1.0) * (*Lki)["ki"] * (*Tpphh)["abkj"],

      // Contract Coulomb integrals with T1 amplitudes
      (*cTDressedGammaGph)["Gai"] <<= (*cTGammaGph)["Gai"],
      (*cTDressedGammaGph)["Gai"] += (-1.0) * (*cTGammaGhh)["Gki"] * (*Tph)["ak"],
      (*Rpphh)["abij"] += (*cTDressedGammaGph)["Gai"] * (*GammaGpp)["Gbc"] * (*Tph)["cj"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhh)["akij"] * (*Tph)["bk"],
      (*Rpphh)["abij"] += (-1.0) * (*Vphhp)["akic"] * (*Tph)["cj"] * (*Tph)["bk"],

      // Build Xakic
      (*cTDressedGammaGph)["Gai"] += ( 1.0) * (*cTGammaGpp)["Gad"] * (*Tph)["di"],
      (*cTDressedGammaGph)["Gai"] += (-0.5) * (*cTGammaGhp)["Gld"] * (*Yabij)["dail"],
      (*Xakic)["akic"] <<= (*cTDressedGammaGph)["Gai"] * (*GammaGhp)["Gkc"],
      (*Yabij)["dclk"] <<= (1.0) * (*Vhhpp)["lkdc"],
      (*Yabij)["dclk"] += (-0.5) * (*Vhhpp)["lkcd"],
      (*Xakic)["akic"] += (*Yabij)["dclk"] * (*Tpphh)["adil"],
      (*Yabij)["cbkj"] <<= (2.0) * (*Tpphh)["cbkj"],
      (*Yabij)["cbkj"] += (-1.0) * (*Tpphh)["bckj"],
      (*Rpphh)["abij"] += (*Xakic)["akic"] * (*Yabij)["cbkj"],

      // Build Xakci
      (*cTDressedGammaGpp)["Gab"] <<= (*cTGammaGpp)["Gab"],
      (*cTDressedGammaGpp)["Gac"] += (-1.0) * (*cTGammaGhp)["Glc"] * (*Tph)["al"],
      (*dressedGammaGhh)["Gij"] <<= (*GammaGhh)["Gij"],
      (*dressedGammaGhh)["Gki"]  += (*GammaGhp)["Gkd"] * (*Tph)["di"],
      // Xakci = Vakci - Vlkci * Tal + Vakcd * Tdi - Vlkcd * Tdail
      (*Xakci)["akci"] <<= (*cTDressedGammaGpp)["Gac"] * (*dressedGammaGhh)["Gki"],
      (*Xakci)["akci"] += (-0.5) * (*Vhhpp)["lkcd"] * (*Tpphh)["dail"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tpphh)["ackj"],

      // Symmetrize Rpphh by applying permutation operator
      (*Rpphh)["abij"] += (*Rpphh)["baji"]

    )->execute();

    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rpphh that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////
    if (ppl) {
      COMPILE(
        (*Rpphh)["abij"] += (*Vpphh)["abij"]
      )->execute();
    }
    COMPILE(
      // Build Xklij intermediate
      (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
      (*Xklij)["klij"]  += (*Vhhhp)["klic"] * (*Tph)["cj"],
      (*Xklij)["klij"]  += (*Vhhhp)["lkjc"] * (*Tph)["ci"],
      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
      (*Rpphh)["abij"]  += (*Xklij)["klij"] * (*Xabij)["abkl"],
      // Construct last term
      //TODO if (!distinguishable) {
      (*Xklij)["klij"] <<= (*Vhhpp)["klcd"] * (*Xabij)["cdij"],
      (ppl) ? (
        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Tpphh)["abkl"]
      ) : (
        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Xabij)["abkl"]
      )
    )->execute();

    if (ppl) {
//      OUT() << "\tAdding Particle-particle contraction"  << std::endl;
      Natural<> NG(cTGammaGpp->getLen(0));
      Natural<> Nv(Rpphh->inspect()->getLen(0));
      Natural<> No(Rpphh->inspect()->getLen(2));
      Natural<> sliceSize(
        this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
      );
      Natural<> numberSlices(Natural<>(ceil(1.0*Nv/sliceSize)));
      std::vector<Ptr<Tensor<Complex<>, TE>>> cTSlicedGammaGpp;
      std::vector<Ptr<Tensor<Complex<>, TE>>>   SlicedGammaGpp;
      COMPILE(
        (*dressedGammaGpp)["Gab"] <<= (*GammaGpp)["Gab"],
        (*dressedGammaGpp)["Gab"] += (-1.0) * (*GammaGhp)["Gkb"] * (*Tph)["ak"]
      )->execute();
      //Slice GammaGab and store it in a vector
      for (Natural<> v(0); v < numberSlices; v++){
        Natural<> xStart = v*sliceSize;
        Natural<> xEnd = std::min((v+1)*sliceSize,Nv);
        auto dummy(   Tcc<TE>::template tensor<Complex<>>("dummy")   );
        auto dummyct( Tcc<TE>::template tensor<Complex<>>("dummyct") );
        COMPILE(
          (*dummy )["Gxb"]  <<=
            (*(*dressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"],
          (*dummyct)["Gxb"] <<=
            (*(*cTDressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"]
        )->execute();
          SlicedGammaGpp.push_back(dummy);
        cTSlicedGammaGpp.push_back(dummyct);
      }
      // loop over slices
      for (Natural<> m(0); m < numberSlices; m++)
      for (Natural<> n(m); n < numberSlices; n++){
        auto Vxycd( Tcc<TE>::template tensor<Complex<>>("Vxycd") );
        auto Rxyij( Tcc<TE>::template tensor<Complex<>>("Rxyij") );
        auto Ryxji( Tcc<TE>::template tensor<Complex<>>("Ryxji") );
        Natural<> a(n*sliceSize); Natural<> b(m*sliceSize);
        Natural<> Nx(cTSlicedGammaGpp[n]->lens[1]);
        Natural<> Ny(SlicedGammaGpp[m]->lens[1]);
        COMPILE(
          (*Vxycd)["xycd"] <<=
            (*cTSlicedGammaGpp[n])["Gxc"] * (*SlicedGammaGpp[m])["Gyd"],
          (*Rxyij)["xyij"] <<= (*Vxycd)["xycd"] * (*Xabij)["cdij"],
          (*(*Rpphh)({a, b, 0, 0},{a+Nx, b+Ny, No, No}))["xyij"] += (*Rxyij)["xyij"],
          // if a>b: add the same slice at (b,a,j,i):
          (a>b) ? (
            (*Ryxji)["yxji"] <<= (*Rxyij)["xyij"],
            (*(*Rpphh)({b, a, 0, 0},{b+Ny, a+Nx, No, No}))["xyij"] += (*Ryxji)["xyij"]
          ) : (
            Tcc<TE>::sequence()
          )
        )->execute();
      }
//      COMPILE(
//        (*dressedGammaGpp)["Gab"] <<= (*GammaGpp)["Gab"],
//        (*dressedGammaGpp)["Gab"] += (-1.0) * (*GammaGhp)["Gkb"] * (*Tph)["ak"],
//        (*Vxycd)["xycd"] <<= (*cTDressedGammaGpp)["Gxc"] * (*dressedGammaGpp)["Gyd"],
//        (*Rxyij)["xyij"] <<= (*Vxycd)["xycd"] * (*Xabij)["cdij"],
//        (*Rpphh)["abij"] += (*Rxyij)["abij"]
//      )->execute();

    }
    //**************************************************************************
    //***********************  T1 amplitude equations  *************************
    //**************************************************************************
//    OUT() << "\tSolving T1 Amplitude Equations" << std::endl;
    COMPILE(
      // Contract Kac and Kki with T1 amplitudes
      (*Rph)["ai"] <<= ( 1.0) * (*Kac)["ac"] * (*Tph)["ci"],
      (*Rph)["ai"] +=  (-1.0) * (*Kki)["ki"] * (*Tph)["ak"],

      //Build Kck
      (*Kck)["ck"] <<= ( 2.0) * (*Vhhpp)["klcd"] * (*Tph)["dl"],
      (*Kck)["ck"]  += (-1.0) * (*Vhhpp)["kldc"] * (*Tph)["dl"],

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rph)["ai"] += ( 2.0) * (*Kck)["ck"] * (*Tpphh)["caki"],
      (*Rph)["ai"] += (-1.0) * (*Kck)["ck"] * (*Tpphh)["caik"],
      (*Rph)["ai"] += ( 1.0) * (*Kck)["ck"] * (*Tph)["ci"] * (*Tph)["ak"],
      (*Rph)["ai"] += ( 2.0) * (*Vphhp)["akic"] * (*Tph)["ck"],
      (*Rph)["ai"] += (-1.0) * (*Vphph)["akci"] * (*Tph)["ck"],
      (*Rph)["ai"] +=
         (2.0) * (*cTGammaGpp)["Gac"] * (*GammaGhp)["Gkd"] * (*Xabij)["cdik"],
      (*Rph)["ai"] +=
        (-1.0) * (*cTGammaGpp)["Gad"] * (*GammaGhp)["Gkc"] * (*Xabij)["cdik"],

      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Xabij)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhhp)["lkic"] * (*Xabij)["ackl"]
    )->execute();

  }
  return residuum;
}

// instantiate
template class cc4s::Ccsd<Real<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsd<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsd<Real<64>, DefaultTensorEngine>;
template class cc4s::Ccsd<Complex<64>, DefaultTensorEngine>;

