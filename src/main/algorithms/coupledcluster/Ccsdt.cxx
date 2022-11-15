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

#include <algorithms/coupledcluster/Ccsdt.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>

using namespace cc4s;

template <typename F, typename TE>
CoupledClusterMethodRegistrar<
  F,TE,Ccsdt<F,TE>
> Ccsdt<F,TE>::registrar_("Ccsdt");


template <typename F, typename TE>
Ptr<TensorSet<F,TE>> Ccsdt<F,TE>::getResiduum(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {
  // construct residuum. Shape will be assumed upon first use.
  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto Rppphhh ( Tcc<TE>::template tensor<F>("Rppphhh") );
  auto residuum(
    New<TensorSet<F,TE>>(
      std::map<std::string,Ptr<TensorExpression<F,TE>>>(
        { {"ph",Rph}, {"pphh",Rpphh}, {"ppphhh", Rppphhh} }
      )
    )
  );

  auto coulombIntegrals(
    this->arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );
  auto Vpphh(coulombIntegrals->get("pphh"));


  if (!amplitudes) {
    // no previous amplitudes given
    COMPILE(
      (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= (*Vpphh)["abij"],
      (*Rppphhh)["abcijk"] <<= (*Vpphh)["abij"] * (*Rph)["ck"]
    )->execute();
  } else {
    // TODO: check if given amplitudes contain expected parts
    // get amplitude parts
    auto Tph( amplitudes->get("ph") );
    auto Tpphh( amplitudes->get("pphh") );
    auto Tppphhh( amplitudes->get("ppphhh") );

    if (!Tppphhh) {
      Tppphhh = Tcc<TE>::template tensor<F>("Tppphhh");
      COMPILE(
        (*Tppphhh)["abcijk"] <<= (0.0) * (*Vpphh)["abij"] * (*Tph)["ck"]
      )->execute();
      amplitudes->get("ppphhh") = Tppphhh;
    }

    Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");
    Tppphhh->inspect()->setName("Tppphhh");
    auto Vpppp(coulombIntegrals->get("pppp"));
    auto Vphph(coulombIntegrals->get("phph"));
    auto Vhhhh(coulombIntegrals->get("hhhh"));
    auto Vhhhp(coulombIntegrals->get("hhhp"));
    auto Vhphh(coulombIntegrals->get("hphh"));
    auto Vppph(coulombIntegrals->get("ppph"));
    auto Vhhpp(coulombIntegrals->get("hhpp"));
    auto Vpphp(coulombIntegrals->get("pphp"));
    auto Vphhh(coulombIntegrals->get("phhh"));
    auto Vhphp(coulombIntegrals->get("hphp"));
    auto Vhpph(coulombIntegrals->get("hpph"));
    auto Vphhp(coulombIntegrals->get("phhp"));
    auto Vphpp(coulombIntegrals->get("phpp"));
    auto Vhhph(coulombIntegrals->get("hhph"));
    auto Vhppp(coulombIntegrals->get("hppp"));
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
      (*Rpphh)["abij"] <<= ( 1.0) * (*Lac)["ac"] * (*Tpphh)["cbij"],

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
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tpphh)["cdik"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tpphh)["cdik"],
      (*Rph)["ai"] += ( 2.0) * (*Vphpp)["akcd"] * (*Tph)["ci"] * (*Tph)["dk"],
      (*Rph)["ai"] += (-1.0) * (*Vphpp)["akdc"] * (*Tph)["ci"] * (*Tph)["dk"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tpphh)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tpphh)["ackl"],
      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Tph)["ak"] * (*Tph)["cl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhph)["klci"] * (*Tph)["ak"] * (*Tph)["cl"]
    )->execute();

    // T3 equations are taken from Noga & Bartlett JCP 86, 7041 (1987)
    // with consideration of the erratum

    // T3 -> R1
    // from:
    // https://github.com/ValeevGroup/mpqc/blob/master/src/mpqc/chemistry/qc/lcao/cc/ccsdt.h#L446
    // bug3?
    COMPILE(
      (*Rph)["ai"] += (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"],
      (*Rph)["ai"] += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"],
      (*Rph)["ai"] += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"],
      (*Rph)["ai"] += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"]
    )->execute();

    // T3 -> R2
    auto Wphpp( Tcc<TE>::template tensor<F>("Wphpp") );
    auto Whhhp( Tcc<TE>::template tensor<F>("Whhhp") );
    auto Whhpp( Tcc<TE>::template tensor<F>("Whhpp") );
    auto Xpphh( Tcc<TE>::template tensor<F>("Xpphh") );


    // https://github.com/ValeevGroup/mpqc/blob/master/src/mpqc/chemistry/qc/lcao/cc/ccsdt.h#L620
    COMPILE(
      // Pure T3->R2
      (*Wphpp)["akcd"] <<= ( 2.0) * (*Vphpp)["akcd"],
      (*Wphpp)["akcd"]  += (-1.0) * (*Vphpp)["akdc"],
      (*Whhhp)["klic"] <<= ( 2.0) * (*Vhhhp)["klic"],
      (*Whhhp)["klic"]  += (-1.0) * (*Vhhhp)["lkic"],
      (*Whhpp)["klcd"] <<= ( 2.0) * (*Vhhpp)["klcd"],
      (*Whhpp)["klcd"]  += (-1.0) * (*Vhhpp)["lkcd"],

      (*Xpphh)["abij"] <<=          (*Wphpp)["akcd"] * (*Tppphhh)["cbdijk"],
      (*Xpphh)["abij"]  += (-1.0) * (*Vphpp)["akcd"] * (*Tppphhh)["cdbijk"],
      (*Xpphh)["abij"]  += (-1.0) * (*Whhhp)["klic"] * (*Tppphhh)["abckjl"],
      (*Xpphh)["abij"]  +=          (*Vhhhp)["klic"] * (*Tppphhh)["acbkjl"],
      (*Xpphh)["abij"]  += (*Xpphh)["baji"],
      (*Rpphh)["abij"]  += (*Xpphh)["abij"],
      // T1+T3 -> R2
      (*Xpphh)["abij"] <<=         (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] += (*Xpphh)["baji"],
      (*Rpphh)["abij"] += (*Xpphh)["abij"]
    )->execute();
    //T1+T2+T3->R3
    auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
    auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );
    auto Xim  ( Tcc<TE>::template tensor<F>("Xim")   );
    auto Xae  ( Tcc<TE>::template tensor<F>("Xae")   );
    auto Xjkmn( Tcc<TE>::template tensor<F>("Xjkmn") );
    auto Xbcef( Tcc<TE>::template tensor<F>("Xbcef") );
    auto Xamie( Tcc<TE>::template tensor<F>("Xamie") );
    auto Xamei( Tcc<TE>::template tensor<F>("Xamei") );
    auto Xp3h3( Tcc<TE>::template tensor<F>("Xp3h3") );
    auto Fphpp( Tcc<TE>::template tensor<F>("Fphpp") );
    auto Fphhh( Tcc<TE>::template tensor<F>("Fphhh") );
    auto Fphhp( Tcc<TE>::template tensor<F>("Fphhp") );
    auto Fphph( Tcc<TE>::template tensor<F>("Fphph") );
    auto Fhpph( Tcc<TE>::template tensor<F>("Fhpph") );
    auto Fhphp( Tcc<TE>::template tensor<F>("Fhphp") );
    auto Fhp( Tcc<TE>::template tensor<F>("Fhp") );

    // build epsa and epsi for F
    auto eigenEnergies(
    // this->arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
      this->arguments->template getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
    );
    auto epsh(eigenEnergies->get("h"));
    auto epsp(eigenEnergies->get("p"));
    auto Fepsh(Tcc<TE>::template tensor<F>(epsh->inspect()->getLens(), "Fepsh"));
    auto Fepsp(Tcc<TE>::template tensor<F>(epsp->inspect()->getLens(), "Fepsp"));
    // convert to type F (either complex or double)
    auto fromReal( [](Real<> eps) {return F(eps);} );
    COMPILE(
      (*Fepsp)["a"] <<= map<F>(fromReal, (*epsp)["a"]),
      (*Fepsh)["i"] <<= map<F>(fromReal, (*epsh)["i"])
    )->execute();

    ///////////////////////
    // NOGA PAPER EQUATIONS
    ///////////////////////
    COMPILE(

      // II.16
      (*Xpphh)["abij"]  <<= (*Tpphh)["abij"],
      (*Xpphh)["abij"]   += (*Tph)["ai"] * (*Tph)["bj"],

      // II.9
      (*Fphpp)["amef"]  <<= (*Vphpp)["amef"],
      (*Fphpp)["amef"]   += (-1.0) * (*Vhhpp)["nmef"] * (*Tph)["an"],

      // II.10
      (*Fphhh)["eimn"]  <<= (*Vphhh)["eimn"],
      (*Fphhh)["eimn"]   += (*Vhhpp)["mnef"] * (*Tph)["fi"],

      // II.11
      (*Fphhp)["amie"]  <<= (*Vphhp)["amie"],
      (*Fphhp)["amie"]   += (*Vphpp)["amfe"] * (*Tph)["fi"],

      // II.12
      (*Fphph)["amei"]  <<= (*Vphph)["amei"],
      (*Fphph)["amei"]   += (*Vphpp)["amef"] * (*Tph)["fi"],

      // II.13
      (*Fhpph)["ieam"]  <<= (*Vhpph)["ieam"],
      // ERRATUM
      (*Fhpph)["ieam"]   += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["an"],

      // II.14
      (*Fhphp)["iema"]  <<= (*Vhphp)["iema"],
      // ERRATUM
      (*Fhphp)["iema"]   += (-1.0) * (*Vhhhp)["inme"] * (*Tph)["an"],

      // II.15
      // bug 3
      (*Fhp)["me"] <<= (*Whhpp)["mnef"] * (*Tph)["fn"],

      // II.1
      (*Xabie)["abie"]  <<=          (*Vpphp)["abie"],
      (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*Xpphh)["abnm"],
      (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*Tpphh)["afim"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*Tpphh)["afim"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*Tpphh)["afmi"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*Tpphh)["bfmi"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphhp)["amie"] * (*Tph)["bm"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*Tph)["am"],
      (*Xabie)["abie"]   += ( 1.0) * (*Vpppp)["abfe"] * (*Tph)["fi"],
      (*Xabie)["abie"]   += (-2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfimn"],
      (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfnmi"],
      (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfinm"],


      // II.2
      (*Xamij)["amij"]  <<=          (*Vphhh)["amij"],
      (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*Xpphh)["efij"],
      // ERRATUM
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*Tpphh)["aein"],
      // ERRATUM
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*Tpphh)["aein"],
      // ERRATUM
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*Tpphh)["eain"],
      // ERRATUM
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*Tpphh)["eajn"],
      (*Xamij)["amij"]   +=          (*Fhpph)["ieam"] * (*Tph)["ej"],
      (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*Tph)["ei"],
      (*Xamij)["amij"]   += (-1.0) * (*Vhhhh)["ijnm"] * (*Tph)["an"],
      (*Xamij)["amij"]   +=          (*Fhp)["me"] * (*Tpphh)["aeij"],
      (*Xamij)["amij"]   += ( 2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["aefijn"],
      (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["feaijn"],
      (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["afeijn"],

      // II.3
      // NOTE: we antisymmetrize in place Vhphh
      (*Xim)["im"]      <<= ( 2.0) * (*Vhphh)["iemn"] * (*Tph)["en"],
      // We dont need this term as we have 0 on the lhs
      //      (*Xim)["ii"]       += (*Fepsh)["i"],
      (*Xim)["im"]       += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["en"],
      (*Xim)["im"]       +=          (*Whhpp)["mnef"] * (*Xpphh)["efin"],

      // II.4
      // NOTE: we antisymmetrize in place Vphpp
      (*Xae)["ae"]      <<= ( 2.0) * (*Vphpp)["amef"] * (*Tph)["fm"],
      // We dont need this term as we have 0 on the lhs
      //      (*Xae)["aa"]       += (*Fepsp)["a"],
      (*Xae)["ae"]       += (-1.0) * (*Vphpp)["amfe"] * (*Tph)["fm"],
      (*Xae)["ae"]       += (-1.0) * (*Whhpp)["mnef"] * (*Xpphh)["afmn"],

      // II.5
      (*Xjkmn)["jkmn"]  <<=          (*Vhhhh)["jkmn"],
      // ERRATUM
      (*Xjkmn)["jkmn"]   +=          (*Vhhpp)["mnef"] * (*Xpphh)["efjk"],
      (*Xjkmn)["jkmn"]   +=          (*Vhphh)["jemn"] * (*Tph)["ek"],
      (*Xjkmn)["jkmn"]   +=          (*Vphhh)["ekmn"] * (*Tph)["ej"],

      // II.6
      (*Xbcef)["bcef"]  <<=          (*Vpppp)["bcef"],
      // ERRATUM
      (*Xbcef)["bcef"]   +=          (*Vpphh)["efmn"] * (*Xpphh)["bcmn"],
      (*Xbcef)["bcef"]   += (-1.0) * (*Vphpp)["bmef"] * (*Tph)["cm"],
      (*Xbcef)["bcef"]   += (-1.0) * (*Vhppp)["mcef"] * (*Tph)["bm"],

      // II.7
      (*Xamei)["amei"]  <<=          (*Vphph)["amei"],
      (*Xamei)["amei"]   += (-1.0) * (*Vhhpp)["mnfe"] * (*Xpphh)["fain"],
      (*Xamei)["amei"]   += (-1.0) * (*Vhhph)["nmei"] * (*Tph)["an"],
      (*Xamei)["amei"]   +=          (*Vphpp)["amef"] * (*Tph)["fi"],

      // II.8
      (*Xamie)["amie"]  <<=          (*Vphhp)["amie"],
      (*Xamie)["amie"]   +=          (*Whhpp)["mnef"] * (*Xpphh)["afin"],
      (*Xamie)["amie"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Xpphh)["fain"],
      (*Xamie)["amie"]   += (-1.0) * (*Vhhhp)["nmie"] * (*Tph)["an"],
      (*Xamie)["amie"]   +=          (*Vppph)["aefm"] * (*Tph)["fi"],

      // I.3 until I.10 (numbers are counting terms)
      (*Xp3h3)["abcijk"]   <<=          (*Xjkmn)["jkmn"] * (*Tppphhh)["abcimn"],
      (*Xp3h3)["abcijk"]    +=          (*Xbcef)["bcef"] * (*Tppphhh)["aefijk"],
      (*Xp3h3)["abcijk"]    += ( 2.0) * (*Xamie)["amie"] * (*Tppphhh)["ebcmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["becmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["cbemjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["amei"] * (*Tppphhh)["ebcmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["bmei"] * (*Tppphhh)["aecmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["cmei"] * (*Tppphhh)["abemjk"],
      (*Xp3h3)["abcijk"]    +=          (*Xae)["ae"]     * (*Tppphhh)["ebcijk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xim)["im"]     * (*Tppphhh)["abcmjk"],

      // Permute I.3 until I.10 through P(ia/jb, kc)
      (*Rppphhh)["abcijk"] <<= (*Xp3h3)["abcijk"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"],

      // I.1 and I.2
      (*Xp3h3)["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["abcijk"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["acbikj"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["cabkij"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["bcajki"],
      (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"]
    )->execute();


  }
  return residuum;
}

// instantiate
template class cc4s::Ccsdt<Real<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsdt<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsdt<Real<64>, DefaultTensorEngine>;
template class cc4s::Ccsdt<Complex<64>, DefaultTensorEngine>;

