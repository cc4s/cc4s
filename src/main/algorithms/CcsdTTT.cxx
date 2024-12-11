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

#include <algorithms/CcsdTTT.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>
#include <TensorSet.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdTTT)


Ptr<MapNode> CcsdTTT::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  } else {
    using TE = DefaultTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}


template <typename F, typename TE>
Ptr<MapNode> CcsdTTT::run(const Ptr<MapNode> &arguments) {

  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto Rppphhh ( Tcc<TE>::template tensor<F>("Rppphhh") );

  auto amplitudes(
    arguments->template getPtr<TensorSet<F,TE>>("amplitudes")
  );

  ASSERT_LOCATION(
    amplitudes, "Missing amplitudes!",
    arguments->sourceLocation
  );

  auto coulombIntegrals(
    arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );
  ASSERT_LOCATION(
    coulombIntegrals, "Missing amplitudes!",
    arguments->sourceLocation
  );
  auto Vpphh(coulombIntegrals->get("pphh"));

  const bool lambda = (arguments->isGiven("lambdaAmplitudes")) ? true : false;

  Ptr<TensorSet<F,TE>> lAmplitudes
    = (lambda) ? arguments->getPtr<TensorSet<F,TE>>("lambdaAmplitudes") : nullptr;

  ASSERT_LOCATION(
    !lambda || lAmplitudes, "Lambda Amplitudes required but not available",
    arguments->sourceLocation
  );
  Ptr<TensorExpression<F,TE>> Zpphh = (lambda) ? lAmplitudes->get("pphh") : nullptr;
  Ptr<TensorExpression<F,TE>> Zph   = (lambda) ? lAmplitudes->get("ph") : nullptr;



  auto eigenEnergies(
    arguments->template getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  ASSERT_LOCATION(
    eigenEnergies, "Missing EigenEnergies!",
    arguments->sourceLocation
  );


  //t3-mode: plain, ring, ct-full
  //t3t3-mode: plain, dressed
  auto t3mode( arguments->template getValue<std::string>("t3mode", "dressed"));
  auto t3t3mode( arguments->template getValue<std::string>("t3t3mode", "dressed"));
  bool singles = arguments->template getValue<int>("singles", 1) == 1;
  OUT() << "T3-mode: " << t3mode << std::endl;
  OUT() << "T3T3-mode: " << t3t3mode << std::endl;
  OUT() << "Singles: " << singles << std::endl;

  //std::vector<std::string> contributions({"[VT2]","[VT3]" ,"[VT1T2]","[VT2T3]","[VT1T3]","[VT2T1T1]","[VT2T2T1]","[VT2T2T1]", "[VT3T1T1]", "[VT2T1T1T1]";

  ASSERT_LOCATION(
    lambda, "We need Lambda amplitudes!",
    arguments->sourceLocation
  );

  // TODO: check if given amplitudes contain expected parts
  // get amplitude parts
  auto _Tph( amplitudes->get("ph") );
  auto _Tpphh( amplitudes->get("pphh") );
  auto Tph( Tcc<TE>::template tensor<F>("Tph") );
  auto Tpphh( Tcc<TE>::template tensor<F>("Tph") );
  auto Tppphhh( Tcc<TE>::template tensor<F>("Tppphhh") );

  F factor = (singles) ? F(1) : F(0);
  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (0.0) * (*Vpphh)["abij"],
    (*Rppphhh)["abcijk"] <<= (0.0) * (*Vpphh)["abij"] * (*Rph)["ck"],
    (*Tppphhh)["abcijk"] <<= (0.0) * (*Vpphh)["abij"] * (*Rph)["ck"],
    (*Tph)["ai"] <<= factor * (*_Tph)["ai"],
    (*Tpphh)["abij"] <<= (*_Tpphh)["abij"]
  )->execute();


  Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");
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
  auto Wphpp( Tcc<TE>::template tensor<F>("Wphpp") );
  auto Whhhp( Tcc<TE>::template tensor<F>("Whhhp") );
  auto Whhpp( Tcc<TE>::template tensor<F>("Whhpp") );
  auto Xpphh( Tcc<TE>::template tensor<F>("Xpphh") );
  auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
  auto Xjkmn( Tcc<TE>::template tensor<F>("Xjkmn") );
  auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );
  auto Xim  ( Tcc<TE>::template tensor<F>("Xim")   );
  auto Xae  ( Tcc<TE>::template tensor<F>("Xae")   );
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
  auto Dpphh( Tcc<TE>::template tensor<F>("Dpphh") );
  auto Dph( Tcc<TE>::template tensor<F>("Dph") );
  auto Fhp( Tcc<TE>::template tensor<F>("Fhp") );
  auto Rph1(  Tcc<TE>::template tensor<F>("Rph1") );
  auto Ippphhh(  Tcc<TE>::template tensor<F>("Ippphhh") );
  auto denergy( Tcc<TE>::template tensor<F>("E") );
  auto senergy( Tcc<TE>::template tensor<F>("E") );


  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));
  auto Fepsh(Tcc<TE>::template tensor<F>(epsh->inspect()->getLens(), "Fepsh"));
  auto Fepsp(Tcc<TE>::template tensor<F>(epsp->inspect()->getLens(), "Fepsp"));
  auto fromReal( [](Real<> eps) {return F(eps);} );
  COMPILE(
    (*Fepsp)["a"] <<= map<F>(fromReal, (*epsp)["a"]),
    (*Fepsh)["i"] <<= map<F>(fromReal, (*epsh)["i"])
  )->execute();

  // T3 equations are taken from Noga & Bartlett JCP 86, 7041 (1987)
  // with consideration of the erratum

  //define intermediates
  COMPILE(
    (*Wphpp)["akcd"] <<= ( 2.0) * (*Vphpp)["akcd"],
    (*Wphpp)["akcd"]  += (-1.0) * (*Vphpp)["akdc"],
    (*Whhhp)["klic"] <<= ( 2.0) * (*Vhhhp)["klic"],
    (*Whhhp)["klic"]  += (-1.0) * (*Vhhhp)["lkic"],
    (*Whhpp)["klcd"] <<= ( 2.0) * (*Vhhpp)["klcd"],
    (*Whhpp)["klcd"]  += (-1.0) * (*Vhhpp)["lkcd"],
    (*Xpphh)["abij"] <<= (*Tpphh)["abij"],
    (*Xpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"]
  )->execute();

  std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates;

  intermediates["Wphpp"] = Wphpp;
  intermediates["Whhhp"] = Whhhp;
  intermediates["Whhpp"] = Whhpp;
  intermediates["Vhhpp"] = Vhhpp;
  intermediates["Zpphh"] = Zpphh;
  intermediates["Zph"]   = Zph;
  intermediates["Tpphh"] = Tpphh;
  intermediates["Tph"]   = Tph;
  intermediates["Fepsp"] = Fepsp;
  intermediates["Fepsh"] = Fepsh;
  intermediates["Vphpp"] = Vphpp;
  intermediates["Vhhhp"] = Vhhhp;
  intermediates["Xpphh"] = Xpphh;
  intermediates["Vhphh"] = Vhphh;



  // construct intermediates
  COMPILE(
    (*Fphpp)["amef"]  <<= (*Vphpp)["amef"],
    (*Fphpp)["amef"]   += (-1.0) * (*Vhhpp)["nmef"] * (*Tph)["an"],
    (*Fphhh)["eimn"]  <<= (*Vphhh)["eimn"],
    (*Fphhh)["eimn"]   += (*Vhhpp)["mnef"] * (*Tph)["fi"],
    (*Fphhp)["amie"]  <<= (*Vphhp)["amie"],
    (*Fphhp)["amie"]   += (*Vphpp)["amfe"] * (*Tph)["fi"],
    (*Fphph)["amei"]  <<= (*Vphph)["amei"],
    (*Fphph)["amei"]   += (*Vphpp)["amef"] * (*Tph)["fi"],
    (*Fhpph)["ieam"]  <<= (*Vhpph)["ieam"],
    (*Fhpph)["ieam"]   += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["an"],
    (*Fhphp)["iema"]  <<= (*Vhphp)["iema"],
    (*Fhphp)["iema"]   += (-1.0) * (*Vhhhp)["inme"] * (*Tph)["an"],
    (*Fhp)["me"] <<= (*Whhpp)["mnef"] * (*Tph)["fn"]
  )->execute();

  //T1+T2+T3->R3
  COMPILE(
    (*Xabie)["abie"]  <<=          (*Vpphp)["abie"],
    (*Xamij)["amij"]  <<=          (*Vphhh)["amij"],
    (t3mode == "ring") ? (
      (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*Tpphh)["afim"],
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*Tpphh)["aein"]
    ) : (
      Tcc<TE>::sequence()
    )
  )->execute();

  if (t3mode == "dressed") {
    COMPILE(
      (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*Xpphh)["abnm"],
      (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*Tpphh)["afim"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*Tpphh)["afim"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*Tpphh)["afmi"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*Tpphh)["bfmi"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphhp)["amie"] * (*Tph)["bm"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*Tph)["am"],
      (*Xabie)["abie"]   += ( 1.0) * (*Vpppp)["abfe"] * (*Tph)["fi"],

      (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*Xpphh)["efij"],
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*Tpphh)["aein"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*Tpphh)["aein"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*Tpphh)["eain"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*Tpphh)["eajn"],
      (*Xamij)["amij"]   +=          (*Fhpph)["ieam"] * (*Tph)["ej"],
      (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*Tph)["ei"],
      (*Xamij)["amij"]   += (-1.0) * (*Vhhhh)["ijnm"] * (*Tph)["an"],
      (*Xamij)["amij"]   +=          (*Fhp)["me"] * (*Tpphh)["aeij"]
    )->execute();
  }

  COMPILE(
    // Construct the R3 from T2 only
    (*Xp3h3)["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["acbikj"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cabkij"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bcajki"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"]
  )->execute();

  /* (cT) (T-2) energy contributions */
  getEnergy(intermediates, Rppphhh, Xp3h3, "T3-abdghj");
  /* get the first order t3 amplitudes */
  getAmplitudes(intermediates, Rppphhh, Tppphhh);

  /*
  COMPILE(
    (*Xp3h3)["abcijk"] <<= (*Fepsp)["a"],
    (*Xp3h3)["abcijk"]  += (*Fepsp)["b"],
    (*Xp3h3)["abcijk"]  += (*Fepsp)["c"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["i"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["j"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["k"]
  )->execute();
  Xp3h3->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Xp3h3)["abcijk"]
  )->execute();

  // T3 -> R2
  COMPILE(
    // Pure T3->R2
    (*Xpphh)["abij"] <<=          (*Wphpp)["akcd"] * (*Tppphhh)["cbdijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*Vphpp)["akcd"] * (*Tppphhh)["cdbijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*Whhhp)["klic"] * (*Tppphhh)["abckjl"],
    (*Xpphh)["abij"]  +=          (*Vhhhp)["klic"] * (*Tppphhh)["acbkjl"],
    (*Xpphh)["abij"]  += (*Xpphh)["baji"],
    (*Rpphh)["abij"] <<= (*Xpphh)["abij"],

    (mode == "T-n") ? (
      // T1+T3 -> R2
      (*Xpphh)["abij"] <<=         (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] += (*Xpphh)["baji"],
      (*Rpphh)["abij"] += (*Xpphh)["abij"]
    ) : (
      Tcc<TE>::sequence()
    )
    //Tcc<TE>::sequence()
  )->execute();

  // T3 -> R1
  COMPILE(
    (*Rph)["ai"] += (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"] += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"],
    (*Rph)["ai"] += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"] += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"]
  )->execute();



  COMPILE(
    (*denergy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Tpphh)["abij"],
    (*denergy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Tpphh)["baij"],
    (*senergy)[""] <<= 2.0 * (*Rph)["ai"] * (*Tph)["ai"]
  )->execute();
  F Ed(denergy->read());
  F Es(senergy->read());

  std::streamsize ss(std::cout.precision());
  std::cout << setprecision(10);
  if ( mode != "T-n") {
    std::string oline;
    if (r3mode == "ct-full")  oline = "(cT) / [cT] energy: ";
    else if (r3mode == "ct-ring") oline = "(cT)-ring / [cT]-ring energy: ";
    else oline = "(T) / [T] energy: ";
    OUT() << oline << real(Ed+Es) << " / " << real(Ed) << std::endl;
  }

  if (lambda) {
    COMPILE(
      (*denergy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Zpphh)["abij"],
      (*denergy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Zpphh)["baij"],
      (*senergy)[""] <<= 2.0 * (*Rph)["ai"] * (*Zph)["ai"]
    )->execute();
    Ed = denergy->read();
    Es = senergy->read();

    if ( mode == "T-n") {
      OUT() << "(T-2) correlation energy: " << real(Ed+Es) << std::endl;
    } else {
      std::string oline;
      if (r3mode == "ct-full") oline += "Λ(cT) / Λ[cT] energy: ";
      else if (r3mode == "ct-ring" ) oline += "Λ(cT)-ring / Λ[cT]-ring energy: ";
      else oline += "Λ(T) / Λ[T] energy: ";
      OUT() << oline << real(Ed+Es) << " / " << real(Ed) << std::endl;
    }
  }
  ENERGY BLOCK
  */
  /*
   *
  COMPILE(
    (*Dph)["ai"] <<= (0.0) * (*Tph)["ai"],
    (*Dph)["ai"]  += (*Fepsp)["a"],
    (*Dph)["ai"]  -= (*Fepsh)["i"],

    (*Dpphh)["abij"] <<= (0.0) * (*Tpphh)["abij"],
    (*Dpphh)["abij"]  += (*Fepsp)["a"],
    (*Dpphh)["abij"]  += (*Fepsp)["b"],
    (*Dpphh)["abij"]  -= (*Fepsh)["i"],
    (*Dpphh)["abij"]  -= (*Fepsh)["j"]
  )->execute();
  Dpphh->sum(F(1.0), "abij", "abij", [] (F delta) { return F(-1) / delta; });
  Dph->sum(F(1.0),   "ai",   "ai",   [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*denergy)[""] <<= ( 2.0) * (*Rpphh)["abij"] * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*denergy)[""]  += (-1.0) * (*Rpphh)["abij"] * (*Vpphh)["baij"] * (*Dpphh)["abij"],
    (*senergy)[""] <<= ( 2.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["abij"] * (*Dph)["ai"] * (*Dph)["bj"],
    (*senergy)[""]  += (-1.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["baij"] * (*Dph)["ai"] * (*Dph)["bj"],
    (*senergy)[""]  += ( 2.0) * (*Rph)["ai"] * (*Tph)["bj"] * (*Vpphh)["abij"] * (*Dph)["ai"],
    (*senergy)[""]  += (-1.0) * (*Rph)["ai"] * (*Tph)["bj"] * (*Vpphh)["baij"] * (*Dph)["ai"],
    (*senergy)[""]  += ( 2.0) * (*Tph)["ai"] * (*Rph)["bj"] * (*Vpphh)["abij"] * (*Dph)["bj"],
    (*senergy)[""]  += (-1.0) * (*Tph)["ai"] * (*Rph)["bj"] * (*Vpphh)["baij"] * (*Dph)["bj"],
    (*Rph1)["ab"] <<= (*Rph)["ab"]
  )->execute();
  Ed = denergy->read();
  Es = senergy->read();
  {
  if (r3mode == "ct-full") oline += "2(cT) / 2[cT] energy: ";
  else if (r3mode == "ct-ring") oline += "2(cT)-ring / 2[cT]-ring energy: ";
  else oline += "2(T) / 2[T] energy: ";
  OUT() << oline << real(Ed+Es) << " / " << real(Ed) << std::endl;
  }
  */

  COMPILE(
    (*Xjkmn)["jkmn"]  <<=          (*Vhhhh)["jkmn"],
    (*Xbcef)["bcef"]  <<=          (*Vpppp)["bcef"],
    (*Xamei)["amei"]  <<=          (*Vphph)["amei"],
    (*Xamie)["amie"]  <<=          (*Vphhp)["amie"],

    (t3t3mode == "plain") ? (
      Tcc<TE>::sequence()
    ) : (

      (*Xjkmn)["jkmn"]   +=          (*Vhhpp)["mnef"] * (*Xpphh)["efjk"],
      (*Xjkmn)["jkmn"]   +=          (*Vhphh)["jemn"] * (*Tph)["ek"],
      (*Xjkmn)["jkmn"]   +=          (*Vphhh)["ekmn"] * (*Tph)["ej"],

      (*Xbcef)["bcef"]   +=          (*Vpphh)["efmn"] * (*Xpphh)["bcmn"],
      (*Xbcef)["bcef"]   += (-1.0) * (*Vphpp)["bmef"] * (*Tph)["cm"],
      (*Xbcef)["bcef"]   += (-1.0) * (*Vhppp)["mcef"] * (*Tph)["bm"],

      (*Xamei)["amei"]   += (-1.0) * (*Vhhpp)["mnfe"] * (*Xpphh)["fain"],
      (*Xamei)["amei"]   += (-1.0) * (*Vhhph)["nmei"] * (*Tph)["an"],
      (*Xamei)["amei"]   +=          (*Vphpp)["amef"] * (*Tph)["fi"],

      (*Xamie)["amie"]   +=          (*Whhpp)["mnef"] * (*Tpphh)["afin"],
      (*Xamie)["amie"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Xpphh)["fain"],
      (*Xamie)["amie"]   += (-1.0) * (*Vhhhp)["nmie"] * (*Tph)["an"],
      (*Xamie)["amie"]   +=          (*Vppph)["aefm"] * (*Tph)["fi"]
    )
  )->execute();

  COMPILE(
    (*Xp3h3)["abcijk"] <<= (0.0) * (*Tppphhh)["abcijk"],
    (*Xp3h3)["abcijk"]   <<=          (*Xjkmn)["jkmn"] * (*Tppphhh)["abcimn"],
    (*Xp3h3)["abcijk"]    +=          (*Xbcef)["bcef"] * (*Tppphhh)["aefijk"],
    (*Xp3h3)["abcijk"]    += ( 2.0) * (*Xamie)["amie"] * (*Tppphhh)["ebcmjk"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["becmjk"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["cbemjk"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["amei"] * (*Tppphhh)["ebcmjk"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["bmei"] * (*Tppphhh)["aecmjk"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamei)["cmei"] * (*Tppphhh)["abemjk"],
    (*Rppphhh)["abcijk"] <<= (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"]
  )->execute();

  /* T3T3 coupling or (T-3) contributions */
  getEnergy(intermediates, Rppphhh, Xp3h3, "T3-cfi");

  r3_VT3T2_b(intermediates, Tph, Tpphh, Tppphhh, Xp3h3, Rppphhh);



  /*
  COMPILE(
    //  T3 <-- [VT3T2]
    (*Xabie)["abie"]  <<= (-2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfimn"],
    (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfnmi"],
    (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfinm"],
    (*Xamij)["amij"]  <<= ( 2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["aefijn"],
    (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["feaijn"],
    (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["afeijn"],

    (*Xp3h3)["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
    (*Rppphhh)["abcijk"] <<= (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["acbikj"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cabkij"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bcajki"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"],

    (*Xim)["im"]        <<= ( 2.0) * (*Vhphh)["iemn"] * (*Tph)["en"],
    (*Xim)["im"]         += (-1.0) * (*Vhphh)["ienm"] * (*Tph)["en"],
    (*Xim)["im"]         +=          (*Whhpp)["mnef"] * (*Xpphh)["efin"],
    (*Xae)["ae"]        <<= ( 2.0) * (*Vphpp)["amef"] * (*Tph)["fm"],
    (*Xae)["ae"]         += (-1.0) * (*Vphpp)["amfe"] * (*Tph)["fm"],
    (*Xae)["ae"]         += (-1.0) * (*Whhpp)["mnef"] * (*Xpphh)["afmn"],
    (*Xp3h3)["abcijk"]  <<=          (*Xae)["ae"]     * (*Tppphhh)["ebcijk"],
    (*Xp3h3)["abcijk"]   += (-1.0) * (*Xim)["im"]     * (*Tppphhh)["abcmjk"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["bacjik"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["cbakji"]

  )->execute();
  */
  /* T3T3 coupling or (T-3) contributions */
  getEnergy(intermediates, Rppphhh, Xp3h3, "T3-e");


  /*
   * replace
  COMPILE(
    (*Xp3h3)["abcijk"] <<= (*Fepsp)["a"],
    (*Xp3h3)["abcijk"]  += (*Fepsp)["b"],
    (*Xp3h3)["abcijk"]  += (*Fepsp)["c"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["i"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["j"],
    (*Xp3h3)["abcijk"]  -= (*Fepsh)["k"]
  )->execute();
  Xp3h3->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Xp3h3)["abcijk"]
  )->execute();

  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Xp3h3)["abcijk"],
    (*Xpphh)["abij"] <<=          (*Wphpp)["akcd"] * (*Tppphhh)["cbdijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*Vphpp)["akcd"] * (*Tppphhh)["cdbijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*Whhhp)["klic"] * (*Tppphhh)["abckjl"],
    (*Xpphh)["abij"]  +=          (*Vhhhp)["klic"] * (*Tppphhh)["acbkjl"],
    (*Xpphh)["abij"]  += (*Xpphh)["baji"],
    (*Rpphh)["abij"] <<= (*Xpphh)["abij"],

    (mode == "T-n") ? (
      // T1+T3 -> R2
      (*Xpphh)["abij"] <<=         (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"],
      (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"],
      (*Xpphh)["abij"] += (*Xpphh)["baji"],
      (*Rpphh)["abij"] += (*Xpphh)["abij"]
    ) : (
      Tcc<TE>::sequence()
    ),


    (*Rph)["ai"] <<= (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"]  += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"],
    (*Rph)["ai"]  += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"]  += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"],

    (*denergy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Tpphh)["abij"],
    (*denergy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Tpphh)["baij"],
    (*senergy)[""] <<= 2.0 * (*Rph)["ai"] * (*Tph)["ai"]
  )->execute();
  Ed = denergy->read();
  Es = senergy->read();
  std::string r3string("(T)");
  if (r3mode == "ct-full") {
    r3string = "(cT)";
  } else if (r3mode == "ct-ring") {
    r3string = "(cT)-ring";
  }

  if (mode != "T-n") {
    OUT() << "T3-T3 energy ( " << t3t3part << " / " << r3string << " / " << t3t3mode << " ) : " << real(Ed+Es) << " | " << real(Ed) << std::endl;
  }

  if (lambda) {
    COMPILE(
      (*denergy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Zpphh)["abij"],
      (*denergy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Zpphh)["baij"],
      (*senergy)[""] <<= 2.0 *  (*Rph)["ai"] * (*Zph)["ai"]
    )->execute();
    Ed = denergy->read();
    Es = senergy->read();
    std::string r3string("(T)");
    if (r3mode == "ct-full") {
      r3string = "(cT)";
    } else if (r3mode == "ct-ring") {
      r3string = "(cT)-ring";
    }
    if (mode == "T-n") {
      OUT() << "(T-3) correlation energy: " << real(Ed+Es) << std::endl;
    } else {
      OUT() << "Λ T3-T3 energy ( " << t3t3part << " / " << r3string << " / " << t3t3mode << " ) : " << real(Ed+Es) << " | " << real(Ed) << std::endl;
    }

  }
  std::cout << std::setprecision(ss);
  */
  /*
  COMPILE(
    (*denergy)[""] <<= ( 2.0) * (*Rpphh)["abij"] * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*denergy)[""]  += (-1.0) * (*Rpphh)["abij"] * (*Vpphh)["baij"] * (*Dpphh)["abij"],
    (*senergy)[""] <<= ( 2.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["abij"] * (*Dph)["ai"] * (*Dph)["bj"],
    (*senergy)[""]  += (-1.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["baij"] * (*Dph)["ai"] * (*Dph)["bj"],
    (*senergy)[""]  += ( 4.0) * (*Rph)["ai"] * (*Tph)["bj"] * (*Vpphh)["abij"] * (*Dph)["ai"],
    (*senergy)[""]  += (-2.0) * (*Rph)["ai"] * (*Tph)["bj"] * (*Vpphh)["baij"] * (*Dph)["ai"],
    (*senergy)[""]  += ( 4.0) * (*Rph)["ai"] * (*Rph1)["bj"] * (*Vpphh)["abij"] * (*Dph)["ai"] * (*Dph)["bj"],
    (*senergy)[""]  += (-2.0) * (*Rph)["ai"] * (*Rph1)["bj"] * (*Vpphh)["baij"] * (*Dph)["ai"] * (*Dph)["bj"]

  )->execute();
  Ed = denergy->read();
  Es = senergy->read();
  OUT() << "2 T3-T3 energy ( " << t3t3part << " / " << r3string << " / " << t3t3mode << " ) : tot/dou. " << real(Ed+Es) << " | " << real(Ed) << std::endl;
  */



  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setValue("allGood", 1);
  return result;


}


template <typename F, typename TE>
void CcsdTTT::getEnergy(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                        Ptr<Tensor<F,TE>> Rppphhh,
                        Ptr<Tensor<F,TE>> Xppphhh,
                        std::string name) {

  auto Xpphh( Tcc<TE>::template tensor<F>("Xpphh") );
  auto Xph( Tcc<TE>::template tensor<F>("Xph") );
  auto tenergy( Tcc<TE>::template tensor<F>("Et") );
  auto lenergy( Tcc<TE>::template tensor<F>("El") );

  COMPILE(
    (*Xppphhh)["abcijk"] <<= (*intermediates["Fepsp"])["a"],
    (*Xppphhh)["abcijk"]  += (*intermediates["Fepsp"])["b"],
    (*Xppphhh)["abcijk"]  += (*intermediates["Fepsp"])["c"],
    (*Xppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["i"],
    (*Xppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["j"],
    (*Xppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["k"]
  )->execute();
  Xppphhh->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*Xppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Xppphhh)["abcijk"]
  )->execute();


  // T3 -> R2
  COMPILE(
    // Pure T3->R2
    (*Xpphh)["abij"] <<=          (*intermediates["Wphpp"])["akcd"] * (*Xppphhh)["cbdijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*intermediates["Vphpp"])["akcd"] * (*Xppphhh)["cdbijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*intermediates["Whhhp"])["klic"] * (*Xppphhh)["abckjl"],
    (*Xpphh)["abij"]  +=          (*intermediates["Vhhhp"])["klic"] * (*Xppphhh)["acbkjl"],
    (*Xpphh)["abij"]  += (*Xpphh)["baji"],
    (*tenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*intermediates["Tpphh"])["abij"],
    (*tenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*intermediates["Tpphh"])["baij"],
    (*lenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*intermediates["Zpphh"])["abij"],
    (*lenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*intermediates["Zpphh"])["baij"]
  )->execute();

  OUT() << name << "-T2a (T2/L2): " << real(tenergy->read()) << " " << real(lenergy->read()) << std::endl;

  COMPILE(
    // T1+T3 -> R2
    (*Xpphh)["abij"] <<=         (*intermediates["Whhpp"])["klcd"] * (*Xppphhh)["abcijk"] * (*intermediates["Tph"])["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*intermediates["Whhpp"])["klcd"] * (*Xppphhh)["acbijk"] * (*intermediates["Tph"])["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*intermediates["Whhpp"])["klcd"] * (*Xppphhh)["acbikl"] * (*intermediates["Tph"])["dj"],
    (*Xpphh)["abij"] +=          (*intermediates["Vhhpp"])["klcd"] * (*Xppphhh)["cabikl"] * (*intermediates["Tph"])["dj"],
    (*Xpphh)["abij"] += (-1.0) * (*intermediates["Whhpp"])["klcd"] * (*Xppphhh)["adcijk"] * (*intermediates["Tph"])["bl"],
    (*Xpphh)["abij"] +=          (*intermediates["Vhhpp"])["klcd"] * (*Xppphhh)["cdaijk"] * (*intermediates["Tph"])["bl"],
    (*Xpphh)["abij"] += (*Xpphh)["baji"],
    (*tenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*intermediates["Tpphh"])["abij"],
    (*tenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*intermediates["Tpphh"])["baij"],
    (*lenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*intermediates["Zpphh"])["abij"],
    (*lenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*intermediates["Zpphh"])["baij"]
  )->execute();

  OUT() << name << "-T2b (T2/L2): " << real(tenergy->read()) << " " << real(lenergy->read()) << std::endl;
  // T3 -> R1
  COMPILE(
    (*Xph)["ai"] += (-2.0) * (*intermediates["Vhhpp"])["jkbc"] * (*Xppphhh)["bacjki"],
    (*Xph)["ai"] += (+2.0) * (*intermediates["Vhhpp"])["jkbc"] * (*Xppphhh)["bcajki"],
    (*Xph)["ai"] += (+1.0) * (*intermediates["Vhhpp"])["jkcb"] * (*Xppphhh)["bacjki"],
    (*Xph)["ai"] += (-1.0) * (*intermediates["Vhhpp"])["jkcb"] * (*Xppphhh)["bcajki"],
    (*tenergy)[""] <<=  2.0 * (*Xph)["ai"] * (*intermediates["Tph"])["ai"],
    (*lenergy)[""] <<=  2.0 * (*Xph)["ai"] * (*intermediates["Zph"])["ai"]
  )->execute();

  OUT() << name << "-T1 (T1/L1): " << real(tenergy->read()) << " " << real(lenergy->read()) << std::endl;


}


template <typename F, typename TE>
void CcsdTTT::getAmplitudes(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                            Ptr<Tensor<F,TE>> Rppphhh,
                            Ptr<Tensor<F,TE>> Tppphhh) {
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*intermediates["Fepsp"])["a"],
    (*Tppphhh)["abcijk"]  += (*intermediates["Fepsp"])["b"],
    (*Tppphhh)["abcijk"]  += (*intermediates["Fepsp"])["c"],
    (*Tppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["i"],
    (*Tppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["j"],
    (*Tppphhh)["abcijk"]  -= (*intermediates["Fepsh"])["k"]
  )->execute();
  Tppphhh->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Tppphhh)["abcijk"]
  )->execute();
}


template <typename F, typename TE>
void CcsdTTT::r3_VT3T2_b(std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates,
                Ptr<Tensor<F,TE>> Tph,
                Ptr<Tensor<F,TE>> Tpphh,
                Ptr<Tensor<F,TE>> Tppphhh,
                Ptr<Tensor<F,TE>> Xp3h3,
                Ptr<Tensor<F,TE>> Rppphhh) {

  auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
  auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );
  auto Xim( Tcc<TE>::template tensor<F>("Xim") );
  auto Xae( Tcc<TE>::template tensor<F>("Xae") );
  COMPILE(
    //  T3 <-- [VT3T2]
    (*Xabie)["abie"]  <<= (-2.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["abfimn"],
    (*Xabie)["abie"]   += ( 1.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["abfnmi"],
    (*Xabie)["abie"]   += ( 1.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["abfinm"],
    (*Xamij)["amij"]  <<= ( 2.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["aefijn"],
    (*Xamij)["amij"]   += (-1.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["feaijn"],
    (*Xamij)["amij"]   += (-1.0) * (*intermediates["Vhhpp"])["mnef"] * (*Tppphhh)["afeijn"],

    (*Xp3h3)["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*Xp3h3)["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
    (*Rppphhh)["abcijk"] <<= (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["acbikj"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cabkij"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bcajki"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"],

    (*Xim)["im"]        <<= ( 2.0) * (*intermediates["Vhphh"])["iemn"] * (*Tph)["en"],
    (*Xim)["im"]         += (-1.0) * (*intermediates["Vhphh"])["ienm"] * (*Tph)["en"],
    (*Xim)["im"]         +=          (*intermediates["Whhpp"])["mnef"] * (*intermediates["Xpphh"])["efin"],
    (*Xae)["ae"]        <<= ( 2.0) * (*intermediates["Vphpp"])["amef"] * (*Tph)["fm"],
    (*Xae)["ae"]         += (-1.0) * (*intermediates["Vphpp"])["amfe"] * (*Tph)["fm"],
    (*Xae)["ae"]         += (-1.0) * (*intermediates["Whhpp"])["mnef"] * (*intermediates["Xpphh"])["afmn"],
    (*Xp3h3)["abcijk"]  <<=          (*Xae)["ae"]     * (*Tppphhh)["ebcijk"],
    (*Xp3h3)["abcijk"]   += (-1.0) * (*Xim)["im"]     * (*Tppphhh)["abcmjk"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["bacjik"],
    (*Rppphhh)["abcijk"] += (*Xp3h3)["cbakji"]

  )->execute();

}
