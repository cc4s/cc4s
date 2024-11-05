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


  auto eigenEnergies(
    arguments->template getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  ASSERT_LOCATION(
    eigenEnergies, "Missing EigenEnergies!",
    arguments->sourceLocation
  );

  auto r3mode( arguments->template getValue<std::string>("r3-mode", "ct-full"));
  auto t3t3mode( arguments->template getValue<std::string>("t3t3-mode", "full"));


  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (*Vpphh)["abij"],
    (*Rppphhh)["abcijk"] <<= (*Vpphh)["abij"] * (*Rph)["ck"]
  )->execute();

  // TODO: check if given amplitudes contain expected parts
  // get amplitude parts
  auto Tph( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );


  auto Tppphhh( Tcc<TE>::template tensor<F>("Tppphhh") );

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
  auto Fhp( Tcc<TE>::template tensor<F>("Fhp") );


  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));
  auto Fepsh(Tcc<TE>::template tensor<F>(epsh->inspect()->getLens(), "Fepsh"));
  auto Fepsp(Tcc<TE>::template tensor<F>(epsp->inspect()->getLens(), "Fepsp"));
  auto fromReal( [](Real<> eps) {return F(eps);} );
  COMPILE(
    (*Fepsp)["a"] <<= map<F>(fromReal, (*epsp)["a"]),
    (*Fepsh)["i"] <<= map<F>(fromReal, (*epsh)["i"])
  )->execute();

  //intermediates
  COMPILE(
    (*Wphpp)["akcd"] <<= ( 2.0) * (*Vphpp)["akcd"],
    (*Wphpp)["akcd"]  += (-1.0) * (*Vphpp)["akdc"],
    (*Whhhp)["klic"] <<= ( 2.0) * (*Vhhhp)["klic"],
    (*Whhhp)["klic"]  += (-1.0) * (*Vhhhp)["lkic"],
    (*Whhpp)["klcd"] <<= ( 2.0) * (*Vhhpp)["klcd"],
    (*Whhpp)["klcd"]  += (-1.0) * (*Vhhpp)["lkcd"]
  )->execute();


  // T3 equations are taken from Noga & Bartlett JCP 86, 7041 (1987)
  // with consideration of the erratum

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Tph)["ai"],
    (*Xpphh)["abij"]  <<= (*Tpphh)["abij"],
    (*Xpphh)["abij"]   += (*Tph)["ai"] * (*Tph)["bj"],
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
  // We first construct the inital R3 object

  // which means that they will be 6th order or higher in the end
  COMPILE(
    (*Xabie)["abie"]  <<=          (*Vpphp)["abie"],
    (*Xamij)["amij"]  <<=          (*Vphhh)["amij"],
    (r3mode == "ct-ring") ? (
      (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*Tpphh)["afim"],
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*Tpphh)["aein"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    )
  )->execute();

  if (r3mode == "ct-full") {
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
    /*
    // these guys are sixth order when they couple back to R3
    (*Xabie)["abie"]   += (-2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfimn"],
    (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfnmi"],
    (*Xabie)["abie"]   += ( 1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["abfinm"],
    (*Xamij)["amij"]   += ( 2.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["aefijn"],
    (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["feaijn"],
    (*Xamij)["amij"]   += (-1.0) * (*Vhhpp)["mnef"] * (*Tppphhh)["afeijn"],
    */



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
    (*Rpphh)["abij"] <<= (*Xpphh)["abij"]

    /*
    // T1+T3 -> R2
    (*Xpphh)["abij"] <<=         (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"],
    (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"],
    (*Xpphh)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"],
    (*Xpphh)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"],
    (*Xpphh)["abij"] += (*Xpphh)["baji"],
    (*Rpphh)["abij"] += (*Xpphh)["abij"]
    */
  )->execute();

  // THIS IS NOT OUR STORY NOW - neglect for the time beeing
  //   Later we can contract this with the Tph to get the (T)-[T] contribution
  // T3 -> R1
  COMPILE(
    (*Rph)["ai"] += (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"] += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"],
    (*Rph)["ai"] += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"],
    (*Rph)["ai"] += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"]
  )->execute();



  auto denergy( Tcc<TE>::template tensor<F>("E") );
  auto senergy( Tcc<TE>::template tensor<F>("E") );
  COMPILE(
    (*denergy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Tpphh)["abij"],
    (*denergy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Tpphh)["baij"]
  )->execute();

  F Ed(denergy->read());

  COMPILE(
    (*senergy)[""] <<= 2.0 * (*Rph)["ai"] * (*Tph)["ai"]
  )->execute();
  F Es(senergy->read());

  if (r3mode == "ct-full") {
    OUT() << "(cT) / [cT] energy: " << real(Ed+Es) << " / " << real(Ed) << std::endl;
  } else if (r3mode == "ct-ring") {
    OUT() << "(cT)-ring / [cT]-ring energy: " << real(Ed+Es) << " / " << real(Ed) << std::endl;
  } else {
    OUT() << "(T) / [T] energy: " << real(Ed+Es) << " / " << real(Ed) << std::endl;
  }



  auto energy( Tcc<TE>::template tensor<F>("E") );
  // These are the T3-T3 coupling terms which can be of fifth order
  // we have removed all terms which introduce a higher order in perturbation
  // TODO: currently they not included as Tppphhh is zero at this point
  COMPILE(
    (*Xp3h3)["abcijk"] <<= (0.0) * (*Tppphhh)["becmjk"],
    (t3t3mode == "full" || t3t3mode == "ladder") ? (
      (*Xp3h3)["abcijk"]   <<=          (*Vhhhh)["jkmn"] * (*Tppphhh)["abcimn"],
      (*Xp3h3)["abcijk"]    +=          (*Vpppp)["bcef"] * (*Tppphhh)["aefijk"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    ),
    (t3t3mode == "full" || t3t3mode == "ring") ? (
      (*Xp3h3)["abcijk"]    += ( 2.0) * (*Vphhp)["amie"] * (*Tppphhh)["ebcmjk"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    ),
    (t3t3mode == "full" || t3t3mode == "rest") ? (
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Vphhp)["amie"] * (*Tppphhh)["becmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Vphhp)["amie"] * (*Tppphhh)["cbemjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Vphph)["amei"] * (*Tppphhh)["ebcmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Vphph)["bmei"] * (*Tppphhh)["aecmjk"],
      (*Xp3h3)["abcijk"]    += (-1.0) * (*Vphph)["cmei"] * (*Tppphhh)["abemjk"],
      Tcc<TE>::sequence()
    ) : (
      Tcc<TE>::sequence()
    )
  )->execute();
  COMPILE(
    (*Rppphhh)["abcijk"] <<= (*Xp3h3)["abcijk"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["bacjik"],
    (*Rppphhh)["abcijk"]  += (*Xp3h3)["cbakji"]
  )->execute();
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
    (*energy)[""] <<=  2.0 * (*Rpphh)["abij"] * (*Tpphh)["abij"],
    (*energy)[""]  += -1.0 * (*Rpphh)["abij"] * (*Tpphh)["baij"]
  )->execute();
  F E(energy->read());
  std::string r3string("(T)");
  if (r3mode == "ct-full") {
    r3string = "(cT)";
  } else if (r3mode == "ct-ring") {
    r3string = "(cT)-ring";
  }

  OUT() << "T3-T3 energy ( " << t3t3mode << " / " << r3string << " ) : " << real(E) << std::endl;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setValue("allGood", 1);
  return result;


}

