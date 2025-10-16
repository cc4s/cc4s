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

#include <typeinfo>

#include <algorithms/CcsdtNonit.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>
#include <TensorSet.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdtNonit)


/*
 * Implementation Overview:
 *
 * This code implements the (T-n) perturbation series as described by
 * Erikson, Jorgensen, and Gauss in *JCP* **142**, 014102 (2015).
 * The primary objective of the (T-n) series is to achieve convergence
 * towards the CCSDT (Coupled Cluster Singles, Doubles, and Triples) limit.
 *
 * Key Concept:
 * Unlike traditional perturbation methods that correct the Hartree-Fock (HF)
 * reference energy, this approach is a perturbative correction to the CCSD
 * (Coupled Cluster Singles and Doubles) energy.
 *
 * Implementation Details:
 * This implementation offers two operational modes:
 *   1. **Multiplier-based Mode**: Utilizes Lagrange multipliers obtained
 *      from LambdaCCSD equations to refine the left-hand side of the
 *      perturbation terms (original theory as discribed by Erikson et al.)
 *   2. **Amplitude-based Mode**: Employs t-amplitudes directly from
 *      the left-hand side for corrections (CCSD(T)-like energy expression) .
 *   3. **Coupled-cluster Mode**: Employs mp2-amplitudes (ie bare Coulomb
 *      kernel) from left (as done in iterative CC methods)
 *
 */

Ptr<MapNode> CcsdtNonit::run(const Ptr<MapNode> &arguments) {
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
Ptr<MapNode> CcsdtNonit::run(const Ptr<MapNode> &arguments) {

  if (typeid(Complex<>) == typeid(F)) {
    ASSERT_LOCATION(false,
                    "CcsdtNonit does not work for complex numbers",
                    arguments->sourceLocation);
  }

  auto coulombIntegrals(
    arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );

  if (!coulombIntegrals) return nullptr;

  auto Vpphh(coulombIntegrals->get("pphh"));

  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto Rppphhh ( Tcc<TE>::template tensor<F>("Rppphhh") );

  Ptr<TensorSet<F,TE>> amplitudes;
  amplitudes = arguments->getPtr<TensorSet<F,TE>>("amplitudes");

  ASSERT_LOCATION(
    amplitudes, "Missing t amplitudes!",
    arguments->sourceLocation
  );
  auto _Tph( amplitudes->get("ph") );
  auto _Tpphh( amplitudes->get("pphh") );

  context["method"] = "CcsdtNonit";
  context["occ"] = _Tph->inspect()->lens[1];
  context["vir"] = _Tph->inspect()->lens[0];
  context["type"] = TypeTraits<F>::getName();



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

  auto t3t3mode( arguments->template getValue<std::string>("t3t3mode", "dressed"));
  bool singles = arguments->template getValue<int>("singles", 1) == 1;

  if (!singles) OUT() << "Warning: Singles turned off!!\n";
  if (! (t3t3mode == "dressed" || t3t3mode == "3CC") ) OUT() << "Warning: t3t3mode: " << t3t3mode << '\n';

  ASSERT_LOCATION(
    lambda, "We need Lambda amplitudes!",
    arguments->sourceLocation
  );

  // TODO: check if given amplitudes contain expected parts
  // get amplitude parts
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

  auto T2pphh( Tcc<TE>::template tensor<F>("T2pphh") );
  auto T2ph( Tcc<TE>::template tensor<F>("T2ph") );


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
    (*Xpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
    (*Xp3h3)["abcijk"] <<= (0.0) * (*Vpphh)["abij"] * (*Rph)["ck"]
  )->execute();

  std::map<std::string, Ptr<TensorExpression<F,TE>>> intermediates;

  intermediates["Xp3h3"] = Xp3h3;
  intermediates["Wphpp"] = Wphpp;
  intermediates["Whhhp"] = Whhhp;
  intermediates["Whhpp"] = Whhpp;
  intermediates["Vhhpp"] = Vhhpp;
  intermediates["Zpphh"] = Zpphh;
  intermediates["Zph"]   = Zph;
  intermediates["Tpphh"] = Tpphh;
  intermediates["Tph"]   = Tph;
  intermediates["Xpphh"] = Xpphh;
  intermediates["Fepsp"] = Fepsp;
  intermediates["Fepsh"] = Fepsh;
  intermediates["Vphpp"] = Vphpp;
  intermediates["Vhhhp"] = Vhhhp;
  intermediates["Vhphh"] = Vhphh;
  intermediates["Vhhhh"] = Vhhhh;
  intermediates["Vpppp"] = Vpppp;
  intermediates["Vphph"] = Vphph;
  intermediates["Vphhp"] = Vphhp;
  intermediates["Vphhh"] = Vphhh;
  intermediates["Vpphh"] = Vpphh;
  intermediates["Vhppp"] = Vhppp;
  intermediates["Vppph"] = Vppph;
  intermediates["Vhhph"] = Vhhph;
  intermediates["Vhpph"] = Vhpph;
  intermediates["Vhphp"] = Vhphp;
  intermediates["Vpphp"] = Vpphp;


  r3_VT2T1(intermediates, Tph, Tpphh, Rppphhh, t3t3mode);

  /* (cT) (T-2) energy contributions */
  auto Et2 = getEnergy(intermediates, Rppphhh, "T3-abdghj | (T-2)");
  /* get the first order t3 amplitudes */
  getAmplitudes(intermediates, Rppphhh, Tppphhh);

  //This is the 2nd order T1 and T2 amplitude for (T-4)(relax) II
  COMPILE(
    // Pure T3->R2
    (*T2pphh)["abij"] <<=          (*Wphpp)["akcd"] * (*Tppphhh)["cbdijk"],
    (*T2pphh)["abij"]  += (-1.0) * (*Vphpp)["akcd"] * (*Tppphhh)["cdbijk"],
    (*T2pphh)["abij"]  += (-1.0) * (*Whhhp)["klic"] * (*Tppphhh)["abckjl"],
    (*T2pphh)["abij"]  +=          (*Vhhhp)["klic"] * (*Tppphhh)["acbkjl"],
    (*T2pphh)["abij"]  +=          (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*Tph)["dl"],
    (*T2pphh)["abij"]  += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*Tph)["dl"],
    (*T2pphh)["abij"]  += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*Tph)["dj"],
    (*T2pphh)["abij"]  +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*Tph)["dj"],
    (*T2pphh)["abij"]  += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*Tph)["bl"],
    (*T2pphh)["abij"]  +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*Tph)["bl"],
    (*T2pphh)["abij"]  += (*T2pphh)["baji"],
    (*T2ph)["ai"]     <<= (-2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bacjki"],
    (*T2ph)["ai"]      += (+2.0) * (*Vhhpp)["jkbc"] * (*Tppphhh)["bcajki"],
    (*T2ph)["ai"]      += (+1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bacjki"],
    (*T2ph)["ai"]      += (-1.0) * (*Vhhpp)["jkcb"] * (*Tppphhh)["bcajki"],
    (*Dpphh)["abij"]  <<= (0.0) * (*Tpphh)["abij"],
    (*Dpphh)["abij"]   += (*Fepsp)["a"],
    (*Dpphh)["abij"]   += (*Fepsp)["b"],
    (*Dpphh)["abij"]   -= (*Fepsh)["i"],
    (*Dpphh)["abij"]   -= (*Fepsh)["j"],
    (*Dph)["ai"]  <<= (0.0) * (*Tph)["ai"],
    (*Dph)["ai"]   += (*Fepsp)["a"],
    (*Dph)["ai"]   -= (*Fepsh)["i"]
  )->execute();
  Dpphh->sum(F(1.0), "abij", "abij", [] (F delta) { return F(-1) / delta; });
  Dph->sum(F(1.0), "ai", "ai", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*T2pphh)["abij"] <<= (*T2pphh)["abij"] * (*Dpphh)["abij"],
    (*T2ph)["ai"] <<= (*T2ph)["ai"] * (*Dph)["ai"]
  )->execute();


  r3_VT2T1_2nd(intermediates, T2ph, T2pphh, Rppphhh, t3t3mode);
  auto E4x = getEnergy(intermediates, Rppphhh, "T4-abdghj | <0|L3^[1][H,T1^[2]+T2^[2]]|0> ");



  auto E4a( Tcc<TE>::template tensor<F>("E") );
  auto X4a( Tcc<TE>::template tensor<F>("Ex"));
  auto R4a( Tcc<TE>::template tensor<F>("E") );

  // THIS IS THE T-4 story
  COMPILE(
    (*R4a)["abij"] <<=         (*Whhpp)["klcd"] * (*Tppphhh)["abcijk"] * (*T2ph)["dl"],
    (*R4a)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbijk"] * (*T2ph)["dl"],
    (*R4a)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["acbikl"] * (*T2ph)["dj"],
    (*R4a)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cabikl"] * (*T2ph)["dj"],
    (*R4a)["abij"] += (-1.0) * (*Whhpp)["klcd"] * (*Tppphhh)["adcijk"] * (*T2ph)["bl"],
    (*R4a)["abij"] +=          (*Vhhpp)["klcd"] * (*Tppphhh)["cdaijk"] * (*T2ph)["bl"],
    (*R4a)["abij"] += (*R4a)["baji"],

    (*E4a)[""] <<= 2.0 * (*R4a)["abij"] * (*Zpphh)["abij"],
    (*E4a)[""]  -= 1.0 * (*R4a)["abij"] * (*Zpphh)["baij"],
    (*X4a)[""] <<= 2.0 * (*R4a)["abij"] * (*Tpphh)["abij"],
    (*X4a)[""]  -= 1.0 * (*R4a)["abij"] * (*Tpphh)["baij"]

  )->execute();
  //OUT() << "<0|L2[[H,T3^[1]],T1^[2]]|0>:  " << real(E4a->read()) << std::endl;

  /* R3 <--- [VT3T2] contributions. T5afg */
  r3_VT3T1_a(intermediates, Tph, Tpphh, Tppphhh, Rppphhh, t3t3mode);

  // prepare for T4 (trip)
  COMPILE(
    (*Ippphhh)["abcijk"] <<= (*Rppphhh)["abcijk"]
  )->execute();
  /* T3T3 coupling or (T-3) contributions */
  auto Et3a = getEnergy(intermediates, Rppphhh, "T3-cfi");

  /* R3 <--- [VT3T2] contributions. T5b-T5e */
  r3_VT3T2_b(intermediates, Tph, Tpphh, Tppphhh, Rppphhh, t3t3mode);

  // prepare for T4 (trip)
  COMPILE(
    (*Ippphhh)["abcijk"] += (*Rppphhh)["abcijk"]
  )->execute();
  /* T3T3 coupling or (T-3) contributions */
  auto Et3b = getEnergy(intermediates, Rppphhh, "T3-e");

  std::streamsize ss(std::cout.precision());

  // T3[2]
  // Ok this also matches with cfour.
  getAmplitudes(intermediates, Ippphhh, Tppphhh);

  r3_VT3T1_a(intermediates, Tph, Tpphh, Tppphhh, Rppphhh, t3t3mode);
  auto T4b1 = getEnergy(intermediates, Rppphhh, "T4-cfi | <0|L3^[1][H,T3^[2]]|0>:");

  r3_VT3T2_b(intermediates, Tph, Tpphh, Tppphhh, Rppphhh, t3t3mode);
  auto T4b2 = getEnergy(intermediates, Rppphhh, "T4-e | <0|L3^[1][H,T3^[2]]|0>:");
  //OUT() << "<0|L3^[1][H,T3^[2]]|0>: " << T4b1[1] + T4b2[1] << '\n';

  OUT() << std::setprecision(10);
  auto energy = New<MapNode>(SOURCE_LOCATION);


  if (t3t3mode == "3CC") {
    auto e3cc(Et2[0] + Et3a[0] + Et3b[0]);
    auto e3ccl(Et2[1] + Et3a[1] + Et3b[1]);
    OUT() << "(3CC)_T correlation energy: " << e3cc << '\n';
    OUT() << "(3CC) correlation energy: " << e3ccl << '\n';
    energy->setValue("(3CC)_T", e3cc);
    energy->setValue("(3CC)", e3ccl);
  } else if (t3t3mode == "3CC-C1" || t3t3mode == "3CC-C2" || t3t3mode == "3CC-D1" || t3t3mode == "3CC-T2") {
    auto e3cc(Et2[0] + Et3a[0] + Et3b[0]);
    OUT() << t3t3mode << " correlation energy: " << e3cc << '\n';
    energy->setValue("partial-(3CC)_nonit", e3cc);
  } else {
    auto et2(Et2[0]);
    auto et3(et2 + Et3a[0] + Et3b[0]);
    auto et4(et3 + E4x[0] + real(X4a->read()) + T4b1[0] + T4b2[0]);
    OUT() << "(T-2)_T correlation energy: " << et2 << '\n';
    OUT() << "(T-3)_T correlation energy: " << et3 << '\n';
    OUT() << "(T-4)_T correlation energy: " << et4 << '\n';
    auto et2l(Et2[1]);
    auto et3l(et2 + Et3a[1] + Et3b[1]);
    auto et4l(et3 + E4x[1] + real(E4a->read()) + T4b1[1] + T4b2[1]);

    OUT() << "(T-2) correlation energy: " << et2l << '\n';
    OUT() << "(T-3) correlation energy: " << et3l << '\n';
    OUT() << "(T-4) correlation energy: " << et4l << std::endl;

    energy->setValue("(T-2)_T", et2);
    energy->setValue("(T-3)_T", et3);
    energy->setValue("(T-4)_T", et4);
    energy->setValue("(T-2)", et2l);
    energy->setValue("(T-3)", et3l);
    energy->setValue("(T-4)", et4l);
  }

  OUT() << std::setprecision(ss) << defaultfloat;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;


}



template <typename F, typename TE>
std::vector<Real<>> CcsdtNonit::getEnergy(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                                          Ptr<Tensor<F,TE>> Rppphhh,
                                          std::string name) {

  double _tenergy(0.0), _lenergy(0.0);
  auto Dph( Tcc<TE>::template tensor<F>("Dph") );
  auto Dpphh( Tcc<TE>::template tensor<F>("Dpphh") );
  auto Xpphh( Tcc<TE>::template tensor<F>("Xpphh") );
  auto Xph( Tcc<TE>::template tensor<F>("Xph") );
  auto tenergy( Tcc<TE>::template tensor<F>("Et") );
  auto lenergy( Tcc<TE>::template tensor<F>("El") );
  auto venergy( Tcc<TE>::template tensor<F>("Ev") );

  COMPILE(
    (*Dph)["ai"] <<= 0.0 * (*ints["Tpphh"])["aaii"],
    (*Dpphh)["abij"] <<= 0.0 * (*ints["Tpphh"])["abij"],
    (*Dph)["ai"]  += (*ints["Fepsp"])["a"],
    (*Dph)["ai"]  -= (*ints["Fepsh"])["i"],
    (*Dpphh)["abij"]  += (*ints["Fepsp"])["a"],
    (*Dpphh)["abij"]  += (*ints["Fepsp"])["b"],
    (*Dpphh)["abij"]  -= (*ints["Fepsh"])["i"],
    (*Dpphh)["abij"]  -= (*ints["Fepsh"])["j"],
    (*ints["Xp3h3"])["abcijk"] <<= (*ints["Fepsp"])["a"],
    (*ints["Xp3h3"])["abcijk"]  += (*ints["Fepsp"])["b"],
    (*ints["Xp3h3"])["abcijk"]  += (*ints["Fepsp"])["c"],
    (*ints["Xp3h3"])["abcijk"]  -= (*ints["Fepsh"])["i"],
    (*ints["Xp3h3"])["abcijk"]  -= (*ints["Fepsh"])["j"],
    (*ints["Xp3h3"])["abcijk"]  -= (*ints["Fepsh"])["k"]
  )->execute();
  ints["Xp3h3"]->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  Dpphh->sum(F(1.0), "abij", "abij", [](F delta) { return F(-1) / delta; });
  Dph->sum(F(1.0), "ai", "ai", [](F delta) { return F(-1) / delta; });
  COMPILE(
    (*ints["Xp3h3"])["abcijk"] <<= (*Rppphhh)["abcijk"] * (*ints["Xp3h3"])["abcijk"]
  )->execute();


  // T3 -> R2
  COMPILE(
    // Pure T3->R2
    (*Xpphh)["abij"] <<=          (*ints["Wphpp"])["akcd"] * (*ints["Xp3h3"])["cbdijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*ints["Vphpp"])["akcd"] * (*ints["Xp3h3"])["cdbijk"],
    (*Xpphh)["abij"]  += (-1.0) * (*ints["Whhhp"])["klic"] * (*ints["Xp3h3"])["abckjl"],
    (*Xpphh)["abij"]  +=          (*ints["Vhhhp"])["klic"] * (*ints["Xp3h3"])["acbkjl"],
    (*Xpphh)["abij"]  += (*Xpphh)["baji"],
    (*tenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Tpphh"])["abij"],
    (*tenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Tpphh"])["baij"],
    (*lenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Zpphh"])["abij"],
    (*lenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Zpphh"])["baij"],
    (*venergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Vpphh"])["abij"] * (*Dpphh)["abij"],
    (*venergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Vpphh"])["baij"] * (*Dpphh)["abij"]
  )->execute();

  _tenergy += real(tenergy->read());
  _lenergy += real(lenergy->read());
  if (0) {
  OUT() << name << "-T2a (T2/L2/V2): "
        << real(tenergy->read()) << " " << real(lenergy->read()) << " " << real(venergy->read())
        << std::endl;
  }
  COMPILE(
    // T1+T3 -> R2
    (*Xpphh)["abij"] <<=         (*ints["Whhpp"])["klcd"] * (*ints["Xp3h3"])["abcijk"] * (*ints["Tph"])["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*ints["Whhpp"])["klcd"] * (*ints["Xp3h3"])["acbijk"] * (*ints["Tph"])["dl"],
    (*Xpphh)["abij"] += (-1.0) * (*ints["Whhpp"])["klcd"] * (*ints["Xp3h3"])["acbikl"] * (*ints["Tph"])["dj"],
    (*Xpphh)["abij"] +=          (*ints["Vhhpp"])["klcd"] * (*ints["Xp3h3"])["cabikl"] * (*ints["Tph"])["dj"],
    (*Xpphh)["abij"] += (-1.0) * (*ints["Whhpp"])["klcd"] * (*ints["Xp3h3"])["adcijk"] * (*ints["Tph"])["bl"],
    (*Xpphh)["abij"] +=          (*ints["Vhhpp"])["klcd"] * (*ints["Xp3h3"])["cdaijk"] * (*ints["Tph"])["bl"],
    (*Xpphh)["abij"] += (*Xpphh)["baji"],
    (*tenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Tpphh"])["abij"],
    (*tenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Tpphh"])["baij"],
    (*lenergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Zpphh"])["abij"],
    (*lenergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Zpphh"])["baij"],
    (*venergy)[""] <<=  2.0 * (*Xpphh)["abij"] * (*ints["Vpphh"])["abij"] * (*Dpphh)["abij"],
    (*venergy)[""]  += -1.0 * (*Xpphh)["abij"] * (*ints["Vpphh"])["baij"] * (*Dpphh)["abij"]
  )->execute();

  _tenergy += real(tenergy->read());
  _lenergy += real(lenergy->read());
  if (0) {
  OUT() << name << "-T2b (T2/L2/V2): "
        << real(tenergy->read()) << " " << real(lenergy->read()) << " " << real(venergy->read())
        << std::endl;
  }
  // T3 -> R1
  COMPILE(
    (*Xph)["ai"] += (-2.0) * (*ints["Vhhpp"])["jkbc"] * (*ints["Xp3h3"])["bacjki"],
    (*Xph)["ai"] += (+2.0) * (*ints["Vhhpp"])["jkbc"] * (*ints["Xp3h3"])["bcajki"],
    (*Xph)["ai"] += (+1.0) * (*ints["Vhhpp"])["jkcb"] * (*ints["Xp3h3"])["bacjki"],
    (*Xph)["ai"] += (-1.0) * (*ints["Vhhpp"])["jkcb"] * (*ints["Xp3h3"])["bcajki"],
    (*tenergy)[""] <<=  2.0 * (*Xph)["ai"] * (*ints["Tph"])["ai"],
    (*lenergy)[""] <<=  2.0 * (*Xph)["ai"] * (*ints["Zph"])["ai"],
    (*venergy)[""] <<=  2.0 * (*Xph)["ai"] * (*Dph)["ai"] * (*Xph)["bj"] * (*Dph)["bj"] * (*ints["Vpphh"])["abij"],
    (*venergy)[""]  += -1.0 * (*Xph)["ai"] * (*Dph)["ai"] * (*Xph)["bj"] * (*Dph)["bj"] * (*ints["Vpphh"])["baij"]
  )->execute();

  _tenergy += real(tenergy->read());
  _lenergy += real(lenergy->read());

  if (0) {
  OUT() << name << "-T1 (T1/L1/V1): "
        << real(tenergy->read()) << " " << real(lenergy->read()) << " " << real(venergy->read())
        << std::endl;
  }

  return {_tenergy, _lenergy};
}


template <typename F, typename TE>
void CcsdtNonit::getAmplitudes(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                            Ptr<Tensor<F,TE>> Rppphhh,
                            Ptr<Tensor<F,TE>> Tppphhh) {
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*ints["Fepsp"])["a"],
    (*Tppphhh)["abcijk"]  += (*ints["Fepsp"])["b"],
    (*Tppphhh)["abcijk"]  += (*ints["Fepsp"])["c"],
    (*Tppphhh)["abcijk"]  -= (*ints["Fepsh"])["i"],
    (*Tppphhh)["abcijk"]  -= (*ints["Fepsh"])["j"],
    (*Tppphhh)["abcijk"]  -= (*ints["Fepsh"])["k"]
  )->execute();
  Tppphhh->sum(F(1.0), "abcijk", "abcijk", [] (F delta) { return F(-1) / delta; });
  COMPILE(
    (*Tppphhh)["abcijk"] <<= (*Rppphhh)["abcijk"] * (*Tppphhh)["abcijk"]
  )->execute();
}


template <typename F, typename TE>
void CcsdtNonit::r3_VT3T2_b(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                Ptr<Tensor<F,TE>> Tph,
                Ptr<Tensor<F,TE>> Tpphh,
                Ptr<Tensor<F,TE>> Tppphhh,
                Ptr<Tensor<F,TE>> Rppphhh,
                std::string &t3t3mode) {

  auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
  auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );
  auto Xim( Tcc<TE>::template tensor<F>("Xim") );
  auto Xae( Tcc<TE>::template tensor<F>("Xae") );
  COMPILE(
    //  T3 <-- [VT3T2]
    (t3t3mode == "3CC" || t3t3mode == "3CC-C2") ? (
      (*Xabie)["abie"]  <<= ( 0.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["abfimn"],
      Tcc<TE>::sequence()
    ) : (
      (*Xabie)["abie"]  <<= (-2.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["abfimn"],
      (*Xabie)["abie"]   += ( 1.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["abfnmi"],
      (*Xabie)["abie"]   += ( 1.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["abfinm"],
      Tcc<TE>::sequence()
    ),
    (*Xamij)["amij"]  <<= ( 2.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["aefijn"],
    (*Xamij)["amij"]   += (-1.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["feaijn"],
    (*Xamij)["amij"]   += (-1.0) * (*ints["Vhhpp"])["mnef"] * (*Tppphhh)["afeijn"],

    (*ints["Xp3h3"])["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
    (*Rppphhh)["abcijk"] <<= (*ints["Xp3h3"])["abcijk"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["acbikj"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cabkij"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cbakji"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bcajki"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bacjik"],

    (*Xim)["im"]        <<= ( 2.0) * (*ints["Vhphh"])["iemn"] * (*Tph)["en"],
    (*Xim)["im"]         += (-1.0) * (*ints["Vhphh"])["ienm"] * (*Tph)["en"],
    (*Xim)["im"]         +=          (*ints["Whhpp"])["mnef"] * (*Tpphh)["efin"],
    (*Xim)["im"]         +=          (*ints["Whhpp"])["mnef"] * (*Tph)["ei"] * (*Tph)["fn"],
    (*Xae)["ae"]        <<= ( 2.0) * (*ints["Vphpp"])["amef"] * (*Tph)["fm"],
    (*Xae)["ae"]         += (-1.0) * (*ints["Vphpp"])["amfe"] * (*Tph)["fm"],
    (t3t3mode == "3CC" || t3t3mode == "3CC-C1") ? (
      Tcc<TE>::sequence()
    ) : (
      (*Xae)["ae"]         += (-1.0) * (*ints["Whhpp"])["mnef"] * (*Tpphh)["afmn"],
      Tcc<TE>::sequence()
    ),
    (*Xae)["ae"]         += (-1.0) * (*ints["Whhpp"])["mnef"] * (*Tph)["am"] * (*Tph)["fn"],
    (*ints["Xp3h3"])["abcijk"]  <<=          (*Xae)["ae"]     * (*Tppphhh)["ebcijk"],
    (*ints["Xp3h3"])["abcijk"]   += (-1.0) * (*Xim)["im"]     * (*Tppphhh)["abcmjk"],
    (*Rppphhh)["abcijk"] += (*ints["Xp3h3"])["abcijk"],
    (*Rppphhh)["abcijk"] += (*ints["Xp3h3"])["bacjik"],
    (*Rppphhh)["abcijk"] += (*ints["Xp3h3"])["cbakji"]

  )->execute();

}


template <typename F, typename TE>
void CcsdtNonit::r3_VT3T1_a(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                Ptr<Tensor<F,TE>> Tph,
                Ptr<Tensor<F,TE>> Tpphh,
                Ptr<Tensor<F,TE>> Tppphhh,
                Ptr<Tensor<F,TE>> Rppphhh,
                std::string &t3t3mode) {

  auto Xjkmn( Tcc<TE>::template tensor<F>("Xjkmn") );
  auto Xbcef( Tcc<TE>::template tensor<F>("Xbcef") );
  auto Xamei( Tcc<TE>::template tensor<F>("Xamei") );
  auto Xamie( Tcc<TE>::template tensor<F>("Xamie") );


  COMPILE(
    (*Xjkmn)["jkmn"]  <<=          (*ints["Vhhhh"])["jkmn"],
    (*Xbcef)["bcef"]  <<=          (*ints["Vpppp"])["bcef"],
    (*Xamei)["amei"]  <<=          (*ints["Vphph"])["amei"],
    (*Xamie)["amie"]  <<=          (*ints["Vphhp"])["amie"],

    (t3t3mode == "plain") ? (
      Tcc<TE>::sequence()
    ) : (
      (*Xjkmn)["jkmn"]   +=          (*ints["Vhhpp"])["mnef"] * (*Tpphh)["efjk"],
      (*Xjkmn)["jkmn"]   +=          (*ints["Vhhpp"])["mnef"] * (*Tph)["ej"] * (*Tph)["fk"],
      (*Xjkmn)["jkmn"]   +=          (*ints["Vhphh"])["jemn"] * (*Tph)["ek"],
      (*Xjkmn)["jkmn"]   +=          (*ints["Vphhh"])["ekmn"] * (*Tph)["ej"],

      (*Xbcef)["bcef"]   +=          (*ints["Vpphh"])["efmn"] * (*Tpphh)["bcmn"],
      (*Xbcef)["bcef"]   +=          (*ints["Vpphh"])["efmn"] * (*Tph)["bm"] * (*Tph)["cn"],
      (*Xbcef)["bcef"]   += (-1.0) * (*ints["Vphpp"])["bmef"] * (*Tph)["cm"],
      (*Xbcef)["bcef"]   += (-1.0) * (*ints["Vhppp"])["mcef"] * (*Tph)["bm"],
      (t3t3mode == "3CC" || t3t3mode == "3CC-D1") ? (
        Tcc<TE>::sequence()
      ) : (
        (*Xamei)["amei"]   += (-1.0) * (*ints["Vhhpp"])["mnfe"] * (*Tpphh)["fain"],
        Tcc<TE>::sequence()
      ),
      (*Xamei)["amei"]   += (-1.0) * (*ints["Vhhpp"])["mnfe"] * (*Tph)["fi"] * (*Tph)["an"],
      (*Xamei)["amei"]   += (-1.0) * (*ints["Vhhph"])["nmei"] * (*Tph)["an"],
      (*Xamei)["amei"]   +=          (*ints["Vphpp"])["amef"] * (*Tph)["fi"],

      (*Xamie)["amie"]   +=          (*ints["Whhpp"])["mnef"] * (*Tpphh)["afin"],
      (t3t3mode == "3CC" || t3t3mode == "3CC-D1") ? (
        Tcc<TE>::sequence()
      ) : (
        (*Xamie)["amie"]   += (-1.0) * (*ints["Vhhpp"])["mnef"] * (*Tpphh)["fain"],
        Tcc<TE>::sequence()
      ),
      (*Xamie)["amie"]   += (-1.0) * (*ints["Vhhpp"])["mnef"] * (*Tph)["fi"] * (*Tph)["an"],
      (*Xamie)["amie"]   += (-1.0) * (*ints["Vhhhp"])["nmie"] * (*Tph)["an"],
      (*Xamie)["amie"]   +=          (*ints["Vppph"])["aefm"] * (*Tph)["fi"]
    ),
    (*ints["Xp3h3"])["abcijk"]   <<=          (*Xjkmn)["jkmn"] * (*Tppphhh)["abcimn"],
    (*ints["Xp3h3"])["abcijk"]    +=          (*Xbcef)["bcef"] * (*Tppphhh)["aefijk"],
    (*ints["Xp3h3"])["abcijk"]    += ( 2.0) * (*Xamie)["amie"] * (*Tppphhh)["ebcmjk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["becmjk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamie)["amie"] * (*Tppphhh)["cbemjk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamei)["amei"] * (*Tppphhh)["ebcmjk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamei)["bmei"] * (*Tppphhh)["aecmjk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamei)["cmei"] * (*Tppphhh)["abemjk"],
    (*Rppphhh)["abcijk"] <<= (*ints["Xp3h3"])["abcijk"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bacjik"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cbakji"]
  )->execute();

}

template <typename F, typename TE>
void CcsdtNonit::r3_VT2T1(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                Ptr<Tensor<F,TE>> Tph,
                Ptr<Tensor<F,TE>> Tpphh,
                Ptr<Tensor<F,TE>> Rppphhh,
                std::string &t3t3mode) {

  auto Fphpp( Tcc<TE>::template tensor<F>("Fphpp") );
  auto Fphhh( Tcc<TE>::template tensor<F>("Fphhh") );
  auto Fphhp( Tcc<TE>::template tensor<F>("Fphhp") );
  auto Fphph( Tcc<TE>::template tensor<F>("Fphph") );
  auto Fhpph( Tcc<TE>::template tensor<F>("Fhpph") );
  auto Fhphp( Tcc<TE>::template tensor<F>("Fhphp") );
  auto Fhp( Tcc<TE>::template tensor<F>("Fhp") );
  auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
  auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );

  // construct intermediates
  COMPILE(
    (*Fphpp)["amef"]  <<= (*ints["Vphpp"])["amef"],
    (*Fphpp)["amef"]   += (-1.0) * (*ints["Vhhpp"])["nmef"] * (*Tph)["an"],
    (*Fphhh)["eimn"]  <<= (*ints["Vphhh"])["eimn"],
    (*Fphhh)["eimn"]   += (*ints["Vhhpp"])["mnef"] * (*Tph)["fi"],
    (*Fphhp)["amie"]  <<= (*ints["Vphhp"])["amie"],
    (*Fphhp)["amie"]   += (*ints["Vphpp"])["amfe"] * (*Tph)["fi"],
    (*Fphph)["amei"]  <<= (*ints["Vphph"])["amei"],
    (*Fphph)["amei"]   += (*ints["Vphpp"])["amef"] * (*Tph)["fi"],
    (*Fhpph)["ieam"]  <<= (*ints["Vhpph"])["ieam"],
    (*Fhpph)["ieam"]   += (-1.0) * (*ints["Vhphh"])["ienm"] * (*Tph)["an"],
    (*Fhphp)["iema"]  <<= (*ints["Vhphp"])["iema"],
    (*Fhphp)["iema"]   += (-1.0) * (*ints["Vhhhp"])["inme"] * (*Tph)["an"],
    (*Fhp)["me"] <<= (*ints["Whhpp"])["mnef"] * (*Tph)["fn"]
  )->execute();

  //T1+T2+T3->R3
  COMPILE(
    (*Xabie)["abie"]  <<=          (*ints["Vpphp"])["abie"],
    (*Xamij)["amij"]  <<=          (*ints["Vphhh"])["amij"],
    (t3t3mode == "ring") ? (
      (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afim"],
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["aein"]
    ) : (
      Tcc<TE>::sequence()
    )
  )->execute();

  if (t3t3mode == "dressed" || t3t3mode == "3CC" || t3t3mode == "3CC-T2" || t3t3mode == "3CC-D1" || t3t3mode == "3CC-C1" || t3t3mode == "3CC-C2") {
    COMPILE(
      (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*ints["Xpphh"])["abnm"],
      ( t3t3mode == "3CC" || t3t3mode == "3CC-T2" ) ? (
        (*Fphpp)["amef"]  <<= (-1.0) * (*ints["Vhhpp"])["nmef"] * (*ints["Tph"])["an"],
        (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afim"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*ints["Tpphh"])["afim"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afmi"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*ints["Tpphh"])["bfmi"],
        (*Fphpp)["amef"]   += (*ints["Vphpp"])["amef"],
        Tcc<TE>::sequence()
      ) : (
        (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afim"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*ints["Tpphh"])["afim"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afmi"],
        (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*ints["Tpphh"])["bfmi"],
        Tcc<TE>::sequence()
      ),
      (*Xabie)["abie"]   += (-1.0) * (*Fphhp)["amie"] * (*ints["Tph"])["bm"],
      (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*ints["Tph"])["am"],
      (*Xabie)["abie"]   += ( 1.0) * (*ints["Vpppp"])["abfe"] * (*Tph)["fi"],

      (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*ints["Xpphh"])["efij"],
      (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["aein"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*ints["Tpphh"])["aein"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["eain"],
      (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*ints["Tpphh"])["eajn"],
      (*Xamij)["amij"]   +=          (*Fhpph)["ieam"] * (*ints["Tph"])["ej"],
      (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*ints["Tph"])["ei"],
      (*Xamij)["amij"]   += (-1.0) * (*ints["Vhhhh"])["ijnm"] * (*Tph)["an"]
      //(*Xamij)["amij"]   +=          (*Fhp)["me"] * (*ints["Tpphh"])["aeij"]
    )->execute();
  }

  COMPILE(
    // Construct the R3 from T2 only
    (*ints["Xp3h3"])["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Fhp)["me"] * (*ints["Tpphh"])["aeij"] * (*Tpphh)["bcmk"],
    (*Rppphhh)["abcijk"] <<= (*ints["Xp3h3"])["abcijk"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["acbikj"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cabkij"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cbakji"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bcajki"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bacjik"]
  )->execute();

}



// This is used for T1^[2] && T2^[2] in (T-4)-theory
// Now the input arguments are second order amplitudes.
// Dressing is still done with the convention CCSD amplitudes
template <typename F, typename TE>
void CcsdtNonit::r3_VT2T1_2nd(std::map<std::string, Ptr<TensorExpression<F,TE>>> ints,
                Ptr<Tensor<F,TE>> Tph,
                Ptr<Tensor<F,TE>> Tpphh,
                Ptr<Tensor<F,TE>> Rppphhh,
                std::string &t3t3mode) {

  auto Fphpp( Tcc<TE>::template tensor<F>("Fphpp") );
  auto Fphhh( Tcc<TE>::template tensor<F>("Fphhh") );
  auto Fphhp( Tcc<TE>::template tensor<F>("Fphhp") );
  auto Fphph( Tcc<TE>::template tensor<F>("Fphph") );
  auto Fhpph( Tcc<TE>::template tensor<F>("Fhpph") );
  auto Fhphp( Tcc<TE>::template tensor<F>("Fhphp") );
  auto Fhp( Tcc<TE>::template tensor<F>("Fhp") );
  auto Xabie( Tcc<TE>::template tensor<F>("Xabie") );
  auto Xamij( Tcc<TE>::template tensor<F>("Xamij") );

  // construct intermediates
  COMPILE(
    (*Fphpp)["amef"]  <<= (*ints["Vphpp"])["amef"],
    (*Fphpp)["amef"]   += (-1.0) * (*ints["Vhhpp"])["nmef"] * (*ints["Tph"])["an"],
    (*Fphhh)["eimn"]  <<= (*ints["Vphhh"])["eimn"],
    (*Fphhh)["eimn"]   += (*ints["Vhhpp"])["mnef"] * (*ints["Tph"])["fi"],
    (*Fphhp)["amie"]  <<= (*ints["Vphhp"])["amie"],
    (*Fphhp)["amie"]   += (*ints["Vphpp"])["amfe"] * (*ints["Tph"])["fi"],
    (*Fphph)["amei"]  <<= (*ints["Vphph"])["amei"],
    (*Fphph)["amei"]   += (*ints["Vphpp"])["amef"] * (*ints["Tph"])["fi"],
    (*Fhpph)["ieam"]  <<= (*ints["Vhpph"])["ieam"],
    (*Fhpph)["ieam"]   += (-1.0) * (*ints["Vhphh"])["ienm"] * (*ints["Tph"])["an"],
    (*Fhphp)["iema"]  <<= (*ints["Vhphp"])["iema"],
    (*Fhphp)["iema"]   += (-1.0) * (*ints["Vhhhp"])["inme"] * (*ints["Tph"])["an"],
    (*Fhp)["me"] <<= (*ints["Whhpp"])["mnef"] * (*ints["Tph"])["fn"]
  )->execute();

  //T1+T2+T3->R3
  COMPILE(
    (*Xabie)["abie"]  <<=          (*ints["Vpphp"])["abie"],
    (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*ints["Xpphh"])["abnm"],
    (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*ints["Tpphh"])["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afmi"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*ints["Tpphh"])["bfmi"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphhp)["amie"] * (*ints["Tph"])["bm"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*ints["Tph"])["am"],
    (*Xabie)["abie"]   += ( 1.0) * (*ints["Vpppp"])["abfe"] * (*ints["Tph"])["fi"],

    (*Xamij)["amij"]  <<=          (*ints["Vphhh"])["amij"],
    (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*ints["Xpphh"])["efij"],
    (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*ints["Tpphh"])["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["eain"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*ints["Tpphh"])["eajn"],
    (*Xamij)["amij"]   +=          (*Fhpph)["ieam"] * (*ints["Tph"])["ej"],
    (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*ints["Tph"])["ei"],
    (*Xamij)["amij"]   += (-1.0) * (*ints["Vhhhh"])["ijnm"] * (*ints["Tph"])["an"],
    (*Xamij)["amij"]   +=          (*Fhp)["me"] * (*ints["Tpphh"])["aeij"],

    // add that two body part Hbar: contract T2^[2] with intermediate
    (*ints["Xp3h3"])["abcijk"]   <<= ( 1.0) * (*Xabie)["abie"] * (*Tpphh)["cekj"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*Tpphh)["bcmk"],

    // three body part of Hbar with T2^[2]
    (*Xabie)["abie"]  <<=          (*Fphhh)["eimn"] * (*Tpphh)["abnm"],
    (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*Tpphh)["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*Tpphh)["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*Tpphh)["afmi"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*Tpphh)["bfmi"],

    (*Xamij)["amij"]  <<=          (*Fphpp)["amef"] * (*Tpphh)["efij"],
    (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*Tpphh)["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*Tpphh)["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*Tpphh)["eain"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*Tpphh)["eajn"],
    (*Xamij)["amij"]   +=          (*Fhp)["me"] * (*Tpphh)["aeij"],

    (*ints["Xp3h3"])["abcijk"]    += ( 1.0) * (*Xabie)["abie"] * (*ints["Tpphh"])["cekj"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*ints["Tpphh"])["bcmk"],


    // now we need the three body part containing T1^[2]:
    //

    (*Xabie)["abie"]   <<= (-1.0) * (*Fphhp)["amie"] * (*Tph)["bm"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*Tph)["am"],
    (*Xabie)["abie"]   += ( 1.0) * (*ints["Vpppp"])["abfe"] * (*Tph)["fi"],
//        (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*ints["Xpphh"])["abnm"],
     (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*Tph)["an"] * (*ints["Tph"])["bm"],
     (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*Tph)["bm"] * (*ints["Tph"])["an"],


    (*Xamij)["amij"]   <<=          (*Fhpph)["ieam"] * (*Tph)["ej"],
    (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*Tph)["ei"],
    (*Xamij)["amij"]   += (-1.0) * (*ints["Vhhhh"])["ijnm"] * (*Tph)["an"],
//        (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*ints["Xpphh"])["efij"],
    (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*Tph)["ei"] * (*ints["Tph"])["fj"],
    (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*Tph)["fj"] * (*ints["Tph"])["ei"],




    (*Fphpp)["amef"]   <<= (-1.0) * (*ints["Vhhpp"])["nmef"] * (*Tph)["an"],
    (*Fphhh)["eimn"]   <<= (*ints["Vhhpp"])["mnef"] * (*Tph)["fi"],
    (*Fphhp)["amie"]   <<= (*ints["Vphpp"])["amfe"] * (*Tph)["fi"],
    (*Fphph)["amei"]   <<= (*ints["Vphpp"])["amef"] * (*Tph)["fi"],
    (*Fhpph)["ieam"]   <<= (-1.0) * (*ints["Vhphh"])["ienm"] * (*Tph)["an"],
    (*Fhphp)["iema"]   <<= (-1.0) * (*ints["Vhhhp"])["inme"] * (*Tph)["an"],
    (*Fhp)["me"]       <<= (*ints["Whhpp"])["mnef"] * (*Tph)["fn"],

    (*Xabie)["abie"]   +=          (*Fphhh)["eimn"] * (*ints["Xpphh"])["abnm"],
    (*Xabie)["abie"]   += ( 2.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmfe"] * (*ints["Tpphh"])["afim"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["bmef"] * (*ints["Tpphh"])["afmi"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphpp)["amfe"] * (*ints["Tpphh"])["bfmi"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphhp)["amie"] * (*ints["Tph"])["bm"],
    (*Xabie)["abie"]   += (-1.0) * (*Fphph)["bmei"] * (*ints["Tph"])["am"],


    (*Xamij)["amij"]   +=          (*Fphpp)["amef"] * (*ints["Xpphh"])["efij"],
    (*Xamij)["amij"]   += ( 2.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejmn"] * (*ints["Tpphh"])["aein"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["ejnm"] * (*ints["Tpphh"])["eain"],
    (*Xamij)["amij"]   += (-1.0) * (*Fphhh)["eimn"] * (*ints["Tpphh"])["eajn"],
    (*Xamij)["amij"]   +=          (*Fhpph)["ieam"] * (*ints["Tph"])["ej"],
    (*Xamij)["amij"]   +=          (*Fhphp)["jema"] * (*ints["Tph"])["ei"],
    (*Xamij)["amij"]   +=          (*Fhp)["me"] * (*ints["Tpphh"])["aeij"],


    (*ints["Xp3h3"])["abcijk"]    += ( 1.0) * (*Xabie)["abie"] * (*ints["Tpphh"])["cekj"],
    (*ints["Xp3h3"])["abcijk"]    += (-1.0) * (*Xamij)["amij"] * (*ints["Tpphh"])["bcmk"],

    // Construct the R3 from T2 only
    (*Rppphhh)["abcijk"] <<= (*ints["Xp3h3"])["abcijk"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["acbikj"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cabkij"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["cbakji"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bcajki"],
    (*Rppphhh)["abcijk"]  += (*ints["Xp3h3"])["bacjik"]
  )->execute();

}

