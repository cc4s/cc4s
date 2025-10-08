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

#include <algorithms/coupledcluster/LambdaCcsd.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>


//////////////////////////////////////////////////////////////////////
// Restricted lambda ccsd eqs from PySCF
// https://github.com/pyscf/pyscf/blob/master/pyscf/cc/rccsd_lambda.py
//////////////////////////////////////////////////////////////////////

namespace cc4s {


template <typename F, typename TE>
CoupledClusterMethodRegistrar<
  F,TE,LambdaCcsd<F,TE>
> LambdaCcsd<F,TE>::registrar_("LambdaCcsd");


template <typename F, typename TE>
Ptr<TensorSet<F,TE>> LambdaCcsd<F,TE>::getIntermediates(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {

  auto coulombIntegrals(
    this->arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );
  auto tAmplitudes(
    this->arguments->template getPtr<TensorSet<F,TE>>("tAmplitudes")
  );
  auto eigenEnergies(
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
  // t amplitudes
  auto Tph( tAmplitudes->get("ph") );
  auto Tpphh( tAmplitudes->get("pphh") );
  Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");
  auto Thp( Tcc<TE>::template tensor<F>("Thp") );
  auto Thhpp( Tcc<TE>::template tensor<F>("Thhpp") );
  COMPILE(
    (*Thp)["ia"] <<= map<F>(conj<F>, (*Tph)["ai"]),
    (*Thhpp)["ijab"] <<= map<F>(conj<F>, (*Tpphh)["abij"])
  )->execute();

  auto Vhhhh(coulombIntegrals->get("hhhh"));
  auto Vpppp(coulombIntegrals->get("pppp"));
  auto Vhppp(coulombIntegrals->get("hppp"));
  auto Vhpph(coulombIntegrals->get("hpph"));
  auto Vhhpp(coulombIntegrals->get("hhpp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vpphh(coulombIntegrals->get("pphh"));
  auto Vhphp(coulombIntegrals->get("hphp"));

  // Switch to chemists notation
  auto Whhhh(Tcc<TE>::template tensor<F>("Whhhh"));
  auto Whhpp(Tcc<TE>::template tensor<F>("Whhpp"));
  auto Wpppp(Tcc<TE>::template tensor<F>("Wpppp"));
  auto Whppp(Tcc<TE>::template tensor<F>("Whppp"));
  auto Whpph(Tcc<TE>::template tensor<F>("Whpph"));
  auto Whphp(Tcc<TE>::template tensor<F>("Whphp"));
  auto Whphh(Tcc<TE>::template tensor<F>("Whphh"));

  COMPILE(
    (*Whhhh)["ijkl"] <<= map<F>(conj<F>, (*Vhhhh)["ikjl"]),
    (*Wpppp)["abcd"] <<= map<F>(conj<F>, (*Vpppp)["acbd"]),
    (*Whhpp)["ijab"] <<= map<F>(conj<F>, (*Vhphp)["iajb"]),
    (*Whppp)["iabc"] <<= map<F>(conj<F>, (*Vhppp)["ibac"]),
    (*Whpph)["iabj"] <<= map<F>(conj<F>, (*Vhpph)["ibaj"]),
    (*Whphp)["iajb"] <<= map<F>(conj<F>, (*Vhhpp)["ijab"]),
    (*Whphh)["iajk"] <<= map<F>(conj<F>, (*Vhhph)["ijak"])
  )->execute();

  // effective intermediate
  auto woooo( Tcc<TE>::template tensor<F>("woooo") );
  auto wovvo( Tcc<TE>::template tensor<F>("wovvo") );
  auto woVVo( Tcc<TE>::template tensor<F>("woVVo") );
  auto woovo( Tcc<TE>::template tensor<F>("woovo") );
  auto wvvvo( Tcc<TE>::template tensor<F>("wvvvo") );
  auto v1( Tcc<TE>::template tensor<F>("v1") );
  auto v2( Tcc<TE>::template tensor<F>("v2") );
  auto w3( Tcc<TE>::template tensor<F>("w3") );
  auto v4( Tcc<TE>::template tensor<F>("v4") );

  // helpes
  auto tau( Tcc<TE>::template tensor<F>("tau") );
  auto tmp( Tcc<TE>::template tensor<F>("tmp") );
  auto theta( Tcc<TE>::template tensor<F>("theta") );

  // helpers for intermediates
  auto Whphp1( Tcc<TE>::template tensor<F>("Whphp1") );
  auto Whphh1( Tcc<TE>::template tensor<F>("Whphh1") );
  auto v5( Tcc<TE>::template tensor<F>("v5") );
  auto v4OVvo( Tcc<TE>::template tensor<F>("v4OVvo") );
  auto v4oVVo( Tcc<TE>::template tensor<F>("v4oVVo") );
  auto v4ovvo( Tcc<TE>::template tensor<F>("v4ovvo") );
  auto wOVvo( Tcc<TE>::template tensor<F>("wOVvo") );
  auto ovvv( Tcc<TE>::template tensor<F>("ovvv") );

  auto intermediates(
    New<TensorSet<F,TE>>(
      std::map<std::string,Ptr<TensorExpression<F,TE>>>({
        {"v1",v1},
        {"v2",v2},
        {"w3",w3},
        {"v4",v4},
        {"woooo",woooo},
        {"wovvo",wovvo},
        {"woVVo",woVVo},
        {"woovo",woovo},
        {"wvvvo",wvvvo},
        })
    )
  );

  // intermediates
  COMPILE(
    (*tau)["ijab"] <<= (*Thhpp)["ijab"],
    (*tau)["ijab"] +=  (*Thp)["ia"] * (*Thp)["jb"],
    (*Whphp1)["iajb"] <<= (2.0) * (*Whphp)["iajb"],
    (*Whphp1)["iajb"] += (-1.0) * (*Whphp)["ibja"],
    (*Whphh1)["iajk"] <<= (2.0) * (*Whphh)["iajk"],
    (*Whphh1)["iajk"] += (-1.0) * (*Whphh)["jaik"],
    (*v1)["ba"] <<= (-1.0) * (*Whphp1)["jakc"] * (*tau)["jkbc"],
    (*v1)["cc"] += (*Fepsp)["c"],
    (*v2)["ij"] <<= (1.0) * (*Whphp1)["ibkc"] * (*tau)["jkbc"],
    (*v2)["ij"] += (1.0) * (*Whphh1)["kbij"] * (*Thp)["kb"],
    (*v2)["ii"] += (*Fepsh)["i"],
    (*v4)["jb"] <<= (1.0) * (*Whphp1)["jbkc"] * (*Thp)["kc"],
    (*v5)["bj"] <<= (1.0) * (*v4)["kc"] * (*Thp)["kb"] * (*Thp)["jc"],
    (*v5)["bj"] += (-1.0) * (*Whphh1)["lckj"] * (*Thhpp)["klbc"],
    (*woooo)["ikjl"] <<= (1.0) * (*Whphh)["icjl"] * (*Thp)["kc"],
    (*woooo)["iljk"] += (1.0) * (*Whphh)["jcil"] * (*Thp)["kc"],
    (*woooo)["ijkl"] += (*Whhhh)["ijkl"],
    (*woooo)["ikjl"] += (1.0) * (*Whphp)["icjd"] * (*tau)["klcd"],
    (*theta)["ijab"] <<= (2.0) * (*Thhpp)["ijab"],
    (*theta)["ijab"] += (-1.0) * (*Thhpp)["ijba"],
    (*v4OVvo)["jbck"] <<= (1.0) * (*Whphp1)["ldjb"] * (*Thhpp)["klcd"],
    (*v4OVvo)["jbck"] += (-1.0) * (*Whphp)["ldjb"] * (*Thhpp)["kldc"],
    (*v4OVvo)["jbck"] += (*Whpph)["jbck"],
    (*v4oVVo)["jbck"] <<= (1.0) * (*Whphp)["jdlb"] * (*Thhpp)["kldc"],
    (*v4oVVo)["jbck"] += (-1.0) * (*Whhpp)["jkcb"],
    (*v4ovvo)["jbck"] <<= (2.0) * (*v4OVvo)["jbck"],
    (*v4ovvo)["jbck"] += (1.0) * (*v4oVVo)["jbck"],
    (*w3)["ck"] <<= (1.0) * (*v4ovvo)["jbck"] * (*Thp)["jb"],
    (*woovo)["ijck"] <<= (1.0) * (*v4ovvo)["ibck"] * (*Thp)["jb"],
    (*woovo)["ijck"] += (-1.0) * (*woovo)["ikcj"],
    (*woovo)["ikcj"] += (1.0) * (*v4OVvo)["ibck"] * (*Thp)["jb"],
    (*woovo)["ikcj"] += (-1.0) * (*v4oVVo)["ibck"] * (*Thp)["jb"],
    (*woovo)["ikcj"] += map<F>(conj<F>, (*Whphh1)["jcki"]),
    (*woovo)["ikbj"] += (1.0) * (*Whphh1)["lcik"] * (*theta)["jlbc"],
    (*woovo)["ijbk"] += (-1.0) * (*Whphh1)["lcik"] * (*Thhpp)["jlbc"],
    (*woovo)["ijbk"] += (-1.0) * (*Whphh1)["iclk"] * (*Thhpp)["ljbc"],
    (*wvvvo)["back"] <<= (1.0) * (*v4ovvo)["jack"] * (*Thp)["jb"],
    (*wvvvo)["back"] += (-1.0) * (*wvvvo)["cabk"],
    (*wvvvo)["cabk"] += (1.0) * (*v4OVvo)["jack"] * (*Thp)["jb"],
    (*wvvvo)["cabk"] += (-1.0) * (*v4oVVo)["jack"] * (*Thp)["jb"],
    (*wvvvo)["cabk"] += (-1.0) * (*Whphh1)["lajk"] * (*tau)["jlbc"],
    (*wOVvo)["iabj"] <<= (*v4OVvo)["iabj"],
    (*woVVo)["iabj"] <<= (*v4oVVo)["iabj"],
    (*wOVvo)["jbck"] += (-1.0) * (*Whphp)["jbld"] * (*Thp)["kd"] * (*Thp)["lc"],
    (*woVVo)["jbck"] += (1.0) * (*Whphp)["jdlb"] * (*Thp)["kd"] * (*Thp)["lc"],
    (*wOVvo)["jbck"] += (-1.0) * (*Whphh)["jblk"] * (*Thp)["lc"],
    (*woVVo)["jbck"] += (1.0) * (*Whphh)["lbjk"] * (*Thp)["lc"],
    (*wvvvo)["bacj"] += (1.5) * (*Whppp)["kacd"] * (*Thhpp)["kjbd"] ,
    (*wOVvo)["jbck"] += (1.0) * (*Whppp)["jbcd"] * (*Thp)["kd"],
    (*woVVo)["jbck"] += (-1.0) * (*Whppp)["jdcb"] * (*Thp)["kd"],
    (*ovvv)["iabc"] <<= (2.0) * (*Whppp)["iabc"],
    (*ovvv)["iabc"] += (-1.0) * (*Whppp)["icba"],
    (*v1)["ba"] += (1.0) * (*ovvv)["jcba"] * (*Thp)["jc"],
    (*v5)["bj"] += (1.0) * (*ovvv)["kdbc"] * (*Thhpp)["jkcd"],
    (*woovo)["ijck"] += (1.0) * (*ovvv)["idcb"] * (*tau)["jkdb"],
    (*tmp)["cabj"] <<= (1.0) * (*ovvv)["kdca"] * (*theta)["jkbd"],
    (*wvvvo)["abci"] += (-1.0) * (*tmp)["abci"],
    (*wvvvo)["abci"] += (0.5) * (*tmp)["cbai"],
    (*wvvvo)["abci"] += (-1.0) * map<F>(conj<F>, (*ovvv)["icba"]),
    (*w3)["cj"] += (*v5)["cj"],
    (*w3)["cj"] += (1.0) * (*v1)["cb"] * (*Thp)["jb"],
    (*w3)["bk"] += (-1.0) * (*v2)["jk"] * (*Thp)["jb"],
    (*wovvo)["jbck"] <<= (2.0) * (*wOVvo)["jbck"],
    (*wovvo)["jbck"] += (*woVVo)["jbck"]

  )->execute();

  return intermediates;

}

template <typename F, typename TE>
Ptr<TensorSet<F,TE>> LambdaCcsd<F,TE>::getResiduum(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {

  // get intermediates
  if (!intermediates) {
      intermediates = this->getIntermediates(amplitudes);
  }

  auto v1 = intermediates->get("v1");
  auto v2 = intermediates->get("v2");
  auto w3 = intermediates->get("w3");
  auto v4 = intermediates->get("v4");
  auto woooo = intermediates->get("woooo");
  auto wovvo = intermediates->get("wovvo");
  auto woVVo = intermediates->get("woVVo");
  auto woovo = intermediates->get("woovo");
  auto wvvvo = intermediates->get("wvvvo");

  // construct residuum. Shape will be assumed upon first use.
  auto Rhp( Tcc<TE>::template tensor<F>("Rhp") );
  auto Rhhpp( Tcc<TE>::template tensor<F>("Rhhpp") );
  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto residuum(
    New<TensorSet<F,TE>>(
      std::map<std::string,Ptr<TensorExpression<F,TE>>>(
        {{"ph",Rph}, {"pphh",Rpphh}}
      )
    )
  );

  auto coulombIntegrals(
    this->arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );
  auto tAmplitudes(
    this->arguments->template getPtr<TensorSet<F,TE>>("tAmplitudes")
  );
  auto eigenEnergies(
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
  // t amplitudes
  auto Tph( tAmplitudes->get("ph") );
  auto Tpphh( tAmplitudes->get("pphh") );
  Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");
  auto Thp( Tcc<TE>::template tensor<F>("Thp") );
  auto Thhpp( Tcc<TE>::template tensor<F>("Thhpp") );
  COMPILE(
    (*Thp)["ia"] <<= map<F>(conj<F>, (*Tph)["ai"]),
    (*Thhpp)["ijab"] <<= map<F>(conj<F>, (*Tpphh)["abij"])
  )->execute();

  // energy difference
  auto Dhp = Tcc<TE>::template tensor<F>("Dhp");
  auto Dhhpp = Tcc<TE>::template tensor<F>("Dhhpp");

  auto Vhhhh(coulombIntegrals->get("hhhh"));
  auto Vpppp(coulombIntegrals->get("pppp"));
  auto Vhppp(coulombIntegrals->get("hppp"));
  auto Vhpph(coulombIntegrals->get("hpph"));
  auto Vhhpp(coulombIntegrals->get("hhpp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vpphh(coulombIntegrals->get("pphh"));
  auto Vhphp(coulombIntegrals->get("hphp"));

  // Switch to chemists notation
  auto Whhhh(Tcc<TE>::template tensor<F>("Whhhh"));
  auto Whhpp(Tcc<TE>::template tensor<F>("Whhpp"));
  auto Wpppp(Tcc<TE>::template tensor<F>("Wpppp"));
  auto Whppp(Tcc<TE>::template tensor<F>("Whppp"));
  auto Whpph(Tcc<TE>::template tensor<F>("Whpph"));
  auto Whphp(Tcc<TE>::template tensor<F>("Whphp"));
  auto Whphh(Tcc<TE>::template tensor<F>("Whphh"));

  COMPILE(
    (*Whhhh)["ijkl"] <<= map<F>(conj<F>, (*Vhhhh)["ikjl"]),
    (*Wpppp)["abcd"] <<= map<F>(conj<F>, (*Vpppp)["acbd"]),
    (*Whhpp)["ijab"] <<= map<F>(conj<F>, (*Vhphp)["iajb"]),
    (*Whppp)["iabc"] <<= map<F>(conj<F>, (*Vhppp)["ibac"]),
    (*Whpph)["iabj"] <<= map<F>(conj<F>, (*Vhpph)["ibaj"]),
    (*Whphp)["iajb"] <<= map<F>(conj<F>, (*Vhhpp)["ijab"]),
    (*Whphh)["iajk"] <<= map<F>(conj<F>, (*Vhhph)["ijak"])
  )->execute();


  auto Zhp( Tcc<TE>::template tensor<F>("Zhp") );
  auto Zhhpp( Tcc<TE>::template tensor<F>("Zhhpp") );
  if (!amplitudes) {
    COMPILE(
      (*Zhp)["ia"] <<= (0.0) * (*Thp)["ia"],
      (*Zhhpp)["ijab"] <<= (0.0) * (*Thhpp)["ijab"]
    )->execute();
  } else {
    auto _Zph(amplitudes->get("ph"));
    auto _Zpphh(amplitudes->get("pphh"));
    COMPILE(
      (*Zhp)["ia"] <<= (*_Zph)["ai"],
      (*Zhhpp)["ijab"] <<= (*_Zpphh)["abij"]
    )->execute();
  }
  // helpes
  auto tau( Tcc<TE>::template tensor<F>("tau") );
  auto tmp1( Tcc<TE>::template tensor<F>("tmp1") );
  auto tmp2( Tcc<TE>::template tensor<F>("tmp2") );
  auto theta( Tcc<TE>::template tensor<F>("theta") );

  // helpers fpr residuum eqs
  auto moo( Tcc<TE>::template tensor<F>("moo") );
  auto moo1( Tcc<TE>::template tensor<F>("moo1") );
  auto mvv( Tcc<TE>::template tensor<F>("mvv") );
  auto mvv1( Tcc<TE>::template tensor<F>("mvv1") );
  auto m3( Tcc<TE>::template tensor<F>("m3") );
  auto ovov( Tcc<TE>::template tensor<F>("ovov") );
  auto ZhhppThp( Tcc<TE>::template tensor<F>("ZhhppThp") );
  auto Zhhpptheta( Tcc<TE>::template tensor<F>("Zhhpptheta") );
  auto Zhhpptau( Tcc<TE>::template tensor<F>("Zhhpptau") );

  // lambda ccsd residuum eqs
  COMPILE(

    (*tau)["ijab"] <<= (*Thhpp)["ijab"],
    (*tau)["ijab"] +=  (*Thp)["ia"] * (*Thp)["jb"],
    (*theta)["ijab"] <<= (2.0) * (*Thhpp)["ijab"],
    (*theta)["ijab"] += (-1.0) * (*Thhpp)["ijba"],
    (*mvv)["ba"] <<= (1.0) * (*Zhhpp)["klca"] * (*theta)["klcb"],
    (*moo)["ij"] <<= (1.0) * (*Zhhpp)["kicd"] * (*theta)["kjcd"],
    (*mvv1)["bc"] <<= (1.0) * (*Zhp)["jc"] * (*Thp)["jb"],
    (*mvv1)["bc"] += (1.0) * (*mvv)["bc"],
    (*moo1)["ik"] <<= (1.0) * (*Zhp)["ic"] * (*Thp)["kc"],
    (*moo1)["ik"] += (1.0) * (*moo)["ik"],
    (*m3)["ijcd"] <<= (1.0) * (*Zhhpp)["ijab"] * (*Wpppp)["acbd"],
    (*m3)["ijab"] += (1.0) * (*Zhhpp)["klab"] * (*woooo)["ikjl"],
    (*m3)["ijab"] <<= (0.5) * (*m3)["ijab"],
    (*Zhhpptau)["ijkl"] <<= (1.0) * (*Zhhpp)["ijcd"] * (*tau)["klcd"],
    (*m3)["ijab"] += (0.5) * (*Whphp)["kalb"] * (*Zhhpptau)["ijkl"],
    (*Rhhpp)["ijab"] <<= (0.5) * (*Whphp)["iajb"],
    (*Rhhpp)["ijab"] += (1.0) * (*Zhhpp)["ijac"] * (*v1)["cb"],
    (*Rhhpp)["ijab"] += (-1.0) * (*Zhhpp)["ikab"] * (*v2)["jk"],
    (*Rhhpp)["ijab"] += (-1.0) * (*mvv1)["ca"] * (*Whphp)["icjb"],
    (*Rhhpp)["ijab"] += (-1.0) * (*moo1)["ik"] * (*Whphp)["kajb"],
    (*ovov)["iajb"] <<= (2.0) * (*Whphp)["iajb"],
    (*ovov)["iajb"] += (-1.0) * (*Whphp)["ibja"],
    (*Rhp)["ia"] <<= (-1.0) * (*moo)["ik"] * (*v4)["ka"],
    (*Rhp)["ia"] += (-1.0) * (*mvv)["ca"] * (*v4)["ic"],
    (*Rhhpp)["ijab"] += (1.0) * (*Zhp)["ia"] * (*v4)["jb"],
    (*tmp1)["jb"] <<= (*Thp)["jb"],
    (*tmp1)["jb"] += (*Zhp)["kc"] * (*theta)["kjcb"],
    (*tmp1)["jb"] += (-1.0) * (*mvv1)["bd"] * (*Thp)["jd"],
    (*tmp1)["jb"] += (-1.0) * (*moo)["lj"] * (*Thp)["lb"],
    (*Rhp)["ia"] += (1.0) * (*ovov)["jbia"] * (*tmp1)["jb"],
    (*Rhp)["ia"] += (2.0) * (*Whppp)["iacb"] * (*mvv1)["bc"],
    (*Rhp)["ia"] += (-1.0) * (*Whppp)["ibca"] * (*mvv1)["bc"],
    (*Rhhpp)["jiba"] += (1.0) * (*Zhp)["ic"] * (*Whppp)["jbca"],
    (*ZhhppThp)["ijck"] <<= (1.0) * (*Zhhpp)["ijcd"] * (*Thp)["kd"],
    (*m3)["ijab"] += (-1.0) * (*Whppp)["kbca"] * (*ZhhppThp)["ijck"],
    (*Rhhpp)["ijab"] += (*m3)["ijab"],
    (*Rhp)["ia"] += (2.0) * (*m3)["ijab"] * (*Thp)["jb"],
    (*Rhp)["ia"] += (2.0) * (*m3)["jiba"] * (*Thp)["jb"],
    (*Rhp)["ia"] += (-1.0) * (*m3)["ijba"] * (*Thp)["jb"],
    (*Rhp)["ia"] += (-1.0) * (*m3)["jiab"] * (*Thp)["jb"],
    (*Rhp)["ia"] += (-2.0) * (*Whphh)["iajk"] * (*moo1)["kj"],
    (*Rhp)["ia"] += (1.0) * (*Whphh)["jaik"] * (*moo1)["kj"],
    (*Rhhpp)["ijab"] += (-1.0) * (*Zhp)["ka"] * (*Whphh)["jbik"],
    (*Zhhpptheta)["ijab"] <<= (2.0) * (*Zhhpp)["ijab"],
    (*Zhhpptheta)["ijab"] +=  (-1.0) * (*Zhhpp)["ijba"],
    (*Rhhpp)["ijab"] += (0.5) * (*Zhhpptheta)["ikac"] * (*wovvo)["jbck"],
    (*tmp2)["ijab"] <<= (1.0) * (*Zhhpp)["ikca"] * (*woVVo)["jbck"],
    (*Rhhpp)["ijab"] += (0.5) * (*tmp2)["ijab"],
    (*Rhhpp)["ijab"] += (*tmp2)["jiab"],
    (*Rhp)["ia"] += (1.0) * (*Zhp)["ib"] * (*v1)["ba"],
    (*Rhp)["ia"] += (-1.0) * (*Zhp)["ja"] * (*v2)["ij"],
    (*Rhp)["ia"] += (2.0) * (*Zhp)["jb"] * (*Whpph)["iabj"],
    (*Rhp)["ia"] += (-1.0) * (*Zhp)["jb"] * (*Whhpp)["ijba"],
    (*Rhp)["ia"] += (-1.0) * (*Zhhpp)["ijbc"] * (*wvvvo)["bacj"],
    (*Rhp)["ia"] += (-1.0) * (*Zhhpp)["kjca"] * (*woovo)["ijck"],
    (*Rhp)["ia"] += (2.0) * (*Zhhpp)["ijab"] * (*w3)["bj"],
    (*Rhp)["ia"] += (-1.0) * (*Zhhpp)["ijba"] * (*w3)["bj"],
    (*Rhhpp)["ijab"] += (*Rhhpp)["jiba"]

  )->execute();

  COMPILE(

    (*Dhhpp)["ijab"] <<= 0.0 * (*Thhpp)["ijab"],
    (*Dhhpp)["ijab"] -= (*Fepsp)["a"],
    (*Dhhpp)["ijab"] -= (*Fepsp)["b"],
    (*Dhhpp)["ijab"] += (*Fepsh)["i"],
    (*Dhhpp)["ijab"] += (*Fepsh)["j"],
    (*Dhp)["ia"] <<= 0.0 * (*Thp)["ia"],
    (*Dhp)["ia"] -= (*Fepsp)["a"],
    (*Dhp)["ia"] += (*Fepsh)["i"],
    (*Rhp)["ia"] += (*Zhp)["ia"] * (*Dhp)["ia"],
    (*Rhhpp)["ijab"] += (*Zhhpp)["ijab"] * (*Dhhpp)["ijab"],
    (*Rph)["ai"] <<= (*Rhp)["ia"],
    (*Rpphh)["abij"] <<= (*Rhhpp)["ijab"]

  )->execute();


  return residuum;
}

// instantiate
template class cc4s::LambdaCcsd<Real<64>, DefaultDryTensorEngine>;
template class cc4s::LambdaCcsd<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::LambdaCcsd<Real<64>, DefaultTensorEngine>;
template class cc4s::LambdaCcsd<Complex<64>, DefaultTensorEngine>;
} //namespace

