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

#include <algorithms/coupledcluster/LambdaCcsdReference.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>

using namespace cc4s;

template <typename F, typename TE>
CoupledClusterMethodRegistrar<
  F,TE,LambdaCcsdReference<F,TE>
> LambdaCcsdReference<F,TE>::registrar_("LambdaCcsdReference");


template <typename F, typename TE>
Ptr<TensorSet<F,TE>> LambdaCcsdReference<F,TE>::getResiduum(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {
  // construct residuum. Shape will be assumed upon first use.
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
  auto Tph( tAmplitudes->get("ph") );
  auto Tpphh( tAmplitudes->get("pphh") );
  Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");
  auto Vhhhh(coulombIntegrals->get("hhhh"));
  auto Vhhhp(coulombIntegrals->get("hhhp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vhhpp(coulombIntegrals->get("hhpp"));
  auto Vhphh(coulombIntegrals->get("hphh"));
  auto Vhphp(coulombIntegrals->get("hphp"));
  auto Vhpph(coulombIntegrals->get("hpph"));
  auto Vhppp(coulombIntegrals->get("hppp"));
  auto Vphhh(coulombIntegrals->get("phhh"));
  auto Vphhp(coulombIntegrals->get("phhp"));
  auto Vphph(coulombIntegrals->get("phph"));
  auto Vphpp(coulombIntegrals->get("phpp"));
  auto Vpphh(coulombIntegrals->get("pphh"));
  auto Vpphp(coulombIntegrals->get("pphp"));
  auto Vppph(coulombIntegrals->get("ppph"));
  auto Vpppp(coulombIntegrals->get("pppp"));

  auto Zph( Tcc<TE>::template tensor<F>("Zph") );
  auto Zpphh( Tcc<TE>::template tensor<F>("Zpphh") );

  if (!amplitudes) {
    COMPILE(
      (*Zph)["ai"] <<= (*Tph)["ai"],
      (*Zpphh)["abij"] <<= (*Tpphh)["abij"]
    )->execute();
  } else {
    auto _Zph(amplitudes->get("ph"));
    auto _Zpphh(amplitudes->get("pphh"));
    COMPILE(
      (*Zph)["ai"] <<= (*_Zph)["ai"],
      (*Zpphh)["abij"] <<= (*_Zpphh)["abij"]
    )->execute();

  }

  COMPILE(
    (*Rph)["ai"] <<= (0.0) * (*Tph)["ai"],
    (*Rpphh)["abij"] <<= (0.0) * (*Tpphh)["abij"],
    //(*Rph)["ai"] += (1.00) * f["ia"],
    //(*Rph)["ai"] += (1.00) * (*Zph)["pi"] * f["pa"],
    //(*Rph)["ai"] += (-1.00) * (*Zph)["aq"] * f["iq"],
    //(*Rph)["ai"] += (2.00) * (*Zpphh)["apiq"] * f["pq"],
    //(*Rph)["ai"] += (-1.00) * (*Zpphh)["paiq"] * f["pq"],
    (*Rph)["ai"] += (-1.00) * (*Zph)["ps"] * (*Vhphp)["ipsa"],
    (*Rph)["ai"] += (2.00) * (*Zph)["qs"] * (*Vhpph)["iqas"],
    (*Rph)["ai"] += (2.00) * (*Vhhpp)["ijab"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-1.00) * (*Vhhpp)["ijba"] * (*Tph)["bj"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["pqis"] * (*Vhppp)["spqa"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["qpis"] * (*Vhppp)["spqa"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["apsr"] * (*Vhhph)["rips"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pasr"] * (*Vhhph)["rips"],
    //(*Rph)["ai"] += (-1.00) * (*Zph)["bi"] * f["ja"] * (*Tph)["bj"],
    //(*Rph)["ai"] += (-1.00) * (*Zph)["aj"] * f["ib"] * (*Tph)["bj"],
    //(*Rph)["ai"] += (2.00) * (*Zpphh)["apik"] * f["pc"] * (*Tph)["ck"],
    //(*Rph)["ai"] += (-1.00) * (*Zpphh)["paik"] * f["pc"] * (*Tph)["ck"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["aciq"] * f["kq"] * (*Tph)["ck"],
    //(*Rph)["ai"] += (1.00) * (*Zpphh)["caiq"] * f["kq"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zph)["pi"] * (*Vhppp)["jpba"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-1.00) * (*Zph)["pi"] * (*Vhppp)["jpab"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-1.00) * (*Zph)["pj"] * (*Vhppp)["ipba"] * (*Tph)["bj"],
    (*Rph)["ai"] += (2.00) * (*Zph)["pj"] * (*Vhppp)["ipab"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["as"] * (*Vhhph)["ijbs"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["as"] * (*Vhhph)["jibs"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["bs"] * (*Vhhph)["ijas"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["bs"] * (*Vhhph)["jias"] * (*Tph)["bj"],
    //(*Rph)["ai"] += (1.00) * (*Zpphh)["daik"] * f["lc"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["dcik"] * f["la"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["qpik"] * (*Vpppp)["pqac"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["qpik"] * (*Vpppp)["pqca"] * (*Tph)["ck"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["acil"] * f["kd"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (4.00) * (*Zpphh)["adil"] * f["kc"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["dail"] * f["kc"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (1.00) * (*Zpphh)["dcil"] * f["ka"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["apis"] * (*Vhphp)["kpsc"] * (*Tph)["ck"],
    (*Rph)["ai"] += (4.00) * (*Zpphh)["aqis"] * (*Vhpph)["kqcs"] * (*Tph)["ck"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["cpis"] * (*Vhphp)["kpsa"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["cqis"] * (*Vhpph)["kqas"] * (*Tph)["ck"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pais"] * (*Vhphp)["kpsc"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pcis"] * (*Vhphp)["kpsa"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["qais"] * (*Vhpph)["kqcs"] * (*Tph)["ck"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["qcis"] * (*Vhpph)["kqas"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["bi"] * (*Vhhpp)["jkac"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["bi"] * (*Vhhpp)["jkca"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["aj"] * (*Vhhpp)["ikbc"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["aj"] * (*Vhhpp)["ikcb"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (4.00) * (*Zph)["bj"] * (*Vhhpp)["ikac"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["bj"] * (*Vhhpp)["ikca"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["cj"] * (*Vhhpp)["ikab"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["cj"] * (*Vhhpp)["ikba"] * (*Tpphh)["cbkj"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["apks"] * (*Vhphp)["ipsc"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["aqks"] * (*Vhpph)["iqcs"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["paks"] * (*Vhphp)["ipsc"] * (*Tph)["ck"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["qaks"] * (*Vhpph)["iqcs"] * (*Tph)["ck"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["aclk"] * f["id"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (1.00) * (*Zpphh)["adlk"] * f["ic"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["acsr"] * (*Vhhhh)["ikrs"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["acsr"] * (*Vhhhh)["kirs"] * (*Tph)["ck"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["dpik"] * (*Vhppp)["lpca"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["dpik"] * (*Vhppp)["lpac"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["apil"] * (*Vhppp)["kpdc"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (4.00) * (*Zpphh)["apil"] * (*Vhppp)["kpcd"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["dpil"] * (*Vhppp)["kpca"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["dpil"] * (*Vhppp)["kpac"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pail"] * (*Vhppp)["kpdc"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pail"] * (*Vhppp)["kpcd"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pcil"] * (*Vhppp)["kpda"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pcil"] * (*Vhppp)["kpad"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (4.00) * (*Zpphh)["pdil"] * (*Vhppp)["kpca"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pdil"] * (*Vhppp)["kpac"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["adis"] * (*Vhhph)["klcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adis"] * (*Vhhph)["lkcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dais"] * (*Vhhph)["klcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["dais"] * (*Vhhph)["lkcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["dcis"] * (*Vhhph)["klas"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dcis"] * (*Vhhph)["lkas"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["aplk"] * (*Vhppp)["ipdc"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["aplk"] * (*Vhppp)["ipcd"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["cplk"] * (*Vhppp)["ipda"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["cplk"] * (*Vhppp)["ipad"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["dplk"] * (*Vhppp)["ipca"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (4.00) * (*Zpphh)["dplk"] * (*Vhppp)["ipac"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["acls"] * (*Vhhph)["ikds"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["acls"] * (*Vhhph)["kids"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["adls"] * (*Vhhph)["ikcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adls"] * (*Vhhph)["kics"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["cals"] * (*Vhhph)["ikds"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cals"] * (*Vhhph)["kids"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cdls"] * (*Vhhph)["ikas"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["cdls"] * (*Vhhph)["kias"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dals"] * (*Vhhph)["ikcs"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["dals"] * (*Vhhph)["kics"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["dcls"] * (*Vhhph)["ikas"] * (*Tpphh)["cdkl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dcls"] * (*Vhhph)["kias"] * (*Tpphh)["cdkl"],
    //(*Rph)["ai"] += (1.00) * (*Zpphh)["daik"] * f["lc"] * (*Tph)["ck"] * (*Tph)["dl"],
    //(*Rph)["ai"] += (-2.00) * (*Zpphh)["acil"] * f["kd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["bi"] * (*Vhhpp)["jkac"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["bi"] * (*Vhhpp)["jkca"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["aj"] * (*Vhhpp)["ikbc"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["aj"] * (*Vhhpp)["ikcb"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (-2.00) * (*Zph)["cj"] * (*Vhhpp)["ikab"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zph)["cj"] * (*Vhhpp)["ikba"] * (*Tph)["ck"] * (*Tph)["bj"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["dpik"] * (*Vhppp)["lpca"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["dpik"] * (*Vhppp)["lpac"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["apil"] * (*Vhppp)["kpdc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (4.00) * (*Zpphh)["apil"] * (*Vhppp)["kpcd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pail"] * (*Vhppp)["kpdc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pail"] * (*Vhppp)["kpcd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["pcil"] * (*Vhppp)["kpda"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["pcil"] * (*Vhppp)["kpad"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["adis"] * (*Vhhph)["klcs"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adis"] * (*Vhhph)["lkcs"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dais"] * (*Vhhph)["klcs"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["dais"] * (*Vhhph)["lkcs"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["dcis"] * (*Vhhph)["klas"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["dcis"] * (*Vhhph)["lkas"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-2.00) * (*Zpphh)["aplk"] * (*Vhppp)["ipdc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (1.00) * (*Zpphh)["aplk"] * (*Vhppp)["ipcd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["acls"] * (*Vhhph)["ikds"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["acls"] * (*Vhhph)["kids"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["cals"] * (*Vhhph)["ikds"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cals"] * (*Vhhph)["kids"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["aeik"] * (*Vhhpp)["mlcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["aeik"] * (*Vhhpp)["mldc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["eaik"] * (*Vhhpp)["mlcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["eaik"] * (*Vhhpp)["mldc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["edik"] * (*Vhhpp)["mlac"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["edik"] * (*Vhhpp)["mlca"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["eail"] * (*Vhhpp)["kmcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["eail"] * (*Vhhpp)["kmdc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["ecil"] * (*Vhhpp)["kmad"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["ecil"] * (*Vhhpp)["kmda"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["edil"] * (*Vhhpp)["kmac"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["edil"] * (*Vhhpp)["kmca"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["acim"] * (*Vhhpp)["klde"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["acim"] * (*Vhhpp)["kled"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["adim"] * (*Vhhpp)["klce"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adim"] * (*Vhhpp)["klec"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (8.00) * (*Zpphh)["aeim"] * (*Vhhpp)["klcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["aeim"] * (*Vhhpp)["kldc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["caim"] * (*Vhhpp)["klde"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["caim"] * (*Vhhpp)["kled"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cdim"] * (*Vhhpp)["klae"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["cdim"] * (*Vhhpp)["klea"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["ceim"] * (*Vhhpp)["klad"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["ceim"] * (*Vhhpp)["klda"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["eaim"] * (*Vhhpp)["klcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["eaim"] * (*Vhhpp)["kldc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["ecim"] * (*Vhhpp)["klad"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["ecim"] * (*Vhhpp)["klda"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["edim"] * (*Vhhpp)["klac"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["edim"] * (*Vhhpp)["klca"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["eakl"] * (*Vhhpp)["imcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["eakl"] * (*Vhhpp)["imdc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["edkl"] * (*Vhhpp)["imac"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["edkl"] * (*Vhhpp)["imca"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adkm"] * (*Vhhpp)["ilce"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["adkm"] * (*Vhhpp)["ilec"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["aekm"] * (*Vhhpp)["ilcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["aekm"] * (*Vhhpp)["ildc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["eakm"] * (*Vhhpp)["ilcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["eakm"] * (*Vhhpp)["ildc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["edkm"] * (*Vhhpp)["ilac"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["edkm"] * (*Vhhpp)["ilca"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["acml"] * (*Vhhpp)["ikde"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["acml"] * (*Vhhpp)["iked"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adml"] * (*Vhhpp)["ikce"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["adml"] * (*Vhhpp)["ikec"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["aeml"] * (*Vhhpp)["ikcd"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["aeml"] * (*Vhhpp)["ikdc"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["cdml"] * (*Vhhpp)["ikae"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cdml"] * (*Vhhpp)["ikea"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["ceml"] * (*Vhhpp)["ikad"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["ceml"] * (*Vhhpp)["ikda"] * (*Tpphh)["delm"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-4.00) * (*Zpphh)["aeik"] * (*Vhhpp)["mlcd"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["aeik"] * (*Vhhpp)["mldc"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["caim"] * (*Vhhpp)["klde"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["caim"] * (*Vhhpp)["kled"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["cdim"] * (*Vhhpp)["klae"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["cdim"] * (*Vhhpp)["klea"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (2.00) * (*Zpphh)["adkm"] * (*Vhhpp)["ilce"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rph)["ai"] += (-1.00) * (*Zpphh)["adkm"] * (*Vhhpp)["ilec"] * (*Tph)["dl"] * (*Tph)["em"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Vhhpp)["ijab"],
    //(*Rpphh)["abij"] += (1.00) * (*Zph)["ai"] * f["jb"],
    //(*Rpphh)["abij"] += (1.00) * (*Zph)["bj"] * f["ia"],
    //(*Rpphh)["abij"] += (1.00) * (*Zpphh)["apij"] * f["pb"],
    //(*Rpphh)["abij"] += (1.00) * (*Zpphh)["pbij"] * f["pa"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["abiq"] * f["jq"],
    (*Rpphh)["abij"] += (1.00) * (*Zph)["pi"] * (*Vhppp)["jpba"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["bajq"] * f["iq"],
    (*Rpphh)["abij"] += (1.00) * (*Zph)["pj"] * (*Vhppp)["ipab"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["as"] * (*Vhhph)["jibs"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["bs"] * (*Vhhph)["ijas"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["qpij"] * (*Vpppp)["pqba"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["apis"] * (*Vhphp)["jpsb"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["aqis"] * (*Vhpph)["jqbs"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pbis"] * (*Vhphp)["jpsa"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["qais"] * (*Vhpph)["jqbs"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["bpjs"] * (*Vhphp)["ipsa"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["bqjs"] * (*Vhpph)["iqas"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pajs"] * (*Vhphp)["ipsb"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["qbjs"] * (*Vhpph)["iqas"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["absr"] * (*Vhhhh)["jirs"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["acij"] * f["kb"] * (*Tph)["ck"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["cbij"] * f["ka"] * (*Tph)["ck"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["abik"] * f["jc"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (2.00) * (*Zph)["ai"] * (*Vhhpp)["jlbd"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["ai"] * (*Vhhpp)["jldb"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["di"] * (*Vhhpp)["jlba"] * (*Tph)["dl"],
    //(*Rpphh)["abij"] += (-1.00) * (*Zpphh)["bajk"] * f["ic"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (2.00) * (*Zph)["bj"] * (*Vhhpp)["ilad"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["bj"] * (*Vhhpp)["ilda"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["dj"] * (*Vhhpp)["ilab"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["al"] * (*Vhhpp)["ijdb"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-1.00) * (*Zph)["bl"] * (*Vhhpp)["ijad"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["apij"] * (*Vhppp)["kpcb"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["apij"] * (*Vhppp)["kpbc"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["cpij"] * (*Vhppp)["kpab"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["pbij"] * (*Vhppp)["kpca"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pbij"] * (*Vhppp)["kpac"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pcij"] * (*Vhppp)["kpba"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["apik"] * (*Vhppp)["jpcb"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["apik"] * (*Vhppp)["jpbc"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["paik"] * (*Vhppp)["jpbc"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pbik"] * (*Vhppp)["jpca"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["abis"] * (*Vhhph)["jkcs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["abis"] * (*Vhhph)["kjcs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["acis"] * (*Vhhph)["jkbs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["acis"] * (*Vhhph)["kjbs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["cais"] * (*Vhhph)["jkbs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["cbis"] * (*Vhhph)["kjas"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["bpjk"] * (*Vhppp)["ipca"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (2.00) * (*Zpphh)["bpjk"] * (*Vhppp)["ipac"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pajk"] * (*Vhppp)["ipcb"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-1.00) * (*Zpphh)["pbjk"] * (*Vhppp)["ipac"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bajs"] * (*Vhhph)["ikcs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bajs"] * (*Vhhph)["kics"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bcjs"] * (*Vhhph)["ikas"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bcjs"] * (*Vhhph)["kias"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["cajs"] * (*Vhhph)["kibs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["cbjs"] * (*Vhhph)["ikas"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["abks"] * (*Vhhph)["ijcs"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["baks"] * (*Vhhph)["jics"] * (*Tph)["ck"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["adij"] * (*Vhhpp)["lkbc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["adij"] * (*Vhhpp)["lkcb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dbij"] * (*Vhhpp)["lkac"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbij"] * (*Vhhpp)["lkca"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dcij"] * (*Vhhpp)["lkab"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["daik"] * (*Vhhpp)["jlbc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbik"] * (*Vhhpp)["jlca"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dcik"] * (*Vhhpp)["jlba"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["abil"] * (*Vhhpp)["jkcd"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["abil"] * (*Vhhpp)["jkdc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["acil"] * (*Vhhpp)["jkbd"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["acil"] * (*Vhhpp)["jkdb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (4.00) * (*Zpphh)["adil"] * (*Vhhpp)["jkbc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["adil"] * (*Vhhpp)["jkcb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dail"] * (*Vhhpp)["jkbc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dail"] * (*Vhhpp)["jkcb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dcil"] * (*Vhhpp)["jkba"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dajk"] * (*Vhhpp)["ilcb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbjk"] * (*Vhhpp)["ilac"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dcjk"] * (*Vhhpp)["ilab"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bajl"] * (*Vhhpp)["ikcd"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bajl"] * (*Vhhpp)["ikdc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bcjl"] * (*Vhhpp)["ikad"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bcjl"] * (*Vhhpp)["ikda"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (4.00) * (*Zpphh)["bdjl"] * (*Vhhpp)["ikac"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bdjl"] * (*Vhhpp)["ikca"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dbjl"] * (*Vhhpp)["ikac"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbjl"] * (*Vhhpp)["ikca"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dcjl"] * (*Vhhpp)["ikab"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["ablk"] * (*Vhhpp)["ijdc"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["aclk"] * (*Vhhpp)["ijdb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["adlk"] * (*Vhhpp)["ijcb"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bclk"] * (*Vhhpp)["ijad"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bdlk"] * (*Vhhpp)["ijac"] * (*Tpphh)["cdkl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["adij"] * (*Vhhpp)["lkbc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["adij"] * (*Vhhpp)["lkcb"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["dbij"] * (*Vhhpp)["lkac"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbij"] * (*Vhhpp)["lkca"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dcij"] * (*Vhhpp)["lkab"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["daik"] * (*Vhhpp)["jlbc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbik"] * (*Vhhpp)["jlca"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["abil"] * (*Vhhpp)["jkcd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["abil"] * (*Vhhpp)["jkdc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["acil"] * (*Vhhpp)["jkbd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["acil"] * (*Vhhpp)["jkdb"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dajk"] * (*Vhhpp)["ilcb"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["dbjk"] * (*Vhhpp)["ilac"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bajl"] * (*Vhhpp)["ikcd"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bajl"] * (*Vhhpp)["ikdc"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (-2.00) * (*Zpphh)["bcjl"] * (*Vhhpp)["ikad"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["bcjl"] * (*Vhhpp)["ikda"] * (*Tph)["ck"] * (*Tph)["dl"],
    (*Rpphh)["abij"] += (1.00) * (*Zpphh)["ablk"] * (*Vhhpp)["ijdc"] * (*Tph)["ck"] * (*Tph)["dl"]
  )->execute();
/*
// implekdcaation based on JCP 87, 5361 (1987)
// Employed interemediates:
// Iph, gpp, khh, Ipphh, Ehhhp, Fphpp, Hppph, rphhp, qphhp, jhphh
// implementation removed
*/
  return residuum;
}

// instantiate
template class cc4s::LambdaCcsdReference<Real<64>, DefaultDryTensorEngine>;
template class cc4s::LambdaCcsdReference<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::LambdaCcsdReference<Real<64>, DefaultTensorEngine>;
template class cc4s::LambdaCcsdReference<Complex<64>, DefaultTensorEngine>;

