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

#include <algorithms/CcsdFocalPointBasisSetCorrection.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdFocalPointBasisSetCorrection)

Ptr<MapNode> CcsdFocalPointBasisSetCorrection::run(
  const Ptr<MapNode> &arguments
) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  // multiplex calls to template methods
  bool success(false);
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    success =
      run<Real<>,TE>(arguments, result) ||
      run<Complex<>,TE>(arguments, result);
  } else {
    using TE = DefaultTensorEngine;
    success =
      run<Real<>,TE>(arguments, result) ||
      run<Complex<>,TE>(arguments, result);
  }
  ASSERT(
    success, "unsupported orbitals type in amplitudes"
  );
  return result;
}



template <typename F, typename TE>
bool CcsdFocalPointBasisSetCorrection::run(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {
  using Tr = TensorExpression<Real<>, TE>;
  using T = TensorExpression<F, TE>;
  auto amplitudes( arguments->getPtr<TensorUnion<F,TE>>("amplitudes") );
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  auto DabijNode(arguments->getMap("deltaIntegrals"));
  auto Dabij(DabijNode->getPtr<Tr>("data"));
  auto nijNode(arguments->getMap("nij"));
  auto nij(nijNode->getPtr<Tr>("data"));

  // these should be the cbs estimates for the mp2 pair energies
  auto mp2PairEnergiesNode(arguments->getMap("mp2PairEnergies"));
  auto mp2PairEnergiesCbs(mp2PairEnergiesNode->getPtr<Tr>("data"));

  // these are the mp2 pair energies in the fno/finite/cc4s basis
  auto mp2PairEnergiesFno( Tcc<TE>::template tensor<F>("mp2PairEnergiesFno"));
  auto eMp2( Tcc<TE>::template tensor<F>("Emp2"));
  auto eMp2Cbs( Tcc<TE>:: template tensor<Real<>>("Emp2Cbs"));

  auto eCcsd( Tcc<TE>::template tensor<F>("Eccsd"));

  auto Tabij( Tcc<TE>::template tensor<F>("Tabij"));

  auto gijccd( Tcc<TE>::template tensor<F>("gijccd"));
  auto gijmp2( Tcc<TE>::template tensor<F>("gijmp2"));
  auto geff  ( Tcc<TE>::template tensor<Real<>>("geff"));
  auto deltaEppl( Tcc<TE>::template tensor<Real<>>("deltaEppl"));
  //mp2 amplitudes
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vabij(coulombSlices->getPtr<T>("pphh"));

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<Tr>("h"));
  auto epsa(energySlices->getPtr<Tr>("p"));

  auto No(epsi->inspect()->getLen(0));
  auto Nv(epsa->inspect()->getLen(0));
  auto Mabij(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Mabij")
  );

  auto fromReal( [](Real<> x) { return F(x); } );
  auto inverse( [](F x) { return 1.0 / x; } );
  COMPILE(
  //reconstruct mp2 amplitudes on-the-fly
    (*Mabij)["abij"] <<= map<F>(fromReal, (*epsi)["i"]),
    (*Mabij)["abij"] +=  map<F>(fromReal, (*epsi)["j"]),
    (*Mabij)["abij"] -=  map<F>(fromReal, (*epsa)["a"]),
    (*Mabij)["abij"] -=  map<F>(fromReal, (*epsa)["b"]),

    (*Mabij)["abij"] <<=
      map<F>(conj<F>, (*Vabij)["abij"]) *
      map<F>(inverse, (*Mabij)["abij"]),
  // calculate mp2 pair energies in fno basis
    (*mp2PairEnergiesFno)["ij"] <<= ( 2.0) * (*Mabij)["abij"] * (*Vabij)["abij"],
    (*mp2PairEnergiesFno)["ij"]  += (-1.0) * (*Mabij)["abij"] * (*Vabij)["abji"],
    (*eMp2)[""] <<= (*mp2PairEnergiesFno)["ij"],
    (*eMp2Cbs)[""] <<= (*mp2PairEnergiesCbs)["ij"],

  //ccsd amplitudes
    (*Tabij)["abij"] <<= (*Tpphh)["abij"],
    (*Tabij)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
  // reevalaute ccsd energy
    (*eCcsd)[""] <<= ( 2.0) * (*Tabij)["abij"] * (*Vabij)["abij"],
    (*eCcsd)[""]  += (-1.0) * (*Tabij)["abij"] * (*Vabij)["abji"],
  // evaluate nominator
    (*gijccd)["ij"] <<= map<F>(fromReal, (*Dabij)["abij"]) * (*Tabij)["abij"],
    (*gijmp2)["ij"] <<= map<F>(fromReal, (*Dabij)["abij"]) * (*Mabij)["abij"],
  // divide by <ij|\delta|ij>
    (*gijccd)["ij"] <<= (*gijccd)["ij"] * map<F>(inverse, (*nij)["ij"]),
    (*gijmp2)["ij"] <<= (*gijmp2)["ij"] * map<F>(inverse, (*nij)["ij"]),

  // final multiplicative factor for the cbsEigenEnergies
    (*geff)["ij"] <<= map<Real<>>(real<F>, (*gijccd)["ij"]),
    (*geff)["ij"]  += map<Real<>>(real<F>, (*gijmp2)["ij"]),
    (*geff)["ij"]  += map<Real<>>(real<F>, (*gijmp2)["ij"])
                    * map<Real<>>(real<F>, (*gijccd)["ij"]),

  // construct \Delta Emp2 and scale with geff
    (*mp2PairEnergiesCbs)["ij"] +=
      (-1.0) * map<Real<>>(real<F>, (*mp2PairEnergiesFno)["ij"]),
    (*deltaEppl)[""] <<= (*geff)["ij"] * (*mp2PairEnergiesCbs)["ij"]

  )->execute();

  Real<> Emp2Cbs(eMp2Cbs->read());
  Real<> deltaPsPpl(deltaEppl->read());
  Real<> Emp2(real<F>(eMp2->read()));
  Real<> Eccsd(real<F>(eCcsd->read()));

  OUT() << "ccsd fno:         " << Eccsd << "\n";
  OUT() << "mp2 fno:          " << Emp2 << "\n";
  OUT() << "Delta mp2:        " << Emp2Cbs - Emp2 << "\n";
  OUT() << "ps-ppl:           " << deltaPsPpl << "\n";
  OUT() << "corrected Energy: " << Eccsd - Emp2 + Emp2Cbs + deltaPsPpl << std::endl;

//  result->setPtr("geff", geff);
//  result->setPtr("gijccd", nij);
//  result->setPtr("mp2p", mp2PairEnergies);

  return true;
}
