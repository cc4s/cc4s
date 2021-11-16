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

#include <algorithms/StructureFactor.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <math/MathFunctions.hpp>

using namespace cc4s;

template <typename F> inline
F projectReal(Complex<> d) {
  return F(std::real(d));
}

ALGORITHM_REGISTRAR_DEFINITION(StructureFactor)

Ptr<MapNode> StructureFactor::run(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));

  auto slicedCoulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto momentumType(
    slicedCoulombVertex->getMap(
      "indices"
    )->getMap("momentum")->getValue<std::string>("type")
  );



  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (momentumType == "halfGrid") {
      return calculateStructureFactor<Real<>, TE>(arguments);
    } else if (momentumType == "fullGrid") {
      return calculateStructureFactor<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (momentumType == "halfGrid") {
      return calculateStructureFactor<Real<>,TE>(arguments);
    } else if (momentumType == "fullgrid") {
      return calculateStructureFactor<Complex<>,TE>(arguments);
    }
  }
  ASSERT(
    false, "either use 'half Grid' or 'fullGrid' for the momentum"
  );
}

// This algorithm works as follows:
// - undo the SVD of the CoulombVertex and bring it back to reciprocal mesh
// - divide te CoulombVertex by the Coulomb potential to obtain the codensities

template <typename F, typename TE>
Ptr<MapNode> StructureFactor::calculateStructureFactor(
  const Ptr<MapNode> &arguments
) {
  using Tc = TensorExpression<Complex<>, TE>;
  using Tr = TensorExpression<Real<>, TE>;

  auto coulombVertexSingularVectors = arguments->getMap("coulombVertexSingularVectors");
  auto singularVectors(coulombVertexSingularVectors->getPtr<Tc>("data"));

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaFph(slices->getPtr<Tc>("ph"));
  auto GammaFhp(slices->getPtr<Tc>("hp"));
  auto GammaFhh(slices->getPtr<Tc>("hh"));
  auto GammaGph( Tcc<TE>::template tensor<Complex<>>("Gph"));
  auto GammaGhp( Tcc<TE>::template tensor<Complex<>>("Ghp"));
  auto GammaGhh( Tcc<TE>::template tensor<Complex<>>("Ghh"));
  COMPILE(
    (*GammaGph)["Gai"] <<= (*GammaFph)["Fai"] * (*singularVectors)["GF"],
    (*GammaGhp)["Gia"] <<= (*GammaFhp)["Fia"] * (*singularVectors)["GF"],
    (*GammaGhh)["Gij"] <<= (*GammaFhh)["Fij"] * (*singularVectors)["GF"]
  )->execute();


  //We have to take out calculate the overlap coefficients Cpq(G) from Γpq(G) by taking
  //out the reciprocal Coulomb kernel
  //Finally the StructureFactor reads: S(G)=Cai(G)*Cbj*(G)*Tabij
  auto CGph   = ( Tcc<TE>::template tensor<Complex<>>("CGph"));
  auto cTCGph = ( Tcc<TE>::template tensor<Complex<>>("cTCGph"));
  auto CGhh   = ( Tcc<TE>::template tensor<Complex<>>("CGhh"));
  auto cTCGhh = ( Tcc<TE>::template tensor<Complex<>>("cTCGhh"));
  auto Nijc   = ( Tcc<TE>::template tensor<Complex<>>("Nijc"));
  auto Nij    = ( Tcc<TE>::template tensor<F>("Nij"));
  auto Dpphhc = ( Tcc<TE>::template tensor<Complex<>>("Dpphhc"));
  auto Dpphh  = ( Tcc<TE>::template tensor<F>("Dpphh"));

  auto CoulombPotential(arguments->getMap("coulombPotential"));
  auto VofG(CoulombPotential->getPtr<Tr>("data"));
  auto invSqrtCoulombPotential( Tcc<TE>::template tensor<Complex<>>
    ("invSqrtCoulombPotential"));

  COMPILE(
    (*invSqrtCoulombPotential)["G"] <<=
      map<Complex<>>(inverseSqrt<Complex<>>, (*VofG)["G"]),
    // PH codensities
    (*CGph)["Gai"]    <<= (*GammaGph)["Gai"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGph)["Gai"]  <<= map<Complex<>>(conj<Complex<>>, (*GammaGhp)["Gia"]),
    (*cTCGph)["Gai"]  <<= (*cTCGph)["Gai"] * (*invSqrtCoulombPotential)["G"],
    (*Dpphhc)["abij"] <<= (*cTCGph)["Gai"] * (*CGph)["Gbj"],
    (*Dpphh)["abij"]  <<= map<F>(projectReal<F>, (*Dpphhc)["abij"]),
    // HH
    (*CGhh)["Gij"]   <<= (*GammaGhh)["Gij"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGhh)["Gji"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGhh)["Gij"]),
    (*cTCGhh)["Gji"] <<= (*cTCGhh)["Gji"] * (*invSqrtCoulombPotential)["G"],
    // Nij
    (*Nijc)["ij"] <<= (*cTCGhh)["Gii"] * (*CGhh)["Gjj"],
    (*Nij)["ij"] <<= map<F>(projectReal<F>, (*Nijc)["ij"])
  )->execute();

  auto StructureFactor( Tcc<TE>::template tensor<Real<>>("StructureFactor"));
  //prepare T amplitudes


  //THESE TWO LINES ARE SEGFAULTING
  auto amplitudes(arguments->getPtr<TensorUnion<F,TE>>("amplitudes"));

  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );
  auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));


  COMPILE(
//    (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
    (*Tabij)["abij"] <<= map<Complex<>>(toComplex<F>, (*Tpphh)["abij"])
  )->execute();

  auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
  COMPILE(
    (*SofG)["G"] <<= ( 2.0) * (*cTCGph)["Gai"] * (*CGph)["Gbj"] * (*Tabij)["abij"],
    (*SofG)["G"]  += (-1.0) * (*cTCGph)["Gaj"] * (*CGph)["Gbi"] * (*Tabij)["abij"],
    (*StructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
  )->execute();

  auto result(New<MapNode>(SOURCE_LOCATION));

  result->setPtr<TensorExpression<Real<>, TE>>("structureFactor", StructureFactor);
  result->setPtr<TensorExpression<F, TE>>("deltaIntegrals", Dpphh);
  result->setPtr<TensorExpression<F, TE>>("nij", Nij);
  return result;
}
