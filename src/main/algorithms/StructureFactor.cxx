#include <algorithms/StructureFactor.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>


using namespace cc4s;

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
    } else {
      return calculateStructureFactor<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (momentumType == "fullGrid") { 
      return calculateStructureFactor<Real<>,TE>(arguments);
    } else {
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
  using TRc = TensorRecipe<Complex<>, TE>;
  using Tc = Tensor<Complex<>, TE>;
  using Tr = Tensor<Real<>, TE>;

  auto coulombVertexSingularVectors = arguments->getMap("coulombVertexSingularVectors");
  auto singularVectors(coulombVertexSingularVectors->getValue<Ptr<Tc>>("data"));

  auto coulombVertex(arguments->getMap("slicedCoulombVertex"));
  auto slices(coulombVertex->getMap("slices"));
  // get input recipes
  auto GammaFph(slices->getValue<Ptr<TRc>>("ph"));
  auto GammaFhp(slices->getValue<Ptr<TRc>>("hp"));
  auto GammaFhh(slices->getValue<Ptr<TRc>>("hh"));
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
  auto Nij    = ( Tcc<TE>::template tensor<Real<>>("Nij"));
  auto Dpphhc = ( Tcc<TE>::template tensor<Complex<>>("Dpphhc"));
  auto Dpphh  = ( Tcc<TE>::template tensor<F>("Dpphh"));

  auto CoulombPotential(arguments->getMap("coulombPotential"));
  auto VofG(CoulombPotential->getValue<Ptr<Tr>>("data"));
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
//TODO: we have to bring in from complex in F 
//    (*Dpphh)["abij"]  <<= map<F>(fromComplex<F>, (*Dpphhc)["abij"]), 
    // HH 
    (*CGhh)["Gij"]   <<= (*GammaGhh)["Gij"] * (*invSqrtCoulombPotential)["G"],
    (*cTCGhh)["Gji"] <<= map<Complex<>>(conj<Complex<>>, (*GammaGhh)["Gij"]),
    (*cTCGhh)["Gji"] <<= (*cTCGhh)["Gji"] * (*invSqrtCoulombPotential)["G"],
    // Nij 
    (*Nijc)["ij"] <<= (*cTCGhh)["Gii"] * (*CGhh)["Gjj"]
//TODO: we have to bring it from complex in F
//    (*Nij)["ij"] <<= map<F>(fromComplex<F>, (*Nijc)["ij"])
  )->execute();

  //prepare T amplitudes
  auto amplitudesNode(
    arguments->get("amplitudes")->toAtom<Ptr<const TensorUnion<F,TE>>>()
  );
  auto StructureFactor( Tcc<TE>::template tensor<Real<>>("StructureFactor"));
  if (amplitudesNode) {
    auto amplitudes(amplitudesNode->value);
    auto Tph( amplitudes->get(0) );
    auto Tpphh( amplitudes->get(1) );
    auto Tabij( Tcc<TE>::template tensor<Complex<>>("Tabij"));

    COMPILE(
      (*Tpphh)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
      (*Tabij)["abij"] <<= map<Complex<>>(toComplex<F>, (*Tpphh)["abij"])
    )->execute();

    auto SofG( Tcc<TE>::template tensor<Complex<>>("SofG"));
    COMPILE(
      (*SofG)["G"] <<= ( 2.0) * (*cTCGph)["Gai"] * (*CGph)["Gbj"] * (*Tabij)["abij"],
      (*SofG)["G"]  += (-1.0) * (*cTCGph)["Gaj"] * (*CGph)["Gbi"] * (*Tabij)["abij"],
      (*StructureFactor)["G"] <<= map<Real<>>(real<Complex<>>, (*SofG)["G"])
    )->execute();
  }


  auto result(New<MapNode>(SOURCE_LOCATION));
  if (amplitudesNode) {
    result->setValue<Ptr<Tensor<Real<>, TE>>>("structureFactor", StructureFactor);
  }
//  result->setValue<Ptr<Tensor<F, TE>>>("deltaIntegrals", Dpphh);
//  result->setValue<Ptr<Tensor<F, TE>>>("nij", Nij);
  return result;
}
