#ifndef PERT_TRIP_STAR_DEFINED
#define PERT_TRIP_STAR_DEFINED

#include <tcc/Tcc.hpp>
#include <TensorSet.hpp>
#include <MathFunctions.hpp>

namespace cc4s {

  /*
   * * What to compute
   *
   *    - MP2 CBS energy :: from the matrix =mp2PairEnergies(ij)=
   *    - MP2     energy :: from the tensors =Vabij=, =epsi= and =epsa=
   *
   *    where
   *
   *      CBS = Complete Basis Set
   *
   *
   * * Inputs
   *
   * | Key                          | Type   | Default |
   * |------------------------------+--------+---------|
   * | coulombIntegrals.slices.pphh | Tensor | No      |
   * | slicedEigenEnergies.slices.h | Tensor | No      |
   * | slicedEigenEnergies.slices.p | Tensor | No      |
   * | mp2PairEnergies.data         | Matrix | No      |
   *
   *
   */

  template <typename F, typename TE>
  Real<> computeCssdPtStar(const Ptr<MapNode> &arguments, Real<> triples) {
    using Tr = TensorExpression<Real<>, TE>;
    using T = TensorExpression<F, TE>;

    auto Vabij
      = arguments
      ->getMap("coulombIntegrals")
      ->getMap("slices")
      ->getPtr<T>("pphh")
      ;

    // these should be the cbs estimates for the mp2 pair energies
    auto mp2PairEnergiesCbs
      = arguments
      ->getMap("mp2PairEnergies")
      ->getPtr<Tr>("data")
      ;

    auto energySlices
      = arguments->getMap("slicedEigenEnergies")
                 ->getMap("slices") ;
    auto epsi(energySlices->getPtr<Tr>("h"))
       , epsa(energySlices->getPtr<Tr>("p"))
       ;
    auto No(epsi->inspect()->getLen(0))
       , Nv(epsa->inspect()->getLen(0))
       ;
    auto Mabij( // Mabij are the two-body amplitudes for MP2
      Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Mabij")
    );

    // these are the mp2 pair energies in the finite/cc4s basis
    auto mp2PairEnergiesNonCbs(Tcc<TE>::template
                               tensor<F>("mp2PairEnergiesNonCbs"));
    auto eMp2( Tcc<TE>::template tensor<F>("Emp2"));
    auto eMp2Cbs( Tcc<TE>:: template tensor<Real<>>("Emp2Cbs"));

    // functions to deal with complex and real transformation
    // and to get the Δε denominator
    auto fromReal( [](Real<> x) { return F(x); } );
    auto inverse( [](F x) { return 1.0 / x; } );

    COMPILE( (*Mabij)["abij"] <<= cc4s::map<F>(fromReal, (*epsi)["i"])
           , (*Mabij)["abij"]  += cc4s::map<F>(fromReal, (*epsi)["j"])
           , (*Mabij)["abij"]  -= cc4s::map<F>(fromReal, (*epsa)["a"])
           , (*Mabij)["abij"]  -= cc4s::map<F>(fromReal, (*epsa)["b"])
           , (*Mabij)["abij"] <<= cc4s::map<F>(cc4s::conj<F>, (*Vabij)["abij"])
                                * cc4s::map<F>(inverse, (*Mabij)["abij"])

           // calculate mp2 pair energies in non-cbs basis
           , (*mp2PairEnergiesNonCbs)["ij"] <<= ( 2.0)
                                              * (*Mabij)["abij"]
                                              * (*Vabij)["abij"]
           , (*mp2PairEnergiesNonCbs)["ij"]  += (-1.0)
                                              * (*Mabij)["abij"]
                                              * (*Vabij)["abji"]
           , (*eMp2)[""]  <<= (*mp2PairEnergiesNonCbs)["ij"]


           // calculate mp2 pair energies in cbs basis
           , (*eMp2Cbs)[""] <<= (*mp2PairEnergiesCbs)["ij"]
           )->execute();

    const Real<>
        mp2Energy = real(eMp2->read())
      , mp2EnergyCBS = real(eMp2Cbs->read())
      , ratio = mp2EnergyCBS / mp2Energy
      , triples_star = ratio * triples
      ;

    return triples_star;
  }

}

#endif
