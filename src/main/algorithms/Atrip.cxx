#include <algorithms/Atrip.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <math/MathFunctions.hpp>
#include <atrip.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(Atrip)

Ptr<MapNode> Atrip::run(const Ptr<MapNode> &arguments) {
  using TE = DefaultTensorEngine;

  auto result(New<MapNode>(SOURCE_LOCATION));

  OUT() << "Atrip init.. " << std::endl;
  atrip::Atrip::init();
  atrip::Atrip::Input in;

#define __V__(_idx)                                   \
    &(arguments                                       \
      ->getMap("coulombIntegrals")                    \
      ->getMap("slices")                              \
      ->getValue<Ptr<TensorRecipe<Real<>,TE>>>(#_idx) \
      ->getResult()                                   \
      ->getMachineTensor()                            \
      ->tensor)
#define __T__(_idx) \
   &(arguments                                                    \
      ->getValue<Ptr<const TensorUnion<Real<>,TE>>>("amplitudes") \
      ->get(_idx)                                                 \
      ->getMachineTensor()                                        \
      ->tensor)
#define __eps__(_idx)                                      \
   &(arguments                                             \
      ->getMap("slicedEigenEnergies")                      \
      ->getMap("slices")                                   \
      ->getValue<Ptr<TensorRecipe<Real<>,TE>>>(#_idx) \
      ->getResult()                                        \
      ->getMachineTensor()                                 \
      ->tensor)

  // this is a hack so that a CTF::World gets created for sure
  CTF::World _w(MPI_COMM_WORLD);
  in
    // setup tensors
    .with_epsilon_i(__eps__(h))
    .with_epsilon_a(__eps__(p))
    .with_Tai(__T__(0))
    .with_Tabij(__T__(1))
    .with_Vabij(__V__(pphh))
    .with_Vijka(__V__(hhhp))
    .with_Vabci(__V__(ppph))
    // some options
    .with_barrier(false)
    .with_iterationMod(100)
    ;
#undef __V__
#undef __T__
#undef __eps__

  auto out = atrip::Atrip::run(in);
  LOG() << "Energy: " << out.energy << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("triples", out.energy);
  result->get("energy") = energy;

  return result;
}
