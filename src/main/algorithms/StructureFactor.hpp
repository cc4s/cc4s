/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STRUCTURE_FACTOR
#define STRUCTURE_FACTOR

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  /**
   * \brief Caclulates Structure Factor and/or Finite Size Corretion
   */
  class StructureFactor: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(StructureFactor)
    /**
     * \brief Calculates Structure Factor
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

    template <typename TE>
    void interpolation(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);
  };

  template <typename F, typename TE>
  class StructureFactorCalculator {
  public:
    Ptr<const TensorUnion<F,TE>> amplitudes;

    static bool run(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);

    void calculate(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);
  protected:
    StructureFactorCalculator(const Ptr<const TensorUnion<F,TE>> &amplitudes);
  };
}

#endif

