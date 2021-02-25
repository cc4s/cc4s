/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STRUCTURE_FACTOR
#define STRUCTURE_FACTOR

#include <algorithms/Algorithm.hpp>

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
  protected:
    template <typename F, typename TE>
    int calculateStructureFactor(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);

    template<typename TE>
    void interpolation(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);

  };
}

#endif

