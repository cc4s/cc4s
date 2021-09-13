#ifndef DELTA_INTEGRALS
#define DELTA_INTEGRALS

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  /**
   * \brief Caclulates Integral with the delta-Kernel
   */
  class StructureFactor: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(StructureFactor)
    /**
     * \brief run routine as always
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    Ptr<MapNode> calculateStructureFactor(const Ptr<MapNode> &arguments);
  
  };
}

#endif

