#ifndef FINITE_SIZE_CORRECTION
#define FINITE_SIZE_CORRECTION

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  /**
   * \brief Caclulates different Finite Size corrections
   */
  class FiniteSizeCorrection: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(FiniteSizeCorrection)
    /**
     * \brief runs main routine
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename TE>
    void interpolation(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);
  };

}

#endif

