#ifndef BASIS_SET_CORRECTION
#define BASIS_SET_CORRECTION

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  /**
   * \brief Caclulates basis set correction scheme called ps-ppl
   */
  class BasisSetCorrection: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(BasisSetCorrection)
    /**
     * \brief run routine 
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

  protected:
    template <typename F, typename TE>
    bool run(const Ptr<MapNode> &arguments, Ptr<MapNode> &result);
  
  };
}

#endif

