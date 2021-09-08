#ifndef _ATRIP_HEADER_FILE
#define _ATRIP_HEADER_FILE

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>

namespace cc4s {
  struct Atrip: public Algorithm {
    ALGORITHM_REGISTRAR_DECLARATION(Atrip)
    /**
     * \brief run routine as always
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  };
}

#endif

