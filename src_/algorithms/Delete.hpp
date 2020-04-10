/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DELETE_DEFINED 
#define DELETE_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class Delete: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Delete);
    Delete(
      std::vector<Argument> const &argumentList
    );
    virtual ~Delete();
    /**
     * \brief Frees all associated resources of the given tensor
     */
    virtual void run();
    /**
     * \brief Frees all associated dry resources of the given tensor
     */
    virtual void dryRun();
  };
}

#endif

