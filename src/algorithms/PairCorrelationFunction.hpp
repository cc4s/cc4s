/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REAL_STRUCTURE_DEFINED
#define REAL_STRUCTURE_DEFINED

#include <math/Vector.hpp>

namespace cc4s {
  class RealStructure: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RealStructure);
    RealStructure(
		  std::vector<Argument> const &argumentList
		  );
    virtual ~RealStructure();
    
    virtual void run();

  protected:
    int NG;

  };

}


#endif
