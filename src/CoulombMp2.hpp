/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_MP2_DEFINED 
#define COULOMB_MP2_DEFINED

#include <Algorithm.hpp>

namespace cc4s {
  class CoulombMp2: public Algorithm {
    public:
      CoulombMp2(
        std::vector<Argument const *> const &argumentList
      );
      virtual ~CoulombMp2();
      virtual void run();


  };
}


#endif