/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FTOD_RANK_DECOMPOSITION_DEFINED
#define FTOD_RANK_DECOMPOSITION_DEFINED

#include "Algorithm.hpp"
#include <ctf.hpp>

class FtodRankDecomposition: public Algorithm {
  public:
    FtodRankDecomposition(std::vector<Argument const *> const &argumentList);
    virtual ~FtodRankDecomposition();
    virtual std::vector<std::string> getDefaultArgumentOrder() {
      std::vector<std::string> argumentOrder;
      argumentOrder.push_back("chi");
      argumentOrder.push_back("x");
      argumentOrder.push_back("gamma");
      argumentOrder.push_back("rank");
//      argumentOrder.push_back("epsilon");
      return argumentOrder;
    }
    virtual void run();
  int64_t rank;
  CTF::Tensor<> *chiR, *chiI;
};

#endif

