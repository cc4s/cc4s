/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REAL_TENSOR_READER_DEFINED 
#define REAL_TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class RealTensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RealTensorReader);
    RealTensorReader(
      std::vector<Argument> const &argumentList
    );
    /**
     * \brief Reads real tensor data into the tensor specified by Data.
     */
    virtual void run();
  protected:
    template <typename F, typename TE>
    PTR(ESC(tcc::Tensor<F,TE>)) read(const std::string &name);
  };
}

#endif

