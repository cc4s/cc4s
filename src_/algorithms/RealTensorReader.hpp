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
    virtual ~RealTensorReader();
    /**
     * \brief Reads real tensor data into the tensor Data.
     */
    virtual void run();
  protected:
    template <typename F>
    CTF::Tensor<F> *read(const std::string &name);
  };
}

#endif

