/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_READER_DEFINED 
#define TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorReader);
    TensorReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorReader();
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

