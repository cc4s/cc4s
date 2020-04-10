/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COMPLEX_TENSOR_READER_DEFINED 
#define COMPLEX_TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class ComplexTensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ComplexTensorReader);
    ComplexTensorReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~ComplexTensorReader();
    /**
     * \brief Reads real tensor data into the tensor Data.
     */
    virtual void run();
  };
}

#endif

