/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_WRITER_DEFINED 
#define TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class TensorWriter: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorWriter);
    TensorWriter(
      std::vector<Argument> const &argumentList
    );
    virtual ~TensorWriter();
    /**
     * \brief Writes the real tensor data given as Data argument to a file.
     */
    virtual void run();

  protected:
    template <typename F>
    void write(const std::string &name);
  };
}

#endif

