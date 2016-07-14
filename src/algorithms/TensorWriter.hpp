/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_WRITER_DEFINED 
#define TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class BinaryHeader {
  public:
    static constexpr int MAX_ORDER = 15;
    BinaryHeader(int order_, int *lens_): order(order_) {
      int i(0);
      for ( ; i < order; ++i) lens[i] = lens_[i];
      for ( ; i < MAX_ORDER; ++i) lens[i] = 0;
    }
    int32_t order;
    int32_t lens[MAX_ORDER];
  };

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
    void writeText();
    void writeBinary();
  };
}

#endif

