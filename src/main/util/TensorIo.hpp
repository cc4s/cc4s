#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <util/Scanner.hpp>
#include <Data.hpp>

namespace cc4s {
  class TensorIo {
  public:
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static PTR(ESC(Tensor<F,TE>)) readBinary(
      const std::string &fileName
    );
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static PTR(ESC(Tensor<F,TE>)) readText(
      const std::string &fileName,
      const std::string &delimiter = " ",
      const size_t bufferSize = 1024*1024*1024
    );

    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static void writeBinary(
      const std::string &fileName,
      const PTR(ESC(Tensor<F,TE>)) &A
    );
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static void writeText(
      const std::string &fileName,
      const PTR(ESC(Tensor<F,TE>)) &A,
      const std::string &rowIndexOrder, const std::string &columnIndexOrder,
      const std::string &delimiter = " "
    );
  protected:
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static PTR(ESC(Tensor<F,TE>)) readBinaryHeader(
      MPI_File &file, size_t &offset
    );
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    static PTR(ESC(Tensor<F,TE>)) readTextHeader(
      Scanner &scanner
    );
  };
}

#endif

