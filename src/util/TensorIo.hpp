#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <ctf.hpp>

namespace cc4s {
  class TensorIo {
  public:
    template <typename F=double>
    static CTF::Tensor<F> *readBinary(std::string const &fileName);
    template <typename F=double>
    static CTF::Tensor<F> *readText(
      std::string const &fileName,
      std::string const &delimiter = " ",
      int64_t const bufferSize = 1024*1024*1024
    );

    template <typename F=double>
    static void writeBinary(std::string const &fileName, CTF::Tensor<F> &A);
    template <typename F=double>
    static void writeText(
      std::string const &fileName, CTF::Tensor<F> &A,
      std::string const &rowIndexOrder, std::string const &columnIndexOrder,
      std::string const &delimiter = " "
    );
  };
}

#endif

