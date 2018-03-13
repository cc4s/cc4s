#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <util/Scanner.hpp>
#include <ctf.hpp>

namespace cc4s {
  class TensorIo {
  public:
    template <typename F=real, typename T=CTF::Tensor<F>>
    static T *readBinary(std::string const &fileName);
    template <typename F=real, typename T=CTF::Tensor<F>>
    static T *readText(
      std::string const &fileName,
      std::string const &delimiter = " ",
      int64_t const bufferSize = 1024*1024*1024
    );

    template <typename F=real, typename T=CTF::Tensor<F>>
    static void writeBinary(std::string const &fileName, T &A);
    template <typename F=real, typename T=CTF::Tensor<F>>
    static void writeText(
      std::string const &fileName, T &A,
      std::string const &rowIndexOrder, std::string const &columnIndexOrder,
      std::string const &delimiter = " "
    );
  protected:
    template <typename F=real, typename T=CTF::Tensor<F>>
    static T *readBinaryHeader(MPI_File &file, int64_t &offset);
    template <typename F=real, typename T=CTF::Tensor<F>>
    static T *readTextHeader(Scanner &scanner);
  };
}

#endif

