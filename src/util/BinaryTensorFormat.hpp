#ifndef BINARY_TENSOR_FORMAT_DEFINED
#define BINARY_TENSOR_FORMAT_DEFINED

#include <math/Complex.hpp>
#include <cstring>
#include <ctf.hpp>

namespace cc4s {
  class BinaryTensorHeaderBase {
  public:
    char magic[4];
    int32_t version;
    char numberType[4];
    int32_t bytesPerNumber;
    int32_t numbersPerElement;
    int32_t order;
    int32_t flags;
    int32_t zero;

    static constexpr char const *MAGIC = "TENS";
    static constexpr int32_t VERSION = 0x10000;

  protected:
    BinaryTensorHeaderBase() { }
    BinaryTensorHeaderBase(
      int32_t bytesPerNum, int32_t numsPerEl, int32_t o
    ):
      version(VERSION),
      bytesPerNumber(bytesPerNum),
      numbersPerElement(numsPerEl),
      order(o)
    {
      std::strncpy(magic, MAGIC, sizeof(magic));
    }
  };

  class BinaryTensorHeader: public BinaryTensorHeaderBase {
  public:
    BinaryTensorHeader(): BinaryTensorHeaderBase() { }
    BinaryTensorHeader(CTF::Tensor<> const &T): BinaryTensorHeaderBase(
      sizeof(double), 1, T.order
    ) {
    }
  };

  class BinaryTensorDimensionHeader {
  public:
  };
}

#endif

