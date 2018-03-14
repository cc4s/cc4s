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
    int32_t reserved;

    static constexpr char const *MAGIC = "TENS";
    static constexpr int32_t VERSION = 0x09000;
    static constexpr char const *IEEE = "IEEE";

  protected:
    BinaryTensorHeaderBase() {
    }
    BinaryTensorHeaderBase(
      int32_t bytesPerNumber_, int32_t numsPerElement_, int32_t order_
    ):
      version(VERSION),
      bytesPerNumber(bytesPerNumber_),
      numbersPerElement(numsPerElement_),
      order(order_),
      flags(0), reserved(0)
    {
      std::strncpy(magic, MAGIC, sizeof(magic));
      std::strncpy(numberType, IEEE, sizeof(numberType));
    }
  };

  class BinaryTensorHeader: public BinaryTensorHeaderBase {
  public:
    BinaryTensorHeader(): BinaryTensorHeaderBase() { }
    template <typename Real>
    BinaryTensorHeader(
      CTF::Tensor<Real> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Real), 1, T.order
    ) {
    }
    template <typename Real>
    BinaryTensorHeader(
      CTF::Tensor<Complex<Real>> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Complex<Real>), 2, T.order
    ) {
    }
  };

  class BinaryTensorDimensionHeader {
  public:
    int32_t length;
    char indexName[1];
    int8_t flags;
    int16_t reserved;

    BinaryTensorDimensionHeader() {
    }
    BinaryTensorDimensionHeader(
      int32_t length_, char indexName_
    ):
      length(length_)
    {
      indexName[0] = indexName_;
    }
  };
}

#endif

