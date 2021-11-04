#ifndef INTEGER_DEFINED
#define INTEGER_DEFINED

// TODO: use configuration for setting default float type sizes in bits
#define DEFAULT_INTEGER_BIT_SIZE 64

namespace cc4s {
  template <int IntegerSize>
  class IntegerTypes;

  template <>
  class IntegerTypes<32> {
  public:
    typedef int32_t signedType;
    typedef uint32_t unsignedType;
  };

  template <>
  class IntegerTypes<64> {
  public:
    typedef int64_t signedType;
    typedef uint64_t unsignedType;
  };

  template <int IntegerSize=DEFAULT_INTEGER_BIT_SIZE>
  using Integer = typename IntegerTypes<IntegerSize>::signedType;

  template <int IntegerSize=DEFAULT_INTEGER_BIT_SIZE>
  using Natural = typename IntegerTypes<IntegerSize>::unsignedType;
}

#endif

