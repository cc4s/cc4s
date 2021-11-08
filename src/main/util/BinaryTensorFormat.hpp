/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef BINARY_TENSOR_FORMAT_DEFINED
#define BINARY_TENSOR_FORMAT_DEFINED

#include <math/Complex.hpp>
#include <cstring>
#include <ctf.hpp>

namespace cc4s {
  class BinaryTensorHeaderBase {
  public:
    char magic[4];
    uint32_t version;
    char numberType[4];
    uint32_t bytesPerNumber;
    uint32_t numbersPerElement;
    uint32_t order;
    uint32_t flags;
    uint32_t reserved;

    static constexpr char const *MAGIC = "TENS";
    static constexpr int32_t VERSION = 0x09000;
    static constexpr char const *IEEE = "IEEE";

  protected:
    BinaryTensorHeaderBase() {
    }
    BinaryTensorHeaderBase(
      uint32_t bytesPerNumber_, uint32_t numsPerElement_, uint32_t order_
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
    BinaryTensorHeader(
      CTF::Tensor<Real<64>> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Real<64>), 1, T.order
    ) {
    }
    BinaryTensorHeader(
      CTF::Tensor<Complex<64>> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Real<64>), 2, T.order
    ) {
    }
    BinaryTensorHeader(
      CTF::Tensor<Real<128>> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Real<128>), 1, T.order
    ) {
    }
    BinaryTensorHeader(
      CTF::Tensor<Complex<128>> const &T
    ): BinaryTensorHeaderBase(
      sizeof(Real<128>), 2, T.order
    ) {
    }
  };

  class BinaryTensorDimensionHeader {
  public:
    uint32_t length;
    char indexName[1];
    uint8_t flags;
    uint16_t reserved;

    BinaryTensorDimensionHeader() {
    }
    BinaryTensorDimensionHeader(
      uint32_t length_, char indexName_
    ):
      length(length_)
    {
      indexName[0] = indexName_;
    }
  };
}

#endif

