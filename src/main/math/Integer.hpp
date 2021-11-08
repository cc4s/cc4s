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

