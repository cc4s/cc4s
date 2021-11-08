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

#ifndef REAL_DEFINED
#define REAL_DEFINED

#include <quadmath.h>
#include <ostream>


// TODO: use configuration for setting default float type sizes in bits
#define DEFAULT_FLOAT_BIT_SIZE 64
#define MACHINE_FLOAT_BIT_SIZE 64

namespace cc4s {
  template <int FloatSize>
  class FloatTypes;

  template <>
  class FloatTypes<32> {
  public:
    typedef float type;
  };

  template <>
  class FloatTypes<64> {
  public:
    typedef double type;
  };

  template <>
  class FloatTypes<128> {
  public:
    typedef __float128 type;
  };

  template <int FloatSize=DEFAULT_FLOAT_BIT_SIZE>
  using Real = typename FloatTypes<FloatSize>::type;

  // define stream output for quadruple precision numbers
  inline std::basic_ostream<char> &operator <<(
    std::basic_ostream<char> &stream, const cc4s::Real<128> x
  ) {
    char buffer[1024];
    quadmath_snprintf(buffer, sizeof(buffer), "%.36Qe", x);
    return stream << buffer;
  }
}

#endif

