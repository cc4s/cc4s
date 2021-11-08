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

#ifndef STATIC_ASSERT_DEFINED
#define STATIC_ASSERT_DEFINED

#include <math/Real.hpp>
#include <math/Complex.hpp>

namespace cc4s {
  template <typename T>
  class StaticAssert {
  public:
    enum { FALSE = false };
  };

  template <typename A, typename B>
  class TypeRelations {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = false };
  };
  template <typename A>
  class TypeRelations<A,A> {
  public:
    enum { EQUALS = true, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <typename A>
  class TypeRelations<A *, A> {
  public:
    enum { EQUALS = false, POINTER_TO = true, CASTABLE_TO = false };
  };

  template <>
  class TypeRelations<int, Real<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<int, Complex<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };
  template <>
  class TypeRelations<Real<64>, Complex<64>> {
  public:
    enum { EQUALS = false, POINTER_TO = false, CASTABLE_TO = true };
  };

  // TODO: 128 bit real and complex
}

#endif

