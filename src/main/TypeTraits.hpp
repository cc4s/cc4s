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

#ifndef TYPE_TRAITS_DEFINED
#define TYPE_TRAITS_DEFINED

#include <Integer.hpp>
#include <Real.hpp>
#include <Complex.hpp>

#include <string>

namespace cc4s {
  /**
   * Traits class for tensor element types used in cc4s.
   * It provides type specific information such as type name to
   * be displayed to the user.
   */
  template <typename F>
  class TypeTraits {
  public:
    static std::string getName() {
      return std::string("unspecified type ") + typeid(F).name();
    }
    enum { IS_PRIMITIVE = false };
  };

  template <>
  class TypeTraits<std::string> {
  public:
    static std::string getName() { return "Text"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<bool> {
  public:
    static std::string getName() { return "Boolean"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Integer<32>> {
  public:
    static std::string getName() { return "Integer32"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Natural<32>> {
  public:
    static std::string getName() { return "Natural32"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Integer<64>> {
  public:
    static std::string getName() { return "Integer64"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Natural<64>> {
  public:
    static std::string getName() { return "Natural64"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Integer<128>> {
  public:
    static std::string getName() { return "Integer128"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Natural<128>> {
  public:
    static std::string getName() { return "Natural128"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Real<64>> {
  public:
    static std::string getName() { return "Real64"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Complex<64>> {
  public:
    static std::string getName() { return "Complex64"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Real<128>> {
  public:
    static std::string getName() { return "Real128"; }
    enum { IS_PRIMITIVE = true };
  };
  template <>
  class TypeTraits<Complex<128>> {
  public:
    static std::string getName() { return "Complex128"; }
    enum { IS_PRIMITIVE = true };
  };
}

#endif

