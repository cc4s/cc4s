// Copyright 2022 Alejandro Gallo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[file:~/cuda/atrip/atrip.org::*Data%20pointer][Data pointer:1]]
#pragma once
#include <atrip/Acc.hpp>
#include <atrip/Complex.hpp>
#include <atrip/Atrip.hpp>

namespace atrip {

template <typename F>
struct DataField;

template <>
struct DataField<float> {
  using type = float;
};

template <>
struct DataField<double> {
  using type = double;
};

#if defined(HAVE_ACC)

template <typename F>
using DataPtr = ACC_DEVICE_PTR;
#  define DataNullPtr 0x00

template <>
struct DataField<Complex> {
  using type = ACC_DOUBLE_COMPLEX;
};

#else

template <typename F>
using DataPtr = F *;
#  define DataNullPtr nullptr

template <>
struct DataField<Complex> {
  using type = Complex;
};

#endif

template <typename F>
using DataFieldType = typename DataField<F>::type;

/// Energy type

template <typename F>
struct EnergyTypeProxy;

#define DEF_ENERGY_TYPE(__F, __type)                                           \
  template <>                                                                  \
  struct EnergyTypeProxy<__F> {                                                \
    using type = __type;                                                       \
  }

DEF_ENERGY_TYPE(float, double);
DEF_ENERGY_TYPE(double, double);
DEF_ENERGY_TYPE(Complex, double);

template <typename F>
using EnergyType = typename EnergyTypeProxy<F>::type;

template <typename F>
using PrecisionType = EnergyType<F>;

#undef DEF_ENERGY_TYPE

} // namespace atrip
// Data pointer:1 ends here
