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
#include <atrip/Complex.hpp>
#include <atrip/Atrip.hpp>

namespace atrip {

template <typename F>
struct DataField;

template <>
struct DataField<double> {
  using type = double;
};

#if defined(HAVE_CUDA)

template <typename F>
using DataPtr = CUdeviceptr;
#  define DataNullPtr 0x00

template <>
struct DataField<Complex> {
  using type = cuDoubleComplex;
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

} // namespace atrip
// Data pointer:1 ends here
