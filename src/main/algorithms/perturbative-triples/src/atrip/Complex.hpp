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

// [[file:~/cuda/atrip/atrip.org::*Complex%20numbers][Complex numbers:1]]
#pragma once

#include <complex>

#include <atrip/mpi.hpp>

#include "config.h"
#if defined(HAVE_CUDA)
#  include <cuComplex.h>
#endif

namespace atrip {

using Complex = std::complex<double>;

template <typename F>
F maybe_conjugate(const F);

#if defined(HAVE_CUDA)
cuDoubleComplex &operator+=(cuDoubleComplex &lz, cuDoubleComplex const &rz);
#endif

namespace traits {

template <typename FF>
bool is_complex();

namespace mpi {
template <typename F>
MPI_Datatype datatype_of(void);
}

} // namespace traits

} // namespace atrip
// Complex numbers:1 ends here
