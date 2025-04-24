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

// [[file:~/cuda/atrip/atrip.org::*Complex%20numbers][Complex numbers:2]]
#include <atrip/Complex.hpp>

namespace atrip {

template <>
float maybe_conjugate(const float a) {
  return a;
}
template <>
double maybe_conjugate(const double a) {
  return a;
}
template <>
Complex maybe_conjugate(const Complex a) {
  return std::conj(a);
}

namespace traits {
template <typename F>
bool is_complex() {
  return false;
}
template <>
bool is_complex<float>() {
  return false;
}
template <>
bool is_complex<double>() {
  return false;
}
template <>
bool is_complex<Complex>() {
  return true;
}
namespace mpi {
template <>
MPI_Datatype datatype_of<float>() {
  return MPI_FLOAT;
}
template <>
MPI_Datatype datatype_of<double>() {
  return MPI_DOUBLE;
}
template <>
MPI_Datatype datatype_of<Complex>() {
  return MPI_DOUBLE_COMPLEX;
}
} // namespace mpi
} // namespace traits

} // namespace atrip
// Complex numbers:2 ends here
