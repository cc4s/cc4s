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

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:1]]
#pragma once

#include <atrip/Atrip.hpp>
#include <atrip/Blas.hpp>
#include <atrip/Utils.hpp>
#include <atrip/Acc.hpp>

namespace atrip {
using ABCTuple = std::array<size_t, 3>;
using PartialTuple = std::array<size_t, 2>;
using ABCTuples = std::vector<ABCTuple>;
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Energy][Energy:1]]
template <typename F = double>
__MAYBE_GLOBAL__ void get_energy_distinct(F const epsabc,
                                          size_t const No,
                                          F *const epsi,
                                          F *const Tijk,
                                          F *const Zijk,
                                          EnergyType<F> *energy);

template <typename F = double>
__MAYBE_GLOBAL__ void get_energy_same(F const epsabc,
                                      size_t const No,
                                      F *const epsi,
                                      F *const Tijk,
                                      F *const Zijk,
                                      EnergyType<F> *energy);
// Energy:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Singles%20contribution][Singles
// contribution:1]]
template <typename F = double>
__MAYBE_GLOBAL__ void singles_contribution(size_t No,
                                           size_t Nv,
                                           size_t a,
                                           size_t b,
                                           size_t c,
                                           DataFieldType<F> *const Tph,
                                           DataFieldType<F> *const VABij,
                                           DataFieldType<F> *const VACij,
                                           DataFieldType<F> *const VBCij,
                                           DataFieldType<F> *Zijk);
// Singles contribution:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Doubles%20contribution][Doubles
// contribution:1]]
template <typename F = double>
void doubles_contribution(size_t const No,
                          size_t const Nv
                          // -- VABCI
                          ,
                          DataPtr<F> const VABph,
                          DataPtr<F> const VACph,
                          DataPtr<F> const VBCph,
                          DataPtr<F> const VBAph,
                          DataPtr<F> const VCAph,
                          DataPtr<F> const VCBph
                          // -- VHHHA
                          ,
                          DataPtr<F> const VhhhA,
                          DataPtr<F> const VhhhB,
                          DataPtr<F> const VhhhC
                          // -- TA
                          ,
                          DataPtr<F> const TAphh,
                          DataPtr<F> const TBphh,
                          DataPtr<F> const TCphh
                          // -- TABIJ
                          ,
                          DataPtr<F> const TABhh,
                          DataPtr<F> const TAChh,
                          DataPtr<F> const TBChh
                          // -- TIJK
                          // , DataPtr<F> Tijk
                          ,
                          DataFieldType<F> *Tijk_,
                          // -- tmp buffers
                          DataFieldType<F> *_t_buffer,
                          DataFieldType<F> *_vhhh);
// Doubles contribution:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:1]]
} // namespace atrip
// Epilog:1 ends here
