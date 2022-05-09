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

// [[file:../../atrip.org::*Macros][Macros:1]]
#pragma once
#include <functional>
#define ATRIP_BENCHMARK
//#define ATRIP_DONT_SLICE
#ifndef ATRIP_DEBUG
#  define ATRIP_DEBUG 1
#endif
//#define ATRIP_WORKLOAD_DUMP
#define ATRIP_USE_DGEMM
//#define ATRIP_PRINT_TUPLES

#ifndef ATRIP_DEBUG
#define ATRIP_DEBUG 1
#endif

#if ATRIP_DEBUG == 4
#  pragma message("WARNING: You have OCD debugging ABC triples "    \
                  "expect GB of output and consult your therapist")
#  include <dbg.h>
#  define HAVE_OCD
#  define OCD_Barrier(com) MPI_Barrier(com)
#  define WITH_OCD
#  define WITH_ROOT if (atrip::Atrip::rank == 0)
#  define WITH_SPECIAL(r) if (atrip::Atrip::rank == r)
#  define WITH_RANK std::cout << atrip::Atrip::rank << ": "
#  define WITH_CRAZY_DEBUG
#  define WITH_DBG
#  define DBG(...) dbg(__VA_ARGS__)
#elif ATRIP_DEBUG == 3
#  pragma message("WARNING: You have crazy debugging ABC triples,"  \
                  " expect GB of output")
#  include <dbg.h>
#  define OCD_Barrier(com)
#  define WITH_OCD if (false)
#  define WITH_ROOT if (atrip::Atrip::rank == 0)
#  define WITH_SPECIAL(r) if (atrip::Atrip::rank == r)
#  define WITH_RANK std::cout << atrip::Atrip::rank << ": "
#  define WITH_CRAZY_DEBUG
#  define WITH_DBG
#  define DBG(...) dbg(__VA_ARGS__)
#elif ATRIP_DEBUG == 2
#  pragma message("WARNING: You have some debugging info for ABC triples")
#  define OCD_Barrier(com)
#  define WITH_OCD if (false)
#  define WITH_ROOT if (atrip::Atrip::rank == 0)
#  define WITH_SPECIAL(r) if (atrip::Atrip::rank == r)
#  define WITH_RANK std::cout << atrip::Atrip::rank << ": "
#  define WITH_CRAZY_DEBUG if (false)
#  define WITH_DBG
#  define DBG(...) dbg(__VA_ARGS__)
#else
#  define OCD_Barrier(com)
#  define WITH_OCD if (false)
#  define WITH_ROOT if (false)
#  define WITH_SPECIAL(r) if (false)
#  define WITH_RANK if (false) std::cout << atrip::Atrip::rank << ": "
#  define WITH_DBG if (false)
#  define WITH_CRAZY_DEBUG if (false)
#  define DBG(...)
#endif
// Macros:1 ends here

// [[file:../../atrip.org::*Macros][Macros:2]]
#ifndef LOG
#define LOG(level, name) if (Atrip::rank == 0) std::cout << name << ": "
#endif
// Macros:2 ends here

// [[file:../../atrip.org::*Macros][Macros:3]]
#ifdef ATRIP_NO_OUTPUT
#  undef LOG
#  define LOG(level, name) if (false) std::cout << name << ": "
#endif
// Macros:3 ends here

// [[file:../../atrip.org::IterationDescriptor][IterationDescriptor]]
namespace atrip {

  struct IterationDescription;
  using IterationDescriptor = std::function<void(IterationDescription const&)>;
  struct IterationDescription {
    static IterationDescriptor descriptor;
    size_t currentIteration;
    size_t totalIterations;
    double currentElapsedTime;
  };

  void registerIterationDescriptor(IterationDescriptor);

}
// IterationDescriptor ends here
