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

// [[file:../../atrip.org::*Header][Header:1]]
#pragma once
#include <sstream>
#include <string>
#include <map>
#include <mpi.h>


#include <atrip/Utils.hpp>

#define ADD_ATTRIBUTE(_type, _name, _default)   \
  _type _name = _default;                       \
  Input& with_ ## _name(_type i) {              \
    _name = i;                                  \
    return *this;                               \
  }

namespace atrip {

  struct Atrip {

    static int rank;
    static int np;
    static MPI_Comm communicator;
    static Timings chrono;
    static void init(MPI_Comm);

    template <typename F=double>
    struct Input {
      CTF::Tensor<F> *ei = nullptr
                   , *ea = nullptr
                   , *Tph = nullptr
                   , *Tpphh = nullptr
                   , *Vpphh = nullptr
                   , *Vhhhp = nullptr
                   , *Vppph = nullptr
                   ;
      Input& with_epsilon_i(CTF::Tensor<F> * t) { ei = t; return *this; }
      Input& with_epsilon_a(CTF::Tensor<F> * t) { ea = t; return *this; }
      Input& with_Tai(CTF::Tensor<F> * t) { Tph = t; return *this; }
      Input& with_Tabij(CTF::Tensor<F> * t) { Tpphh = t; return *this; }
      Input& with_Vabij(CTF::Tensor<F> * t) { Vpphh = t; return *this; }
      Input& with_Vijka(CTF::Tensor<F> * t) { Vhhhp = t; return *this; }
      Input& with_Vabci(CTF::Tensor<F> * t) { Vppph = t; return *this; }

      enum TuplesDistribution {
        NAIVE,
        GROUP_AND_SORT,
      };

      ADD_ATTRIBUTE(bool, deleteVppph, false)
      ADD_ATTRIBUTE(bool, rankRoundRobin, false)
      ADD_ATTRIBUTE(bool, chrono, false)
      ADD_ATTRIBUTE(bool, barrier, false)
      ADD_ATTRIBUTE(int, maxIterations, 0)
      ADD_ATTRIBUTE(int, iterationMod, -1)
      ADD_ATTRIBUTE(int, percentageMod, -1)
      ADD_ATTRIBUTE(TuplesDistribution, tuplesDistribution, NAIVE)
      ADD_ATTRIBUTE(std::string, checkpointPath, "atrip-checkpoint.yaml")
      ADD_ATTRIBUTE(bool, readCheckpointIfExists, true)
      ADD_ATTRIBUTE(bool, writeCheckpoint, true)
      ADD_ATTRIBUTE(float, checkpointAtPercentage, 10)
      ADD_ATTRIBUTE(size_t, checkpointAtEveryIteration, 0)

    };

    struct Output {
      double energy;
    };
    template <typename F=double>
    static Output run(Input<F> const& in);
  };

}

#undef ADD_ATTRIBUTE
// Header:1 ends here
