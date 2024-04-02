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

// [[file:~/cuda/atrip/atrip.org::*Header][Header:1]]
#ifndef ATRIP_HPP_
#define ATRIP_HPP_

#include <sstream>
#include <string>
#include <map>
#include "config.h"

#include <atrip/Chrono.hpp>
#include <atrip/mpi.hpp>

#if defined(HAVE_CUDA)
#  include <cuda.h>
#  define CUBLASAPI
#  include <cublas_api.h>
#  include <cublas_v2.h>
#endif

#include <atrip/Utils.hpp>
#include <atrip/Types.hpp>
#include <atrip/Tuples.hpp>

#define ADD_ATTRIBUTE(_type, _name, _default)                                  \
  _type _name = _default;                                                      \
  Input &with_##_name(_type i) {                                               \
    _name = i;                                                                 \
    return *this;                                                              \
  }

namespace atrip {

struct Atrip {

  static size_t rank;
  static size_t np;
  static ClusterInfo *cluster_info;
  static MPI_Comm communicator;
  static Timings chrono;
  static size_t network_send;
  static size_t local_send;
  static double bytes_sent;
  static size_t ppn;
#if defined(HAVE_CUDA)
  struct CudaContext {
    cublasStatus_t status;
    cublasHandle_t handle;
  };
  static CudaContext cuda;
  static struct KernelDimensions {
    struct {
      size_t blocks, threads;
    } ooo;
  } kernel_dimensions;
#endif

  static void init(MPI_Comm);

  template <typename F = double>
  struct Input {
    CTF::Tensor<F> *ei = nullptr, *ea = nullptr, *Tph = nullptr,
                   *Tpphh = nullptr, *Vpphh = nullptr, *Vhhhp = nullptr,
                   *Vppph = nullptr, *Jppph = nullptr, *Jhhhp = nullptr;
    Input &with_epsilon_i(CTF::Tensor<F> *t) {
      ei = t;
      return *this;
    }
    Input &with_epsilon_a(CTF::Tensor<F> *t) {
      ea = t;
      return *this;
    }
    Input &with_Tai(CTF::Tensor<F> *t) {
      Tph = t;
      return *this;
    }
    Input &with_Tabij(CTF::Tensor<F> *t) {
      Tpphh = t;
      return *this;
    }
    Input &with_Vabij(CTF::Tensor<F> *t) {
      Vpphh = t;
      return *this;
    }
    Input &with_Vijka(CTF::Tensor<F> *t) {
      Vhhhp = t;
      return *this;
    }
    Input &with_Vabci(CTF::Tensor<F> *t) {
      Vppph = t;
      return *this;
    }
    Input &with_Jijka(CTF::Tensor<F> *t) {
      Jhhhp = t;
      return *this;
    }
    Input &with_Jabci(CTF::Tensor<F> *t) {
      Jppph = t;
      return *this;
    }
    enum TuplesDistribution {
      NAIVE,
      GROUP_AND_SORT,
    };

    ADD_ATTRIBUTE(bool, delete_Vppph, false)
    ADD_ATTRIBUTE(bool, rank_round_robin, false)
    ADD_ATTRIBUTE(bool, chrono, false)
    ADD_ATTRIBUTE(bool, barrier, false)
    ADD_ATTRIBUTE(bool, blocking, false)
    ADD_ATTRIBUTE(size_t, max_iterations, 0)
    ADD_ATTRIBUTE(int, iteration_mod, -1)
    ADD_ATTRIBUTE(int, percentage_mod, -1)
    ADD_ATTRIBUTE(TuplesDistribution, tuples_distribution, NAIVE)
    ADD_ATTRIBUTE(std::string, checkpoint_path, "atrip-checkpoint.yaml")
    ADD_ATTRIBUTE(bool, read_checkpoint_if_exists, true)
    ADD_ATTRIBUTE(bool, writeCheckpoint, true)
    ADD_ATTRIBUTE(float, checkpoint_at_percentage, 10)
    ADD_ATTRIBUTE(size_t, checkpoint_at_every_iteration, 0)
    ADD_ATTRIBUTE(bool, ijkabc, 0)
#if defined(HAVE_CUDA)
    ADD_ATTRIBUTE(size_t, ooo_threads, 0)
    ADD_ATTRIBUTE(size_t, ooo_blocks, 0)
#endif
  };

  struct Output {
    double energy;
    double ct_energy;
  };
  template <typename F = double>
  static Output run(Input<F> const &in);
};

} // namespace atrip

#undef ADD_ATTRIBUTE
// Header:1 ends here

#endif
