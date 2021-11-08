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

#ifndef MPI_COMMUNICATOR_DEFINED
#define MPI_COMMUNICATOR_DEFINED

#include "mpi.h"
#include <math/Complex.hpp>
#include <math/Vector.hpp>

namespace cc4s {
  // base template for type traits
  template <typename F>
  class MpiTypeTraits;


  class MpiCommunicator {
  public:
    MpiCommunicator(
      MPI_Comm comm_ = MPI_COMM_WORLD
    ): comm(comm_) {
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &processes);
    }
    MpiCommunicator(
      int rank_, int processes_, MPI_Comm comm_ = MPI_COMM_WORLD
    ): rank(rank_), processes(processes_), comm(comm_) {
    }
    ~MpiCommunicator() {
    }

    void barrier() {
      MPI_Barrier(comm);
    }

    template <typename F>
    void reduce(const F &src, F &dst, int rootRank = 0) {
      MPI_Reduce(
        &src, &dst,
        MpiTypeTraits<F>::elementCount(), MpiTypeTraits<F>::elementType(),
        MPI_SUM, rootRank, comm
      );
    }

    template <typename F>
    void allReduce(const F &src, F &dst) {
      MPI_Allreduce(
        &src, &dst,
        MpiTypeTraits<F>::elementCount(), MpiTypeTraits<F>::elementType(),
        MPI_SUM, comm
      );
    }

    /**
     * \Brief Gathers the src vectors of all ranks together to the dst
     * vector at the given root rank, by default rank 0.
     **/
    template <typename F>
    void gather(
      const std::vector<F> &src, std::vector<F> &dst, int rootRank = 0
    ) {
      if (rank == rootRank) {
        dst.resize(src.size() * processes);
      } else {
        dst.resize(0);
      }
      MPI_Gather(
        src.data(), src.size() * MpiTypeTraits<F>::elementCount(),
        MpiTypeTraits<F>::elementType(),
        dst.data(), src.size() * MpiTypeTraits<F>::elementCount(),
        MpiTypeTraits<F>::elementType(),
        rootRank, comm
      );
    }

    template <typename F>
    void broadcast(
      std::vector<F> &dst, int rootRank = 0
    ) {
      MPI_Bcast(
        dst.data(), dst.size() * MpiTypeTraits<F>::elementCount(),
        MpiTypeTraits<F>::elementType(),
        rootRank,
        comm
      );
    }

    int getRank() const {
      return rank;
    }
    int getProcesses() const {
      return processes;
    }
    MPI_Comm getComm() const {
      return comm;
    }

  protected:
    int rank, processes;
    MPI_Comm comm;
  };


  template <>
  class MpiTypeTraits<int> {
  public:
    static MPI_Datatype elementType() { return MPI_INT; }
    static int elementCount()  { return 1; }
  };

  template <>
  class MpiTypeTraits<int64_t> {
  public:
    static MPI_Datatype elementType() { return MPI_INTEGER8; }
    static int elementCount() { return 1; }
  };

  template <>
  class MpiTypeTraits<uint64_t> {
  public:
    static MPI_Datatype elementType() { return MPI_INTEGER8; }
    static int elementCount() { return 1; }
  };

  template <>
  class MpiTypeTraits<Real<64>> {
  public:
    static MPI_Datatype elementType() { return MPI_REAL8; }
    static int elementCount() { return 1; }
  };

  template <>
  class MpiTypeTraits<Complex<64>> {
  public:
    static MPI_Datatype elementType() { return MPI_DOUBLE_COMPLEX; }
    static int elementCount() { return 1; }
  };

  template <typename F, int D>
  class MpiTypeTraits<Vector<F, D>> {
  public:
    static MPI_Datatype elementType() {
      return MpiTypeTraits<F>::elementType();
    }
    static int elementCount() { return D; }
  };

  // TODO: 128 bit reals and floats
}

#endif

