#ifndef MPI_COMMUNICATOR_DEFINED
#define MPI_COMMUNICATOR_DEFINED

#include "mpi.h"
#include <math/Complex.hpp>

namespace cc4s {
  // base template for type traits
  template <typename F>
  class MpiTypeTraits;


  class MpiCommunicator {
  public:
    MpiCommunicator(
      int rank_, int processes_, int comm_ = MPI_COMM_WORLD
    ): rank(rank_), processes(processes_), comm(comm_) {
    }
    ~MpiCommunicator() {
    }

    void barrier() {
      MPI_Barrier(comm);
    }

    template <typename F>
    void allReduce(const F &src, F &dst) {
      MPI_Allreduce(
        &src, &dst,
        MpiTypeTraits<F>::ElementCount, MpiTypeTraits<F>::ElementType,
        MPI_SUM, comm
      );
    }

    int getRank() const {
      return rank;
    }
    int getProcesses() const {
      return processes;
    }

  protected:
    int rank, processes, comm;
  };


  template <>
  class MpiTypeTraits<uint64_t> {
  public:
    static constexpr MPI_Datatype ElementType = MPI_INTEGER8;
    static constexpr int ElementCount = 1;
  };

  template <>
  class MpiTypeTraits<double> {
  public:
    static constexpr MPI_Datatype ElementType = MPI_REAL8;
    static constexpr int ElementCount = 1;
  };

  template <>
  class MpiTypeTraits<complex> {
  public:
    static constexpr MPI_Datatype ElementType = MPI_DOUBLE_COMPLEX;
    static constexpr int ElementCount = 1;
  };
}

#endif

