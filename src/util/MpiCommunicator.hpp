#ifndef MPI_COMMUNICATOR_DEFINED
#define MPI_COMMUNICATOR_DEFINED

#include "mpi.h"
#include <math/Complex.hpp>
#include <math/Vector.hpp>
#include <ctf.hpp>

namespace cc4s {
  // base template for type traits
  template <typename F>
  class MpiTypeTraits;


  class MpiCommunicator {
  public:
    MpiCommunicator(
      int rank_, int processes_, MPI_Comm comm_ = MPI_COMM_WORLD
    ): rank(rank_), processes(processes_), comm(comm_) {
    }
    MpiCommunicator(
      const CTF::World &world
    ): rank(world.rank), processes(world.np), comm(world.comm) {     
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
        MpiTypeTraits<F>::ElementCount, MpiTypeTraits<F>::ElementType,
        MPI_SUM, rootRank, comm
      );
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
    int rank, processes;
    MPI_Comm comm;
  };


  template <>
  class MpiTypeTraits<int> {
  public:
    static constexpr MPI_Datatype ElementType = MPI_INT;
    static constexpr int ElementCount = 1;
  };

  template <>
  class MpiTypeTraits<int64_t> {
  public:
    static constexpr MPI_Datatype ElementType = MPI_INTEGER8;
    static constexpr int ElementCount = 1;
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

  template <typename F, int D>
  class MpiTypeTraits<Vector<F, D>> {
  public:
    static constexpr MPI_Datatype ElementType = MpiTypeTraits<F>::ElementType;
    static constexpr int ElementCount = D;
  };
}

#endif

