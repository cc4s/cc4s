#ifndef MPI_HPP_ATRIP_
#define MPI_HPP_ATRIP_

#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wcast-function-type"
#endif

#include <mpi.h>
#include <string>
#include <vector>

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif

namespace atrip {

MPI_Comm create_comm(MPI_Comm in_comm, int omp_stide = 1, bool use_round_robin = true);

void print_comm_mapping(MPI_Comm in_comm, MPI_Comm final_comm, const std::string& label = "");

std::vector<int> create_rank_lookup(MPI_Comm in_comm, MPI_Comm final_comm);

} // end namespace trip

#endif
