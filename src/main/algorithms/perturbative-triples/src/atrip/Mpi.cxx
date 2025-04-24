#include <atrip/Mpi.hpp>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>

namespace atrip {

MPI_Comm create_comm(MPI_Comm in_comm, int omp_stride, bool use_round_robin) {
  int glb_rank, glb_size;
  MPI_Comm_rank(in_comm, &glb_rank);
  MPI_Comm_size(in_comm, &glb_size);

  // --- Get hostname
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(hostname, &name_len);
  std::string host_str(hostname, name_len);

  // --- Gather hostnames and ranks
  struct Meta {
    int glb_rank;
    char hostname[MPI_MAX_PROCESSOR_NAME];
  };

  Meta local_meta;
  local_meta.glb_rank = glb_rank;
  std::strncpy(local_meta.hostname, hostname, MPI_MAX_PROCESSOR_NAME);

  std::vector<Meta> all(glb_size);
  MPI_Allgather(&local_meta, sizeof(Meta), MPI_BYTE, all.data(), sizeof(Meta), MPI_BYTE, in_comm);

  // --- Reconstruct node-local maps (preserve insertion order)
  std::vector<std::string> host_order;
  std::unordered_map<std::string, std::vector<int>> node_map;

  for (const auto& meta : all) {
    std::string hname(meta.hostname);
    if (node_map.count(hname) == 0)
      host_order.push_back(hname);
    node_map[hname].push_back(meta.glb_rank);
  }

  // --- Apply omp_stride and prepare per-node filtered ranks
  std::vector<std::vector<int>> per_node_filtered;
  for (const auto& hname : host_order) {
    const auto& ranks = node_map[hname];
    std::vector<int> filtered;
    for (size_t i = 0; i < ranks.size(); ++i) {
      if (omp_stride == 1 || i % omp_stride == 0)
        filtered.push_back(ranks[i]);
    }
    per_node_filtered.push_back(filtered);
  }

  // --- Final flattening (round-robin or node-contiguous)
  std::vector<int> final_ranks;
  if (use_round_robin) {
    size_t max_len = 0;
    for (const auto& v : per_node_filtered)
      max_len = std::max(max_len, v.size());

    for (size_t i = 0; i < max_len; ++i) {
      for (const auto& v : per_node_filtered) {
        if (i < v.size())
          final_ranks.push_back(v[i]);
      }
    }
  } else {
    for (const auto& v : per_node_filtered) {
      final_ranks.insert(final_ranks.end(), v.begin(), v.end());
    }
  }

  // --- Check participation and assign color
  int color = MPI_UNDEFINED;
  int key = 0;
  for (size_t i = 0; i < final_ranks.size(); ++i) {
    if (final_ranks[i] == glb_rank) {
      color = 0;
      key = i;
      break;
    }
  }

  // --- Create communicator
  MPI_Comm out_comm;
  MPI_Comm_split(in_comm, color, key, &out_comm);
  return out_comm;
}

void print_comm_mapping(MPI_Comm in_comm, MPI_Comm final_comm, const std::string& label) {
  int glb_rank, glb_size;
  MPI_Comm_rank(in_comm, &glb_rank);
  MPI_Comm_size(in_comm, &glb_size);

  int final_rank = -1, final_size = -1;
  if (final_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(final_comm, &final_rank);
    MPI_Comm_size(final_comm, &final_size);
  }

  char hostname[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(hostname, &name_len);

  struct Entry {
    int grank, frank;
    char host[MPI_MAX_PROCESSOR_NAME];
  };

  Entry local{glb_rank, final_rank};
  std::memcpy(local.host, hostname, name_len + 1);

  std::vector<Entry> all(glb_size);
  MPI_Allgather(&local, sizeof(Entry), MPI_BYTE, all.data(), sizeof(Entry), MPI_BYTE, in_comm);

  if (glb_rank == 0) {
    std::cout << "\n=== Rank Mapping";
    if (!label.empty()) std::cout << " (" << label << ")";
    std::cout << " ===\n";
    std::cout << "Global  Final  Host\n";
    std::cout << "-----------------------------\n";
    for (const auto& e : all) {
      std::cout << std::setw(6) << e.grank << "  "
                << std::setw(5);
      if (e.frank >= 0)
        std::cout << e.frank;
      else
        std::cout << "-";
      std::cout << "  " << e.host << '\n';
    }
    std::cout << std::flush;
  }
}


std::vector<int> create_rank_lookup(MPI_Comm in_comm, MPI_Comm final_comm) {
  int glb_size, final_size;

  // Get size of the original communicator
  MPI_Comm_size(in_comm, &glb_size);

  // Get size of the final communicator
  MPI_Comm_size(final_comm, &final_size);

  // Precompute rank mapping for each final rank
  std::vector<int> rank_lookup(final_size, -1);
  std::vector<int> all_final_ranks(final_size);

  // Gather all final ranks in the final communicator
  MPI_Allgather(&glb_size, 1, MPI_INT, all_final_ranks.data(), 1, MPI_INT, final_comm);

  // Map final ranks to global ranks from the original communicator
  for (int i = 0; i < final_size; ++i) {
    int final_rank = all_final_ranks[i];
    rank_lookup[i] = final_rank; // Store the corresponding global rank
  }

  return rank_lookup;
}


} // namespace atrip
