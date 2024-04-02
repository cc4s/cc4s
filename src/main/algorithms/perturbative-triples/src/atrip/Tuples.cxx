#include <atrip/Tuples.hpp>
#include <atrip/Atrip.hpp>

namespace atrip {

template <typename A>
static A unique(A const &xs) {
  auto result = xs;
  std::sort(std::begin(result), std::end(result));
  auto const &last = std::unique(std::begin(result), std::end(result));
  result.erase(last, std::end(result));
  return result;
}

std::vector<std::string> get_node_names(MPI_Comm comm) {
  int rank, np;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &np);

  std::vector<std::string> node_list(np);
  char node_name[MPI_MAX_PROCESSOR_NAME];
  char *node_names = (char *)malloc(np * MPI_MAX_PROCESSOR_NAME);
  std::vector<int> name_lengths(np), off(np);
  int name_length;
  MPI_Get_processor_name(node_name, &name_length);
  MPI_Allgather(&name_length,
                1,
                MPI_INT,
                name_lengths.data(),
                1,
                MPI_INT,
                comm);
  for (int i(1); i < np; i++) off[i] = off[i - 1] + name_lengths[i - 1];
  MPI_Allgatherv(node_name,
                 name_lengths[rank],
                 MPI_BYTE,
                 node_names,
                 name_lengths.data(),
                 off.data(),
                 MPI_BYTE,
                 comm);
  for (int i(0); i < np; i++) {
    std::string const s(&node_names[off[i]], name_lengths[i]);
    node_list[i] = s;
  }
  std::free(node_names);
  return node_list;
}

std::vector<RankInfo> get_node_infos(std::vector<string> const &node_names) {
  std::vector<RankInfo> result;
  std::vector<std::string> unique_names;
  for (auto const &name : node_names) {
    if (unique_names.size() == 0
        || unique_names.end()
               == std::find(unique_names.begin(), unique_names.end(), name)) {
      unique_names.push_back(name);
    }
  }
  // auto const unique_names = unique(node_names);
  auto const index = [&unique_names](std::string const &s) {
    auto const &it = std::find(unique_names.begin(), unique_names.end(), s);
    return std::distance(unique_names.begin(), it);
  };
  std::vector<size_t> local_ranks(unique_names.size(), 0);
  size_t global_rank = 0;
  for (auto const &name : node_names) {
    const size_t node_id = index(name);
    result.push_back(
        {name,
         node_id,
         global_rank++,
         local_ranks[node_id]++,
         (size_t)std::count(node_names.begin(), node_names.end(), name)});
  }
  return result;
}

ClusterInfo get_cluster_info(MPI_Comm comm) {
  auto const names = get_node_names(comm);
  auto const rank_infos = get_node_infos(names);

  return ClusterInfo{unique(names).size(),
                     names.size(),
                     rank_infos[0].ranks_per_node,
                     rank_infos};
}

ABCTuples get_tuples_list(size_t Nv, size_t rank, size_t np) {

  const size_t
      // total number of tuples for the problem
      n = Nv * (Nv + 1) * (Nv + 2) / 6 - Nv

      // all ranks should have the same number of tuples_per_rank
      ,
      tuples_per_rank = n / np + size_t(n % np != 0)

      // start index for the global tuples list
      ,
      start = tuples_per_rank * rank

      // end index for the global tuples list
      ,
      end = tuples_per_rank * (rank + 1);

  LOG(1, "Atrip") << "tuples_per_rank = " << tuples_per_rank << "\n";
  WITH_RANK << "start, end = " << start << ", " << end << "\n";
  ABCTuples result(tuples_per_rank, FAKE_TUPLE);

  for (size_t a(0), r(0), g(0); a < Nv; a++)
    for (size_t b(a); b < Nv; b++)
      for (size_t c(b); c < Nv; c++) {
        if (a == b && b == c) continue;
        if (start <= g && g < end) result[r++] = {a, b, c};
        g++;
      }

  return result;
}

ABCTuples get_all_tuples_list(const size_t Nv) {
  const size_t n = Nv * (Nv + 1) * (Nv + 2) / 6 - Nv;
  ABCTuples result(n);

  for (size_t a(0), u(0); a < Nv; a++)
    for (size_t b(a); b < Nv; b++)
      for (size_t c(b); c < Nv; c++) {
        if (a == b && b == c) continue;
        result[u++] = {a, b, c};
      }

  return result;
}

ABCTuples atrip::NaiveDistribution::get_tuples(size_t Nv, MPI_Comm universe) {
  int rank, np;
  MPI_Comm_rank(universe, &rank);
  MPI_Comm_size(universe, &np);
  return get_tuples_list(Nv, (size_t)rank, (size_t)np);
}

namespace group_and_sort {

inline size_t is_on_node(size_t tuple, size_t n_nodes) {
  return tuple % n_nodes;
}

std::vector<size_t> get_tuple_nodes(ABCTuple const &t, size_t n_nodes) {
  std::vector<size_t> n_tuple = {is_on_node(t[0], n_nodes),
                                 is_on_node(t[1], n_nodes),
                                 is_on_node(t[2], n_nodes)};
  return unique(n_tuple);
}

ABCTuples special_distribution(Info const &info, ABCTuples const &all_tuples) {

  ABCTuples node_tuples;
  size_t const n_nodes(info.n_nodes);

  std::vector<ABCTuples> container1d(n_nodes), container2d(n_nodes * n_nodes),
      container3d(n_nodes * n_nodes * n_nodes);

  WITH_DBG if (info.node_id == 0) std::cout
      << "\tGoing through all " << all_tuples.size() << " tuples in " << n_nodes
      << " nodes\n";

  // build container-n-d's
  for (auto const &t : all_tuples) {
    // one which node(s) are the tuple elements located...
    // put them into the right container
    auto const _nodes = get_tuple_nodes(t, n_nodes);

    switch (_nodes.size()) {
    case 1: container1d[_nodes[0]].push_back(t); break;
    case 2: container2d[_nodes[0] + _nodes[1] * n_nodes].push_back(t); break;
    case 3:
      container3d[_nodes[0] + _nodes[1] * n_nodes
                  + _nodes[2] * n_nodes * n_nodes]
          .push_back(t);
      break;
    }
  }

  WITH_DBG if (info.node_id == 0) std::cout << "\tBuilding 1-d containers\n";
  // DISTRIBUTE 1-d containers
  // every tuple which is only located at one node belongs to this node
  {
    auto const &_tuples = container1d[info.node_id];
    node_tuples.resize(_tuples.size(), INVALID_TUPLE);
    std::copy(_tuples.begin(), _tuples.end(), node_tuples.begin());
  }

  WITH_DBG if (info.node_id == 0) std::cout << "\tBuilding 2-d containers\n";
  // DISTRIBUTE 2-d containers
  // the tuples which are located at two nodes are half/half given to these
  // nodes
  for (size_t yx = 0; yx < container2d.size(); yx++) {

    auto const &_tuples = container2d[yx];
    const size_t idx = yx % n_nodes
        // remeber: yx = idy * n_nodes + idx
        ,
                 idy = yx / n_nodes, n_half = _tuples.size() / 2,
                 size = node_tuples.size();

    size_t nbeg, nend;
    if (info.node_id == idx) {
      nbeg = 0 * n_half;
      nend = n_half;
    } else if (info.node_id == idy) {
      nbeg = 1 * n_half;
      nend = _tuples.size();
    } else {
      // either idx or idy is my node
      continue;
    }

    size_t const nextra = nend - nbeg;
    node_tuples.resize(size + nextra, INVALID_TUPLE);
    std::copy(_tuples.begin() + nbeg,
              _tuples.begin() + nend,
              node_tuples.begin() + size);
  }

  WITH_DBG if (info.node_id == 0) std::cout << "\tBuilding 3-d containers\n";
  // DISTRIBUTE 3-d containers
  for (size_t zyx = 0; zyx < container3d.size(); zyx++) {
    auto const &_tuples = container3d[zyx];

    const size_t idx = zyx % n_nodes,
                 idy = (zyx / n_nodes) % n_nodes
        // remember: zyx = idx + idy * n_nodes + idz * n_nodes^2
        ,
                 idz = zyx / n_nodes / n_nodes, n_third = _tuples.size() / 3,
                 size = node_tuples.size();

    size_t nbeg, nend;
    if (info.node_id == idx) {
      nbeg = 0 * n_third;
      nend = 1 * n_third;
    } else if (info.node_id == idy) {
      nbeg = 1 * n_third;
      nend = 2 * n_third;
    } else if (info.node_id == idz) {
      nbeg = 2 * n_third;
      nend = _tuples.size();
    } else {
      // either idx or idy or idz is my node
      continue;
    }

    size_t const nextra = nend - nbeg;
    node_tuples.resize(size + nextra, INVALID_TUPLE);
    std::copy(_tuples.begin() + nbeg,
              _tuples.begin() + nend,
              node_tuples.begin() + size);
  }

  WITH_DBG if (info.node_id == 0) std::cout << "\tswapping tuples...\n";
  /*
   *  sort part of group-and-sort algorithm
   *  every tuple on a given node is sorted in a way that
   *  the 'home elements' are the fastest index.
   *  1:yyy 2:yyn(x) 3:yny(x) 4:ynn(x) 5:nyy 6:nyn(x) 7:nny 8:nnn
   */
  for (auto &nt : node_tuples) {
    if (is_on_node(nt[0], n_nodes) == info.node_id) {   // 1234
      if (is_on_node(nt[2], n_nodes) != info.node_id) { // 24
        size_t const x(nt[0]);
        nt[0] = nt[2]; // switch first and last
        nt[2] = x;
      } else if (is_on_node(nt[1], n_nodes) != info.node_id) { // 3
        size_t const x(nt[0]);
        nt[0] = nt[1]; // switch first two
        nt[1] = x;
      }
    } else {
      if (is_on_node(nt[1], n_nodes) == info.node_id       // 56
          && is_on_node(nt[2], n_nodes) != info.node_id) { // 6
        size_t const x(nt[1]);
        nt[1] = nt[2]; // switch last two
        nt[2] = x;
      }
    }
  }

  WITH_DBG if (info.node_id == 0) std::cout << "\tsorting list of tuples...\n";
  // now we sort the list of tuples
  std::sort(node_tuples.begin(), node_tuples.end());

  WITH_DBG if (info.node_id == 0) std::cout << "\trestoring tuples...\n";
  // we bring the tuples abc back in the order a<b<c
  for (auto &t : node_tuples) std::sort(t.begin(), t.end());

#if ATRIP_DEBUG > 1
  WITH_DBG if (info.node_id == 0) std::cout << "checking for validity of "
                                            << node_tuples.size() << std::endl;
  const bool any_invalid =
      std::any_of(node_tuples.begin(),
                  node_tuples.end(),
                  [](ABCTuple const &t) { return t == INVALID_TUPLE; });
  if (any_invalid) throw "Some tuple is invalid in group-and-sort algorithm";
#endif

  WITH_DBG if (info.node_id == 0) std::cout << "\treturning tuples...\n";
  return node_tuples;
}

std::vector<ABCTuple> main(MPI_Comm universe, size_t Nv) {

  int rank, np;
  MPI_Comm_rank(universe, &rank);
  MPI_Comm_size(universe, &np);

  std::vector<ABCTuple> result;

  // auto const node_names(get_node_names(universe));
  // size_t const n_nodes = unique(node_names).size();
  // auto const node_infos = get_node_infos(node_names);
  auto cluster_info = Atrip::cluster_info;
  auto const node_infos = cluster_info->rank_infos;
  size_t const n_nodes = cluster_info->n_nodes;

  // We want to construct a communicator which only contains of one
  // element per node
  bool const compute_distribution_p = node_infos[rank].local_rank == 0;

  std::vector<ABCTuple> node_tuples =
      compute_distribution_p
          ? special_distribution(Info{n_nodes, node_infos[rank].node_id},
                                 get_all_tuples_list(Nv))
          : std::vector<ABCTuple>();

  LOG(1, "Atrip") << "got node_tuples\n";

  // now we have to send the data from **one** rank on each node
  // to all others ranks of this node
  const int color = node_infos[rank].node_id, key = node_infos[rank].local_rank;

  MPI_Comm INTRA_COMM;
  MPI_Comm_split(universe, color, key, &INTRA_COMM);
  // Main:1 ends here

  // [[file:~/cuda/atrip/atrip.org::*Main][Main:2]]
  size_t const tuples_per_rank_local =
      node_tuples.size() / node_infos[rank].ranks_per_node
      + size_t(node_tuples.size() % node_infos[rank].ranks_per_node != 0);

  size_t tuples_per_rank_global;

  MPI_Reduce(&tuples_per_rank_local,
             &tuples_per_rank_global,
             1,
             MPI_UINT64_T,
             MPI_MAX,
             0,
             universe);

  MPI_Bcast(&tuples_per_rank_global, 1, MPI_UINT64_T, 0, universe);

  LOG(1, "Atrip") << "Tuples per rank: " << tuples_per_rank_global << "\n";
  LOG(1, "Atrip") << "ranks per node " << node_infos[rank].ranks_per_node
                  << "\n";
  LOG(1, "Atrip") << "#nodes " << n_nodes << "\n";
  // Main:2 ends here

  // [[file:~/cuda/atrip/atrip.org::*Main][Main:3]]
  size_t const total_tuples =
      tuples_per_rank_global * node_infos[rank].ranks_per_node;

  if (compute_distribution_p) {
    // pad with FAKE_TUPLEs
    node_tuples.insert(node_tuples.end(),
                       total_tuples - node_tuples.size(),
                       FAKE_TUPLE);
  }
  // Main:3 ends here

  // [[file:~/cuda/atrip/atrip.org::*Main][Main:4]]
  {
    // construct mpi type for abctuple
    MPI_Datatype MPI_ABCTUPLE;
    MPI_Type_vector(node_tuples[0].size(), 1, 1, MPI_UINT64_T, &MPI_ABCTUPLE);
    MPI_Type_commit(&MPI_ABCTUPLE);

    LOG(1, "Atrip") << "scattering tuples \n";

    result.resize(tuples_per_rank_global);
    MPI_Scatter(node_tuples.data(),
                tuples_per_rank_global,
                MPI_ABCTUPLE,
                result.data(),
                tuples_per_rank_global,
                MPI_ABCTUPLE,
                0,
                INTRA_COMM);

    MPI_Type_free(&MPI_ABCTUPLE);
  }

  return result;
}

ABCTuples Distribution::get_tuples(size_t Nv, MPI_Comm universe) {
  return main(universe, Nv);
}

} // namespace group_and_sort
} // namespace atrip
