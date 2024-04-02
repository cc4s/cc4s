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

#include <vector>
#include <array>
#include <numeric>

// TODO: remove some
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <map>
#include <cassert>
#include <chrono>
#include <climits>

#include <atrip/mpi.hpp>
#include <atrip/Utils.hpp>
#include <atrip/Debug.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Tuples%20types][Tuples types:1]]
using ABCTuple = std::array<size_t, 3>;
using PartialTuple = std::array<size_t, 2>;
using ABCTuples = std::vector<ABCTuple>;

constexpr ABCTuple FAKE_TUPLE = {0, 0, 0};
constexpr ABCTuple INVALID_TUPLE = {1, 1, 1};
// Tuples types:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Distributing%20the%20tuples][Distributing the
// tuples:1]]
struct TuplesDistribution {
  virtual ABCTuples get_tuples(size_t Nv, MPI_Comm universe) = 0;
  virtual bool tuple_is_fake(ABCTuple const &t) { return t == FAKE_TUPLE; }
};
// Distributing the tuples:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Node%20information][Node information:1]]
std::vector<std::string> get_node_names(MPI_Comm comm);
// Node information:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Node%20information][Node information:2]]
struct RankInfo {
  const std::string name;
  const size_t node_id;
  const size_t global_rank;
  const size_t local_rank;
  const size_t ranks_per_node;
};

std::vector<RankInfo> get_node_infos(std::vector<string> const &node_names);

struct ClusterInfo {
  const size_t n_nodes, np, ranks_per_node;
  const std::vector<RankInfo> rank_infos;
};

ClusterInfo get_cluster_info(MPI_Comm comm);
// Node information:2 ends here

// [[file:~/cuda/atrip/atrip.org::*Naive%20list][Naive list:1]]
ABCTuples get_tuples_list(size_t Nv, size_t rank, size_t np);
// Naive list:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Naive%20list][Naive list:2]]
ABCTuples get_all_tuples_list(const size_t Nv);
// Naive list:2 ends here

// [[file:~/cuda/atrip/atrip.org::*Naive%20list][Naive list:3]]
struct NaiveDistribution : public TuplesDistribution {
  ABCTuples get_tuples(size_t Nv, MPI_Comm universe) override;
};
// Naive list:3 ends here

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:1]]
namespace group_and_sort {
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Utils][Utils:1]]
// Provides the node on which the slice-element is found
// Right now we distribute the slices in a round robin fashion
// over the different nodes (NOTE: not mpi ranks but nodes)
inline size_t is_on_node(size_t tuple, size_t n_nodes);

// return the node (or all nodes) where the elements of this
// tuple are located
std::vector<size_t> get_tuple_nodes(ABCTuple const &t, size_t n_nodes);

struct Info {
  size_t n_nodes;
  size_t node_id;
};
// Utils:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Distribution][Distribution:1]]
ABCTuples special_distribution(Info const &info, ABCTuples const &all_tuples);
// Distribution:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Main][Main:1]]
std::vector<ABCTuple> main(MPI_Comm universe, size_t Nv);
// Main:5 ends here

// [[file:~/cuda/atrip/atrip.org::*Interface][Interface:1]]
struct Distribution : public TuplesDistribution {
  ABCTuples get_tuples(size_t Nv, MPI_Comm universe) override;
};
// Interface:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:1]]
} // namespace group_and_sort
// Epilog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:1]]
} // namespace atrip
// Epilog:1 ends here
