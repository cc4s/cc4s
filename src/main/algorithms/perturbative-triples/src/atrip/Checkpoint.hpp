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
#include <fstream>
#include <iomanip>

#include <atrip/Atrip.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::checkpoint-definition][checkpoint-definition]]
// template <typename F>
struct Checkpoint {
  size_t no, nv;
  size_t nranks;
  size_t nnodes;
  double global_energy;
  double iteration_energy;
  double global_ct_energy;
  double iteration_ct_energy;
  size_t iteration;
  // TODO
  // Input<F>::TuplesDistribution distribution(GROUP_AND_SORT);
  bool rank_round_robin;
};
// checkpoint-definition ends here

// [[file:~/cuda/atrip/atrip.org::*Input%20and%20output][Input and output:1]]

void write_checkpoint_header(std::string const &filepath);

void write_checkpoint(Checkpoint const &c, std::string const &filepath);

Checkpoint read_checkpoint(std::ifstream &in);

Checkpoint read_checkpoint(std::string const &filepath);
// Input and output:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:1]]
} // namespace atrip
// Epilog:1 ends here
