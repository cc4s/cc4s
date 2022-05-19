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

// [[file:../../atrip.org::*Prolog][Prolog:1]]
#pragma once
#include <fstream>
#include <iomanip>

#include <atrip/Atrip.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:../../atrip.org::checkpoint-definition][checkpoint-definition]]
// template <typename F>
struct Checkpoint {
  size_t no, nv;
  size_t nranks;
  size_t nnodes;
  double energy;
  size_t iteration;
  // TODO
  // Input<F>::TuplesDistribution distribution(GROUP_AND_SORT);
  bool rankRoundRobin;
};
// checkpoint-definition ends here

// [[file:../../atrip.org::*Input and output][Input and output:1]]
void write_checkpoint(Checkpoint const& c, std::string const& filepath) {
  std::ofstream out(filepath);
  out << "No: " << c.no
      << "\n"
      << "Nv: " << c.nv
      << "\n"
      << "Nranks: " << c.nranks
      << "\n"
      << "Nnodes: " << c.nnodes
      << "\n"
      << "Energy: " << std::setprecision(19) << c.energy
      << "\n"
      << "Iteration: " << c.iteration
      << "\n"
      << "RankRoundRobin: " << (c.rankRoundRobin ? "true" : "false")
      << "\n";
}


Checkpoint read_checkpoint(std::ifstream& in) {
  Checkpoint c;
  // trim chars from the string, to be more sure and not use regexes
  auto trim = [](std::string& s, std::string const& chars) {
    s.erase(0, s.find_first_not_of(chars));
    s.erase(s.find_last_not_of(chars) + 1);
    return s;
  };
  for (std::string header, value; std::getline(in, header, ':');) {
    std::getline(in, value, '\n');
    trim(value, " \t"); // trim all whitespaces
    trim(header, " \t");

    /**/ if (header == "No")        c.no = std::atoi(value.c_str());
    else if (header == "Nv")        c.nv = std::atoi(value.c_str());
    else if (header == "Nranks")    c.nranks = std::atoi(value.c_str());
    else if (header == "Nnodes")    c.nnodes = std::atoi(value.c_str());
    else if (header == "Energy")    c.energy = std::atof(value.c_str());
    else if (header == "Iteration") c.iteration = std::atoi(value.c_str());
    else if (header == "RankRoundRobin") c.rankRoundRobin = (value[0] == 't');
  }
  return c;
}


Checkpoint read_checkpoint(std::string const& filepath) {
  std::ifstream in(filepath);
  return read_checkpoint(in);
}
// Input and output:1 ends here

// [[file:../../atrip.org::*Epilog][Epilog:1]]
}
// Epilog:1 ends here
