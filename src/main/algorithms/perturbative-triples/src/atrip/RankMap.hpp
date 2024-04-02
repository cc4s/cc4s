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

// [[file:~/cuda/atrip/atrip.org::*The%20rank%20mapping][The rank mapping:1]]
#pragma once

#include <vector>
#include <algorithm>

#include <atrip/Slice.hpp>
#include <atrip/Tuples.hpp>

namespace atrip {

template <typename F = double>
class RankMap {
public:
  static bool RANK_ROUND_ROBIN;
  std::vector<size_t> const lengths;
  size_t const np, size;
  ClusterInfo const cluster_info;

  RankMap(std::vector<size_t> lens, size_t np_);

  size_t find(typename Slice<F>::Location const &p) const noexcept;

  size_t n_sources() const noexcept;

  bool is_padding_rank(size_t rank) const noexcept;

  bool is_source_padding(const size_t rank, const size_t source) const noexcept;

  typename Slice<F>::Location find(ABCTuple const &abc,
                                   typename Slice<F>::Type slice_type) const;
};

} // namespace atrip
// The rank mapping:1 ends here
