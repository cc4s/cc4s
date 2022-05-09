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

// [[file:../../atrip.org::*The rank mapping][The rank mapping:1]]
#pragma once

#include <vector>
#include <algorithm>

#include <atrip/Slice.hpp>
#include <atrip/Tuples.hpp>

namespace atrip {

  template <typename F=double>
  struct RankMap {

    static bool RANK_ROUND_ROBIN;
    std::vector<size_t> const lengths;
    size_t const np, size;
    ClusterInfo const clusterInfo;

    RankMap(std::vector<size_t> lens, size_t np_, MPI_Comm comm)
      : lengths(lens)
      , np(np_)
      , size(std::accumulate(lengths.begin(), lengths.end(),
                            1UL, std::multiplies<size_t>()))
      , clusterInfo(getClusterInfo(comm))
    { assert(lengths.size() <= 2); }

    size_t find(typename Slice<F>::Location const& p) const noexcept {
      if (RANK_ROUND_ROBIN) {
        return p.source * np + p.rank;
      } else {
        const size_t
          rankPosition = p.source * clusterInfo.ranksPerNode
                       + clusterInfo.rankInfos[p.rank].localRank
                       ;
        return rankPosition * clusterInfo.nNodes
             + clusterInfo.rankInfos[p.rank].nodeId
             ;
      }
    }

    size_t nSources() const noexcept {
      return size / np + size_t(size % np != 0);
    }


    bool isPaddingRank(size_t rank) const noexcept {
      return size % np == 0
          ? false
          : rank > (size % np - 1)
          ;
    }

    bool isSourcePadding(size_t rank, size_t source) const noexcept {
      return source == nSources() && isPaddingRank(rank);
    }

    typename Slice<F>::Location
    find(ABCTuple const& abc, typename Slice<F>::Type sliceType) const {
      // tuple = {11, 8} when abc = {11, 8, 9} and sliceType = AB
      // tuple = {11, 0} when abc = {11, 8, 9} and sliceType = A
      const auto tuple = Slice<F>::subtupleBySlice(abc, sliceType);

      const size_t index
        = tuple[0]
        + tuple[1] * (lengths.size() > 1 ? lengths[0] : 0)
        ;

      size_t rank, source;

      if (RANK_ROUND_ROBIN) {

        rank = index % np;
        source = index / np;

      } else {

        size_t const

          // the node that will be assigned to
            nodeId = index % clusterInfo.nNodes

          // how many times it has been assigned to the node
          , s_n = index / clusterInfo.nNodes

          // which local rank in the node should be
          , localRank = s_n % clusterInfo.ranksPerNode

          // and the local source (how many times we chose this local rank)
          , localSource = s_n / clusterInfo.ranksPerNode
          ;

        // find the localRank-th entry in clusterInfo
        auto const& it =
          std::find_if(clusterInfo.rankInfos.begin(),
                       clusterInfo.rankInfos.end(),
                       [nodeId, localRank](RankInfo const& ri) {
                         return ri.nodeId == nodeId
                             && ri.localRank == localRank
                             ;
                       });
        if (it == clusterInfo.rankInfos.end()) {
          throw "FATAL! Error in node distribution of the slices";
        }

        rank = (*it).globalRank;
        source = localSource;

      }

      return
        { rank
        , source
        };
    }

  };

}
// The rank mapping:1 ends here
