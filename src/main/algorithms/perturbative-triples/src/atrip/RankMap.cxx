#include <atrip/RankMap.hpp>
#include <atrip/Atrip.hpp>
#include <atrip/Types.hpp>

namespace atrip {

template <typename F>
size_t RankMap<F>::find(typename Slice<F>::Location const &p) const noexcept {
  if (RANK_ROUND_ROBIN) {
    return p.source * np + p.rank;
  } else {
    const size_t rank_position = p.source * cluster_info.ranks_per_node
                               + cluster_info.rank_infos[p.rank].local_rank;
    return rank_position * cluster_info.n_nodes
         + cluster_info.rank_infos[p.rank].node_id;
  }
}

template <typename F>
size_t RankMap<F>::n_sources() const noexcept {
  return size / np + size_t(size % np != 0);
}

template <typename F>
bool RankMap<F>::is_padding_rank(size_t rank) const noexcept {
  return size % np == 0 ? false : rank > (size % np - 1);
}

template <typename F>
bool RankMap<F>::is_source_padding(const size_t rank,
                                   const size_t source) const noexcept {
  return source == n_sources() && is_padding_rank(rank);
}

template <typename F>
typename Slice<F>::Location
RankMap<F>::find(ABCTuple const &abc,
                 typename Slice<F>::Type slice_type) const {
  // tuple = {11, 8} when abc = {11, 8, 9} and slice_type = AB
  // tuple = {11, 0} when abc = {11, 8, 9} and slice_type = A
  const auto tuple = Slice<F>::subtuple_by_slice(abc, slice_type);

  const size_t index =
      tuple[0] + tuple[1] * (lengths.size() > 1 ? lengths[0] : 0);

  size_t rank, source;

  if (RANK_ROUND_ROBIN) {

    rank = index % np;
    source = index / np;

  } else {

    size_t const

        // the node that will be assigned to
        node_id = index % cluster_info.n_nodes,

        // how many times it has been assigned to the node
        s_n = index / cluster_info.n_nodes,

        // which local rank in the node should be
        local_rank = s_n % cluster_info.ranks_per_node,

        // and the local source (how many times we chose this local rank)
        local_source = s_n / cluster_info.ranks_per_node;

    // find the local_rank-th entry in cluster_info
    auto const &it = std::find_if(cluster_info.rank_infos.begin(),
                                  cluster_info.rank_infos.end(),
                                  [node_id, local_rank](RankInfo const &ri) {
                                    return ri.node_id == node_id
                                        && ri.local_rank == local_rank;
                                  });
    if (it == cluster_info.rank_infos.end()) {
      throw "FATAL! Error in node distribution of the slices";
    }

    rank = (*it).global_rank;
    source = local_source;
  }

  return {rank, source};
}

template <typename F>
RankMap<F>::RankMap(std::vector<size_t> lens, size_t np_)
    : lengths(lens)
    , np(np_)
    , size(std::accumulate(lengths.begin(),
                           lengths.end(),
                           1UL,
                           std::multiplies<size_t>()))
    , cluster_info(*Atrip::cluster_info) {
  assert(lengths.size() <= 2);
}

template class RankMap<Complex>;
template class RankMap<double>;

} // namespace atrip
