#include <atrip/RankMap.hpp>
#include <atrip/Atrip.hpp>
#include <atrip/Types.hpp>

namespace atrip {

template <typename F>
size_t RankMap<F>::find_element(typename Slice<F>::Location const &p) const noexcept {
  return p.source * np + p.rank;
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
RankMap<F>::find_location(ABCTuple const &abc,
                 typename Slice<F>::Type slice_type) const {
  // tuple = {11, 8} when abc = {11, 8, 9} and slice_type = AB
  // tuple = {11, 0} when abc = {11, 8, 9} and slice_type = A
  const auto tuple = Slice<F>::subtuple_by_slice(abc, slice_type);

  const size_t index =
      tuple[0] + tuple[1] * (lengths.size() > 1 ? lengths[0] : 0);

  int rank = index % Atrip::np,
      source = index / Atrip::np;
   return {rank, source};

}

template <typename F>
RankMap<F>::RankMap(std::vector<size_t> lens)
    : lengths(lens)
    , size(std::accumulate(lengths.begin(),
                           lengths.end(),
                           1UL,
                           std::multiplies<size_t>()))
    , np(Atrip::np) {
  assert(lengths.size() <= 2);
}



template class RankMap<Complex>;
template class RankMap<double>;
template class RankMap<float>;

} // namespace atrip
