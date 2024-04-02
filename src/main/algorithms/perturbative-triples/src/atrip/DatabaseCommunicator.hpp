#pragma once
#include <atrip/Utils.hpp>
#include <atrip/Equations.hpp>
#include <atrip/SliceUnion.hpp>
#include <atrip/Unions.hpp>

namespace atrip {

template <typename F>
using Unions = std::vector<SliceUnion<F> *>;

template <typename F>
typename Slice<F>::Database
naive_database(Unions<F> &unions, size_t nv, size_t np, size_t iteration);

} // namespace atrip
