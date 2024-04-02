#include <atrip/DatabaseCommunicator.hpp>
#include <atrip/Complex.hpp>

namespace atrip {

#if defined(ATRIP_NAIVE_SLOW)
/*
 * This function is really too slow, below are more performant
 * functions to get tuples.
 */
static ABCTuples get_nth_naive_tuples(size_t Nv, size_t np, int64_t i) {

  const size_t
      // total number of tuples for the problem
      n = Nv * (Nv + 1) * (Nv + 2) / 6 - Nv

      // all ranks should have the same number of tuples_per_rank
      ,
      tuples_per_rank = n / np + size_t(n % np != 0);

  ABCTuples result(np);
  if (i < 0) return result;
  std::vector<size_t> rank_indices(np, 0);

  for (size_t a(0), g(0); a < Nv; a++)
    for (size_t b(a); b < Nv; b++)
      for (size_t c(b); c < Nv; c++) {
        if (a == b && b == c) continue;
        for (size_t rank = 0; rank < np; rank++) {

          const size_t
              // start index for the global tuples list
              start = tuples_per_rank * rank

              // end index for the global tuples list
              ,
              end = tuples_per_rank * (rank + 1);

          if (start <= g && g < end) {
            if (rank_indices[rank] == i) { result[rank] = {a, b, c}; }
            rank_indices[rank] += 1;
          }
        }
        g++;
      }

  return result;
}
#endif

static inline size_t a_block_sum_atrip(int64_t T, int64_t nv) {
  const int64_t nv_min_1 = nv - 1, t_plus_1 = T + 1;
  return t_plus_1 * nv_min_1 * nv + nv_min_1 * t_plus_1
       - (nv_min_1 * (T * t_plus_1) / 2)
       - (t_plus_1 * (nv_min_1 * nv) / 2)
       // do not simplify this expression, only the addition of both parts
       // is a pair integer, prepare to endure the consequences of
       // simplifying otherwise
       + (((T * t_plus_1 * (1 + 2 * T)) / 6) - 3 * ((T * t_plus_1) / 2)) / 2;
}

static inline int64_t b_block_sum_atrip(int64_t a, int64_t T, int64_t nv) {
  return nv * ((T - a) + 1) - (T * (T + 1) - a * (a - 1)) / 2 - 1;
}

static std::vector<size_t> a_sums;
static inline ABCTuple nth_atrip(size_t it, size_t nv) {

  // build the sums if necessary
  if (!a_sums.size()) {
    a_sums.resize(nv);
    for (size_t _i = 0; _i < nv; _i++) {
      a_sums[_i] = a_block_sum_atrip(_i, nv);
    }
  }

  int64_t a = -1;
  size_t block_a = 0;
  for (const auto &sum : a_sums) {
    ++a;
    if (sum > it) {
      break;
    } else {
      block_a = sum;
    }
  }

  // build the b_sums
  std::vector<int64_t> b_sums(nv - a);
  for (size_t t = a, i = 0; t < nv; t++) {
    b_sums[i++] = b_block_sum_atrip(a, t, nv);
  }
  int64_t b = a - 1, block_b = block_a;
  for (const auto &sum : b_sums) {
    ++b;
    if (sum + block_a > it) {
      break;
    } else {
      block_b = sum + block_a;
    }
  }

  const int64_t c = b + it - block_b + (a == b);

  return {(size_t)a, (size_t)b, (size_t)c};
}

static inline ABCTuples
nth_atrip_distributed(int64_t it, size_t nv, size_t np) {

  // If we are getting the previous tuples in the first iteration,
  // then just return an impossible tuple, different from the FAKE_TUPLE,
  // because if FAKE_TUPLE is defined as {0,0,0} slices thereof
  // are actually attainable.
  //
  if (it < 0) {
    ABCTuples result(np, {nv, nv, nv});
    return result;
  }

  ABCTuples result(np);

  const size_t
      // total number of tuples for the problem
      n = nv * (nv + 1) * (nv + 2) / 6 - nv

      // all ranks should have the same number of tuples_per_rank
      ,
      tuples_per_rank = n / np + size_t(n % np != 0);

  for (size_t rank = 0; rank < np; rank++) {
    const size_t global_iteration = tuples_per_rank * rank + it;
    result[rank] = nth_atrip(global_iteration, nv);
  }

  return result;
}

template <typename F>
static typename Slice<F>::LocalDatabase
build_local_database_fake(ABCTuple const &abc_prev,
                          ABCTuple const &abc,
                          size_t rank,
                          SliceUnion<F> *u) {

  typename Slice<F>::LocalDatabase result;

  // vector of type x tuple
  auto const needed = u->needed_slices_for_tuple(abc);
  auto const needed_prev = u->needed_slices_for_tuple(abc_prev);

  for (auto const &pair : needed) {
    auto const type = pair.first;
    auto const tuple = pair.second;
    auto const from = u->rank_map.find(abc, type);

    // Try to find in the previously needed slices
    // one that exactly matches the tuple.
    // Not necessarily has to match the type.
    //
    // If we find it, then it means that the fake rank
    // will mark it as recycled. This covers
    // the finding of Ready slices and Recycled slices.
    {
      auto const &it =
          std::find_if(needed_prev.begin(),
                       needed_prev.end(),
                       [&tuple, &type](typename Slice<F>::Ty_x_Tu const &o) {
                         return o.second == tuple;
                       });

      if (it != needed_prev.end()) {
        typename Slice<F>::Info info;
        info.tuple = tuple;
        info.type = type;
        info.from = from;
        info.state = Slice<F>::Recycled;
        result.push_back({u->name, info});
        continue;
      }
    }

    {
      typename Slice<F>::Info info;
      info.type = type;
      info.tuple = tuple;
      info.from = from;

      // Handle self sufficiency
      info.state =
          rank == from.rank ? Slice<F>::SelfSufficient : Slice<F>::Fetch;
      result.push_back({u->name, info});
      continue;
    }
  }

  return result;
}

template <typename F>
typename Slice<F>::Database
naive_database(Unions<F> &unions, size_t nv, size_t np, size_t iteration) {

  using Database = typename Slice<F>::Database;
  Database db;

#ifdef ATRIP_NAIVE_SLOW
  WITH_CHRONO("db:comm:naive:tuples",
              const auto tuples = get_nth_naive_tuples(nv, np, iteration);
              const auto prev_tuples =
                  get_nth_naive_tuples(nv, np, iteration - 1);)
#else
  WITH_CHRONO("db:comm:naive:tuples",
              const auto tuples = nth_atrip_distributed(iteration, nv, np);
              const auto prev_tuples =
                  nth_atrip_distributed(iteration - 1, nv, np);)

#endif

  for (size_t rank = 0; rank < np; rank++) {
    auto abc = tuples[rank];
    typename Slice<F>::LocalDatabase ldb;

    for (auto const &tensor : unions) {
      if (rank == Atrip::rank) {
        auto const &tensor_db = tensor->build_local_database(abc);
        ldb.insert(ldb.end(), tensor_db.begin(), tensor_db.end());
      } else {
        auto const &tensor_db =
            build_local_database_fake(prev_tuples[rank], abc, rank, tensor);
        ldb.insert(ldb.end(), tensor_db.begin(), tensor_db.end());
      }
    }

    db.insert(db.end(), ldb.begin(), ldb.end());
  }

  return db;
}

template typename Slice<double>::Database
naive_database<double>(Unions<double> &unions,
                       size_t nv,
                       size_t np,
                       size_t iteration);

template typename Slice<Complex>::Database
naive_database<Complex>(Unions<Complex> &unions,
                        size_t nv,
                        size_t np,
                        size_t iteration);

} // namespace atrip
