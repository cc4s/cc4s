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
#include <set>
#include <unordered_set>

#include <atrip/Chrono.hpp>
#include <atrip/Types.hpp>
#include <atrip/Debug.hpp>
#include <atrip/Slice.hpp>
#include <atrip/RankMap.hpp>
#include <atrip/Utils.hpp>
#include <atrip/Malloc.hpp>

#if defined(ATRIP_SOURCES_IN_GPU)
#  define SOURCES_DATA(s) (s)
#else
#  define SOURCES_DATA(s) (s).data()
#endif

namespace atrip {
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:2]]
template <typename F = double>
class SliceUnion {
public:
  using Tensor = CTF::Tensor<F>;

  const RankMap<F> rank_map;
  const MPI_Comm world;
  const MPI_Comm universe;
  const std::vector<size_t> slice_length;
  const size_t slice_size;
#if defined(ATRIP_SOURCES_IN_GPU)
  std::vector<DataPtr<F>> sources;
#else
  std::vector<std::vector<F>> sources;
#endif
  std::vector<Slice<F>> slices;
  typename Slice<F>::Name name;
  const std::vector<typename Slice<F>::Type> slice_types;
  std::vector<DataPtr<F>> slice_buffers;
  std::set<DataPtr<F>> free_pointers;
#if defined(ATRIP_MPI_STAGING_BUFFERS)
  std::unordered_set<StagingBufferInfo, StagingBufferInfoHash>
      mpi_staging_buffers;
#endif

#if defined(ATRIP_MPI_STAGING_BUFFERS)
  struct StagingBufferInfo {
    DataPtr<F> data;
    size_t tag;
    MPI_Request *request;
    ABCTuple const abc;

    bool operator==(StagingBufferInfo const &o) const {
      // TODO: think about this more carefully,
      //     can two staging buffers have the same data, it should not be
      //     the case
      return o.data == data;
    }
  };

  class StagingBufferInfoHash {
  public:
    size_t operator()(StagingBufferInfo const &i) const {
      return (size_t)i.data;
    }
  };
#endif /* defined(ATRIP_MPI_STAGING_BUFFERS) */

  virtual void
  slice_into_buffer(size_t iteration, Tensor &to, Tensor const &from) = 0;

  /*
   * This function should enforce an important property of a SliceUnion.
   * Namely, there can be no two Slices of the same nature.
   *
   * This means that there can be at most one slice with a given Ty_x_Tu.
   */
  void check_for_duplicates() const;

  std::vector<typename Slice<F>::Ty_x_Tu>
  needed_slices_for_tuple(ABCTuple const &abc);

  /* build_local_database
   *
   * It should build a database of slices so that we know what is needed
   * to fetch in the next iteration represented by the tuple 'abc'.
   *
   * 1. The algorithm works as follows, we build a database of the all
   * the slice types that we need together with their tuple.
   *
   * 2. Look in the SliceUnion if we already have this tuple,
   * if we already have it mark it (TODO)
   *
   * 3. If we don't have the tuple, look for a (state=acceptor, type=blank)
   * slice and mark this slice as type=Fetch with the corresponding type
   * and tuple.
   *
   * NOTE: The algorithm should certify that we always have enough blank
   * slices.
   *
   */
  typename Slice<F>::LocalDatabase build_local_database(ABCTuple const &abc);

  /*
   * Garbage collect slices not needed for the next iteration.
   *
   * It will throw if it tries to gc a slice that has not been
   * previously unwrapped, as a safety mechanism.
   */
  void clear_unused_slices_for_next_tuple(ABCTuple const &abc);

  DataPtr<F> unwrap_slice(typename Slice<F>::Type type, ABCTuple const &abc);
  void unwrap_all(ABCTuple const &abc);

  ~SliceUnion();

  // static size_t get_size(const std::vector<size_t> &slice_length,
  //                        const std::vector<size_t> &param_length,
  //                        const size_t np);

  // CONSTRUCTOR
  SliceUnion(std::vector<typename Slice<F>::Type> slice_types_,
             std::vector<size_t> slice_length_,
             std::vector<size_t> param_length,
             size_t np,
             MPI_Comm child_world,
             MPI_Comm global_world,
             typename Slice<F>::Name name_,
             size_t n_slice_buffers = 4)
      : rank_map(param_length, np)
      , world(child_world)
      , universe(global_world)
      , slice_length(slice_length_)
      , slice_size(std::accumulate(slice_length.begin(),
                                   slice_length.end(),
                                   1UL,
                                   std::multiplies<size_t>()))

#if defined(ATRIP_SOURCES_IN_GPU)
      , sources(rank_map.n_sources())
#else
      , sources(rank_map.n_sources(),
#  if defined(ATRIP_DRY)
                std::vector<F>(1))
#  else
                std::vector<F>(slice_size))
#  endif /* defined(ATRIP_DRY) */
#endif   /* defined(ATRIP_SOURCES_IN_GPU) */
      , name(name_)
      , slice_types(slice_types_)
      , slice_buffers(n_slice_buffers) { // constructor begin

#if defined(ATRIP_SOURCES_IN_GPU) && defined(HAVE_CUDA)
    for (auto &ptr : sources) {
      MALLOC_DATA_POINTER("SOURCES", &ptr, sizeof(F) * slice_size);
    }
#endif

    for (auto &ptr : slice_buffers) {
      MALLOC_DATA_PTR("Slice Buffer", &ptr, sizeof(F) * slice_size);
    }

    slices = std::vector<Slice<F>>(2 * slice_types.size(), {slice_size});
    // TODO: think exactly    ^------------------- about this number

    // initialize the free_pointers with the pointers to the buffers
    std::transform(slice_buffers.begin(),
                   slice_buffers.end(),
                   std::inserter(free_pointers, free_pointers.begin()),
                   [](DataPtr<F> ptr) { return ptr; });

#if defined(HAVE_CUDA)
    LOG(1, "Atrip") << "warming communication up " << slices.size() << "\n";
    WITH_CHRONO(
        "cuda:warmup", int n_ranks = Atrip::np, request_count = 0;
        int n_sends = slice_buffers.size() * n_ranks;
        MPI_Request *requests =
            (MPI_Request *)malloc(n_sends * 2 * sizeof(MPI_Request));
        MPI_Status *statuses =
            (MPI_Status *)malloc(n_sends * 2 * sizeof(MPI_Status));
        for (int slice_id = 0; slice_id < slice_buffers.size(); slice_id++) {
          for (int rank_id = 0; rank_id < n_ranks; rank_id++) {
            MPI_Isend((void *)SOURCES_DATA(sources[0]),
                      slice_size,
                      traits::mpi::datatype_of<F>(),
                      rank_id,
                      100,
                      universe,
                      &requests[request_count++]);
            MPI_Irecv((void *)slice_buffers[slice_id],
                      slice_size,
                      traits::mpi::datatype_of<F>(),
                      rank_id,
                      100,
                      universe,
                      &requests[request_count++]);
          }
        } MPI_Waitall(n_sends * 2, requests, statuses);
        free(requests);
        free(statuses);)
#endif

    LOG(1, "Atrip") << "#slices " << slices.size() << "\n";
    WITH_RANK << "#slices[0] " << slices[0].size << "\n";
    LOG(1, "Atrip") << "#sources " << sources.size() << "\n";
    WITH_RANK << "#sources[0] " << slice_size << "\n";
    WITH_RANK << "#free_pointers " << free_pointers.size() << "\n";
    LOG(1, "Atrip") << "#slice_buffers " << slice_buffers.size() << "\n";
    LOG(1, "Atrip") << "GB*" << np << " "
                    << double(sources.size() + slice_buffers.size())
                           * slice_size * 8 * np / 1073741824.0
                    << "\n";
  } // constructor ends

  void init(Tensor const &source_tensor);

  /*
   */
  void allocate_free_buffer();

  DataPtr<F> pop_free_pointers();

  /**
   * \brief Send asynchronously only if the state is Fetch
   */
  void send(size_t other_rank,
            typename Slice<F>::LocalDatabaseElement const &el,
            size_t tag,
            ABCTuple abc);

  /**
   * \brief Receive asynchronously only if the state is Fetch
   */
  void receive(typename Slice<F>::Info const &info, size_t tag) noexcept;
}; // class SliceUnion

template <typename F>
SliceUnion<F> &union_by_name(std::vector<SliceUnion<F> *> const &unions,
                             typename Slice<F>::Name name);

} // namespace atrip
