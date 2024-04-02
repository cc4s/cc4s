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

// [[file:~/cuda/atrip/atrip.org::*Unions][Unions:1]]
#pragma once
#include <atrip/SliceUnion.hpp>

namespace atrip {

template <typename F = double>
static void slice_into_vector
#if defined(ATRIP_SOURCES_IN_GPU)
    (DataPtr<F> &source,
#else
    (std::vector<F> &source,
#endif
     size_t slice_size,
     CTF::Tensor<F> &to_slice,
     std::vector<int64_t> const low,
     std::vector<int64_t> const up,
     CTF::Tensor<F> const &origin,
     std::vector<int64_t> const origin_low,
     std::vector<int64_t> const origin_up) {
  // Thank you CTF for forcing me to do this
  struct {
    std::vector<int> up, low;
  } to_slice_ = {{up.begin(), up.end()}, {low.begin(), low.end()}},
    origin_ = {{origin_up.begin(), origin_up.end()},
               {origin_low.begin(), origin_low.end()}};

  WITH_OCD
  WITH_RANK << "slicing into " << pretty_print(to_slice_.up) << ","
            << pretty_print(to_slice_.low) << " from "
            << pretty_print(origin_.up) << "," << pretty_print(origin_.low)
            << "\n";

#if !defined(ATRIP_DONT_SLICE) && !defined(ATRIP_DRY)
  to_slice.slice(to_slice_.low.data(),
                 to_slice_.up.data(),
                 0.0,
                 origin,
                 origin_.low.data(),
                 origin_.up.data(),
                 1.0);

#  if defined(ATRIP_SOURCES_IN_GPU)
  WITH_CHRONO("cuda:sources",
              _CHECK_CUDA_SUCCESS(
                  "copying sources data to device",
                  cuMemcpyHtoD(source, to_slice.data, sizeof(F) * slice_size));)
#  else
  memcpy(source.data(), to_slice.data, sizeof(F) * slice_size);
#  endif /* defined(ATRIP_SOURCES_IN_GPU) */

#else
  IGNORABLE(source);
  IGNORABLE(slice_size);
  IGNORABLE(to_slice);
  IGNORABLE(origin);
#  pragma message("WARNING: COMPILING WITHOUT SLICING THE TENSORS")
#endif /* !defined(ATRIP_DONT_SLICE) && !defined(ATRIP_DRY) */
}

template <typename F = double>
struct TAPHH : public SliceUnion<F> {
  TAPHH(CTF::Tensor<F> const &source_tensor,
        size_t No,
        size_t Nv,
        size_t np,
        MPI_Comm child_world,
        MPI_Comm global_world)
      : SliceUnion<F>({Slice<F>::A, Slice<F>::B, Slice<F>::C},
                      {Nv, No, No} // size of the slices
                      ,
                      {Nv},
                      np,
                      child_world,
                      global_world,
                      Slice<F>::TA,
                      6) {
    this->init(source_tensor);
  }

  void slice_into_buffer(size_t it,
                         CTF::Tensor<F> &to,
                         CTF::Tensor<F> const &from) override {

    const int Nv = this->slice_length[0], No = this->slice_length[1],
              a = this->rank_map.find({static_cast<size_t>(Atrip::rank), it});

    slice_into_vector<F>(this->sources[it],
                         this->slice_size,
                         to,
                         {0, 0, 0},
                         {Nv, No, No},
                         from,
                         {a, 0, 0, 0},
                         {a + 1, Nv, No, No});
  }
};

template <typename F = double>
struct HHHA : public SliceUnion<F> {
  HHHA(CTF::Tensor<F> const &source_tensor,
       size_t No,
       size_t Nv,
       size_t np,
       MPI_Comm child_world,
       MPI_Comm global_world)
      : SliceUnion<F>({Slice<F>::A, Slice<F>::B, Slice<F>::C},
                      {No, No, No} // size of the slices
                      ,
                      {Nv} // size of the parametrization
                      ,
                      np,
                      child_world,
                      global_world,
                      Slice<F>::VIJKA,
                      6) {
    this->init(source_tensor);
  }

  void slice_into_buffer(size_t it,
                         CTF::Tensor<F> &to,
                         CTF::Tensor<F> const &from) override {

    const int No = this->slice_length[0],
              a = this->rank_map.find({static_cast<size_t>(Atrip::rank), it});

    slice_into_vector<F>(this->sources[it],
                         this->slice_size,
                         to,
                         {0, 0, 0},
                         {No, No, No},
                         from,
                         {0, 0, 0, a},
                         {No, No, No, a + 1});
  }
};

template <typename F = double>
struct ABPH : public SliceUnion<F> {
  ABPH(CTF::Tensor<F> const &source_tensor,
       size_t No,
       size_t Nv,
       size_t np,
       MPI_Comm child_world,
       MPI_Comm global_world)
      : SliceUnion<F>({Slice<F>::AB,
                       Slice<F>::BC,
                       Slice<F>::AC,
                       Slice<F>::BA,
                       Slice<F>::CB,
                       Slice<F>::CA},
                      {Nv, No} // size of the slices
                      ,
                      {Nv, Nv} // size of the parametrization
                      ,
                      np,
                      child_world,
                      global_world,
                      Slice<F>::VABCI,
                      2 * 6) {
    this->init(source_tensor);
  }

  void slice_into_buffer(size_t it,
                         CTF::Tensor<F> &to,
                         CTF::Tensor<F> const &from) override {

    const int Nv = this->slice_length[0], No = this->slice_length[1],
              el = this->rank_map.find({static_cast<size_t>(Atrip::rank), it}),
              a = el % Nv, b = el / Nv;

    slice_into_vector<F>(this->sources[it],
                         this->slice_size,
                         to,
                         {0, 0},
                         {Nv, No},
                         from,
                         {a, b, 0, 0},
                         {a + 1, b + 1, Nv, No});
  }
};

template <typename F = double>
struct ABHH : public SliceUnion<F> {
  ABHH(CTF::Tensor<F> const &source_tensor,
       size_t No,
       size_t Nv,
       size_t np,
       MPI_Comm child_world,
       MPI_Comm global_world)
      : SliceUnion<F>({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC},
                      {No, No} // size of the slices
                      ,
                      {Nv, Nv} // size of the parametrization
                      ,
                      np,
                      child_world,
                      global_world,
                      Slice<F>::VABIJ,
                      6) {
    this->init(source_tensor);
  }

  void slice_into_buffer(size_t it,
                         CTF::Tensor<F> &to,
                         CTF::Tensor<F> const &from) override {

    const int Nv = from.lens[0], No = this->slice_length[1],
              el = this->rank_map.find({static_cast<size_t>(Atrip::rank), it}),
              a = el % Nv, b = el / Nv;

    slice_into_vector<F>(this->sources[it],
                         this->slice_size,
                         to,
                         {0, 0},
                         {No, No},
                         from,
                         {a, b, 0, 0},
                         {a + 1, b + 1, No, No});
  }
};

template <typename F = double>
struct TABHH : public SliceUnion<F> {
  TABHH(CTF::Tensor<F> const &source_tensor,
        size_t No,
        size_t Nv,
        size_t np,
        MPI_Comm child_world,
        MPI_Comm global_world)
      : SliceUnion<F>({Slice<F>::AB, Slice<F>::BC, Slice<F>::AC},
                      {No, No} // size of the slices
                      ,
                      {Nv, Nv} // size of the parametrization
                      ,
                      np,
                      child_world,
                      global_world,
                      Slice<F>::TABIJ,
                      6) {
    this->init(source_tensor);
  }

  void slice_into_buffer(size_t it,
                         CTF::Tensor<F> &to,
                         CTF::Tensor<F> const &from) override {
    // TODO: maybe generalize this with ABHH

    const int Nv = from.lens[0], No = this->slice_length[1],
              el = this->rank_map.find({static_cast<size_t>(Atrip::rank), it}),
              a = el % Nv, b = el / Nv;

    slice_into_vector<F>(this->sources[it],
                         this->slice_size,
                         to,
                         {0, 0},
                         {No, No},
                         from,
                         {a, b, 0, 0},
                         {a + 1, b + 1, No, No});
  }
};

} // namespace atrip
// Unions:1 ends here
