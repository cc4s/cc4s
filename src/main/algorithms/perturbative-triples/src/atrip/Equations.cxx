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

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:2]]
#include <cstring>

#include <atrip/Equations.hpp>

#include <atrip/CUDA.hpp>
#include <atrip/Operations.hpp>

namespace atrip {
// Prolog:2 ends here

#if defined(HAVE_CUDA)
#  define FOR_K()                                                              \
    const size_t k = blockIdx.x * blockDim.x + threadIdx.x;                    \
    size_t idx = k * size * size;
#else
#  define FOR_K() for (size_t k = 0, idx = 0; k < size; k++)
#endif
#define _IJK_(i, j, k) i + j *size + k *size *size
#define _REORDER_BODY_(...)                                                    \
  FOR_K()                                                                      \
  for (size_t j = 0; j < size; j++)                                            \
    for (size_t i = 0; i < size; i++, idx++) { __VA_ARGS__ }
#define _MAKE_REORDER_(_enum, ...)                                             \
  template <typename F>                                                        \
  __MAYBE_GLOBAL__ void reorder(reorder_proxy<F, _enum> p,                     \
                                size_t size,                                   \
                                F *to,                                         \
                                F *from) {                                     \
    IGNORABLE(p);                                                              \
    _REORDER_BODY_(__VA_ARGS__)                                                \
  }
#if defined(HAVE_CUDA)
#  define GO(__TO, __FROM) acc::sum_in_place<F>(&__TO, &__FROM);
#else
#  define GO(__TO, __FROM) __TO += __FROM;
#endif

// These are just help structures
// to help with the templating of reorder
// function
enum reordering_t {
  IJK,
  IKJ,
  JIK,
  JKI,
  KIJ,
  KJI
};

/*
 * Please the c++ type checker and template creator
 * in order to have an argument in the signature of
 * the function  that helps the compiler know which
 * instantiation it should take.
 */
template <typename F, reordering_t R>
struct reorder_proxy {};

template <typename F, reordering_t R>
__MAYBE_GLOBAL__ void
reorder(reorder_proxy<F, R> proxy, size_t size, F *to, F *from);

_MAKE_REORDER_(IJK, GO(to[idx], from[_IJK_(i, j, k)]))
_MAKE_REORDER_(IKJ, GO(to[idx], from[_IJK_(i, k, j)]))
_MAKE_REORDER_(JIK, GO(to[idx], from[_IJK_(j, i, k)]))
_MAKE_REORDER_(JKI, GO(to[idx], from[_IJK_(j, k, i)]))
_MAKE_REORDER_(KIJ, GO(to[idx], from[_IJK_(k, i, j)]))
_MAKE_REORDER_(KJI, GO(to[idx], from[_IJK_(k, j, i)]))

#undef LIMS_KS
#undef _MAKE_REORDER
#undef _REORDER_BODY_
#undef _IJK_
#undef GO

#if defined(HAVE_CUDA)
#  define MIN(a, b) min((a), (b))
#else
#  define MIN(a, b) std::min((a), (b))
#endif

#if defined(ATRIP_NEW_ENERGY)

// [[file:~/cuda/atrip/atrip.org::*Energy][Energy:2]]
template <typename F>
__MAYBE_GLOBAL__ void get_energy_distinct(F const epsabc,
                                          size_t const No,
                                          F *const epsi,
                                          F *const Tijk,
                                          F *const Zijk,
                                          double *energy) {
  constexpr size_t block_size = 16;
  F _energy = {0.};
  for (size_t kk = 0; kk < No; kk += block_size) {
    const size_t kend(MIN(No, kk + block_size));
    for (size_t jj(kk); jj < No; jj += block_size) {
      const size_t jend(MIN(No, jj + block_size));
      for (size_t ii(jj); ii < No; ii += block_size) {
        const size_t iend(MIN(No, ii + block_size));
        for (size_t k(kk); k < kend; k++) {
          const F ek(epsi[k]);
          const size_t jstart = jj > k ? jj : k;
          for (size_t j(jstart); j < jend; j++) {
            F const ej(epsi[j]);
            F const facjk = j == k ? F{0.5} : F{1.0};
            size_t istart = ii > j ? ii : j;
            for (size_t i(istart); i < iend; i++) {
              const F ei(epsi[i]),
                  facij = i == j ? F{0.5} : F{1.0},
                  eijk(acc::add(acc::add(ei, ej), ek)),
                  denominator(acc::sub(epsabc, eijk)),
                  U(Zijk[i + No * j + No * No * k]),
                  V(Zijk[i + No * k + No * No * j]),
                  W(Zijk[j + No * i + No * No * k]),
                  X(Zijk[j + No * k + No * No * i]),
                  Y(Zijk[k + No * i + No * No * j]),
                  Z(Zijk[k + No * j + No * No * i]),
                  A(acc::maybe_conjugate_scalar(
                      Tijk[i + No * j + No * No * k])),
                  B(acc::maybe_conjugate_scalar(
                      Tijk[i + No * k + No * No * j])),
                  C(acc::maybe_conjugate_scalar(
                      Tijk[j + No * i + No * No * k])),
                  D(acc::maybe_conjugate_scalar(
                      Tijk[j + No * k + No * No * i])),
                  E(acc::maybe_conjugate_scalar(
                      Tijk[k + No * i + No * No * j])),
                  _F(acc::maybe_conjugate_scalar(
                      Tijk[k + No * j + No * No * i])),
                  AU = acc::prod(A, U), BV = acc::prod(B, V),
                  CW = acc::prod(C, W), DX = acc::prod(D, X),
                  EY = acc::prod(E, Y), FZ = acc::prod(_F, Z),
                  UXY = acc::add(U, acc::add(X, Y)),
                  VWZ = acc::add(V, acc::add(W, Z)),
                  ADE = acc::add(A, acc::add(D, E)),
                  BCF = acc::add(B, acc::add(C, _F))
                  // I just might as well write this in CL
                  ,
                  _first = acc::add(
                      AU,
                      acc::add(BV,
                               acc::add(CW, acc::add(DX, acc::add(EY, FZ))))),
                  _second =
                      acc::prod(acc::sub(UXY, acc::prod(F{-2.0}, VWZ)), ADE),
                  _third =
                      acc::prod(acc::sub(VWZ, acc::prod(F{-2.0}, UXY)), BCF),
                  value = acc::add(acc::prod(F{3.0}, _first),
                                   acc::add(_second, _third)),
                  _loop_energy =
                      acc::prod(acc::prod(F{2.0}, value),
                                acc::div(acc::prod(facjk, facij), denominator));
              acc::sum_in_place(&_energy, &_loop_energy);
            } // i
          }   // j
        }     // k
      }       // ii
    }         // jj
  }           // kk
  const double real_part = acc::real(_energy);
  acc::sum_in_place(energy, &real_part);
}

template <typename F>
__MAYBE_GLOBAL__ void get_energy_same(F const epsabc,
                                      size_t const No,
                                      F *const epsi,
                                      F *const Tijk,
                                      F *const Zijk,
                                      double *energy) {
  constexpr size_t block_size = 16;
  F _energy = F{0.};
  for (size_t kk = 0; kk < No; kk += block_size) {
    const size_t kend(MIN(kk + block_size, No));
    for (size_t jj(kk); jj < No; jj += block_size) {
      const size_t jend(MIN(jj + block_size, No));
      for (size_t ii(jj); ii < No; ii += block_size) {
        const size_t iend(MIN(ii + block_size, No));
        for (size_t k(kk); k < kend; k++) {
          const F ek(epsi[k]);
          const size_t jstart = jj > k ? jj : k;
          for (size_t j(jstart); j < jend; j++) {
            const F facjk(j == k ? F{0.5} : F{1.0});
            const F ej(epsi[j]);
            const size_t istart = ii > j ? ii : j;
            for (size_t i(istart); i < iend; i++) {
              const F ei(epsi[i]), facij(i == j ? F{0.5} : F{1.0}),
                  eijk(acc::add(acc::add(ei, ej), ek)),
                  denominator(acc::sub(epsabc, eijk)),
                  U(Zijk[i + No * j + No * No * k]),
                  V(Zijk[j + No * k + No * No * i]),
                  W(Zijk[k + No * i + No * No * j]),
                  A(acc::maybe_conjugate_scalar(
                      Tijk[i + No * j + No * No * k])),
                  B(acc::maybe_conjugate_scalar(
                      Tijk[j + No * k + No * No * i])),
                  C(acc::maybe_conjugate_scalar(
                      Tijk[k + No * i + No * No * j])),
                  ABC = acc::add(A, acc::add(B, C)),
                  UVW = acc::add(U, acc::add(V, W)), AU = acc::prod(A, U),
                  BV = acc::prod(B, V), CW = acc::prod(C, W),
                  AU_and_BV_and_CW = acc::add(acc::add(AU, BV), CW),
                  value = acc::sub(acc::prod(F{3.0}, AU_and_BV_and_CW),
                                   acc::prod(ABC, UVW)),
                  _loop_energy =
                      acc::prod(acc::prod(F{2.0}, value),
                                acc::div(acc::prod(facjk, facij), denominator));

              acc::sum_in_place(&_energy, &_loop_energy);
            } // i
          }   // j
        }     // k
      }       // ii
    }         // jj
  }           // kk
  const double real_part = acc::real(_energy);
  acc::sum_in_place(energy, &real_part);
}
// Energy:2 ends here

#else

// [[file:~/cuda/atrip/atrip.org::*Energy][Energy:2]]
template <typename F>
__MAYBE_GLOBAL__ void get_energy_distinct(F const epsabc,
                                          size_t const No,
                                          F *const epsi,
                                          F *const Tijk,
                                          F *const Zijk,
                                          double *_energy) {
  constexpr size_t block_size = 16;
  F energy(0.);
  for (size_t kk = 0; kk < No; kk += block_size) {
    const size_t kend(MIN(No, kk + block_size));
    for (size_t jj(kk); jj < No; jj += block_size) {
      const size_t jend(MIN(No, jj + block_size));
      for (size_t ii(jj); ii < No; ii += block_size) {
        const size_t iend(MIN(No, ii + block_size));
        for (size_t k(kk); k < kend; k++) {
          const F ek(epsi[k]);
          const size_t jstart = jj > k ? jj : k;
          for (size_t j(jstart); j < jend; j++) {
            F const ej(epsi[j]);
            F const facjk = j == k ? F(0.5) : F(1.0);
            size_t istart = ii > j ? ii : j;
            for (size_t i(istart); i < iend; i++) {
              const F ei(epsi[i]),
                  facij = i == j ? F(0.5) : F(1.0),
                  denominator(epsabc - ei - ej - ek),
                  U(Zijk[i + No * j + No * No * k]),
                  V(Zijk[i + No * k + No * No * j]),
                  W(Zijk[j + No * i + No * No * k]),
                  X(Zijk[j + No * k + No * No * i]),
                  Y(Zijk[k + No * i + No * No * j]),
                  Z(Zijk[k + No * j + No * No * i]),
                  A(acc::maybe_conjugate_scalar<F>(
                      Tijk[i + No * j + No * No * k])),
                  B(acc::maybe_conjugate_scalar<F>(
                      Tijk[i + No * k + No * No * j])),
                  C(acc::maybe_conjugate_scalar<F>(
                      Tijk[j + No * i + No * No * k])),
                  D(acc::maybe_conjugate_scalar<F>(
                      Tijk[j + No * k + No * No * i])),
                  E(acc::maybe_conjugate_scalar<F>(
                      Tijk[k + No * i + No * No * j])),
                  _F(acc::maybe_conjugate_scalar<F>(
                      Tijk[k + No * j + No * No * i])),
                  value = 3.0 * (A * U + B * V + C * W + D * X + E * Y + _F * Z)
                        + ((U + X + Y) - 2.0 * (V + W + Z)) * (A + D + E)
                        + ((V + W + Z) - 2.0 * (U + X + Y)) * (B + C + _F);
              energy += 2.0 * value / denominator * facjk * facij;
            } // i
          }   // j
        }     // k
      }       // ii
    }         // jj
  }           // kk
  *_energy = acc::real(energy);
}

template <typename F>
__MAYBE_GLOBAL__ void get_energy_same(F const epsabc,
                                      size_t const No,
                                      F *const epsi,
                                      F *const Tijk,
                                      F *const Zijk,
                                      double *_energy) {
  constexpr size_t block_size = 16;
  F energy = F(0.);
  for (size_t kk = 0; kk < No; kk += block_size) {
    const size_t kend(MIN(kk + block_size, No));
    for (size_t jj(kk); jj < No; jj += block_size) {
      const size_t jend(MIN(jj + block_size, No));
      for (size_t ii(jj); ii < No; ii += block_size) {
        const size_t iend(MIN(ii + block_size, No));
        for (size_t k(kk); k < kend; k++) {
          const F ek(epsi[k]);
          const size_t jstart = jj > k ? jj : k;
          for (size_t j(jstart); j < jend; j++) {
            const F facjk(j == k ? F(0.5) : F(1.0));
            const F ej(epsi[j]);
            const size_t istart = ii > j ? ii : j;
            for (size_t i(istart); i < iend; i++) {
              const F ei(epsi[i]), facij(i == j ? F(0.5) : F(1.0)),
                  denominator(epsabc - ei - ej - ek),
                  U(Zijk[i + No * j + No * No * k]),
                  V(Zijk[j + No * k + No * No * i]),
                  W(Zijk[k + No * i + No * No * j]),
                  A(acc::maybe_conjugate_scalar<F>(
                      Tijk[i + No * j + No * No * k])),
                  B(acc::maybe_conjugate_scalar<F>(
                      Tijk[j + No * k + No * No * i])),
                  C(acc::maybe_conjugate_scalar<F>(
                      Tijk[k + No * i + No * No * j])),
                  value = F(3.0) * (A * U + B * V + C * W)
                        - (A + B + C) * (U + V + W);
              energy += F(2.0) * value / denominator * facjk * facij;
            } // i
          }   // j
        }     // k
      }       // ii
    }         // jj
  }           // kk
  *_energy = acc::real(energy);
}
// Energy:2 ends here
#endif /* defined(ATRIP_NEW_ENERGY) */

// [[file:~/cuda/atrip/atrip.org::*Energy][Energy:3]]
// instantiate double
template __MAYBE_GLOBAL__ void
get_energy_distinct(DataFieldType<double> const epsabc,
                    size_t const No,
                    DataFieldType<double> *const epsi,
                    DataFieldType<double> *const Tijk,
                    DataFieldType<double> *const Zijk,
                    DataFieldType<double> *energy);

template __MAYBE_GLOBAL__ void
get_energy_same(DataFieldType<double> const epsabc,
                size_t const No,
                DataFieldType<double> *const epsi,
                DataFieldType<double> *const Tijk,
                DataFieldType<double> *const Zijk,
                DataFieldType<double> *energy);

// instantiate Complex
template __MAYBE_GLOBAL__ void
get_energy_distinct(DataFieldType<Complex> const epsabc,
                    size_t const No,
                    DataFieldType<Complex> *const epsi,
                    DataFieldType<Complex> *const Tijk,
                    DataFieldType<Complex> *const Zijk,
                    DataFieldType<double> *energy);

template __MAYBE_GLOBAL__ void
get_energy_same(DataFieldType<Complex> const epsabc,
                size_t const No,
                DataFieldType<Complex> *const epsi,
                DataFieldType<Complex> *const Tijk,
                DataFieldType<Complex> *const Zijk,
                DataFieldType<double> *energy);
// Energy:3 ends here

// [[file:~/cuda/atrip/atrip.org::*Singles%20contribution][Singles
// contribution:2]]
template <typename F>
__MAYBE_GLOBAL__ void singles_contribution(size_t No,
                                           size_t Nv,
                                           size_t a,
                                           size_t b,
                                           size_t c,
                                           DataFieldType<F> *const Tph,
                                           DataFieldType<F> *const VABij,
                                           DataFieldType<F> *const VACij,
                                           DataFieldType<F> *const VBCij,
                                           DataFieldType<F> *Zijk) {
  const size_t NoNo = No * No;
  // TODO: change order of for loops
  for (size_t k = 0; k < No; k++)
    for (size_t i = 0; i < No; i++)
      for (size_t j = 0; j < No; j++) {
        const size_t ijk = i + j * No + k * NoNo;

#if defined(HAVE_CUDA)

#  define GO(__TPH, __VABIJ)                                                   \
    do {                                                                       \
      const DataFieldType<F> product =                                         \
          acc::prod<DataFieldType<F>>((__TPH), (__VABIJ));                     \
      acc::sum_in_place<DataFieldType<F>>(&Zijk[ijk], &product);               \
    } while (0)

#else

#  define GO(__TPH, __VABIJ) Zijk[ijk] += (__TPH) * (__VABIJ)

#endif // HAVE_CUDA

        GO(Tph[a + i * Nv], VBCij[j + k * No]);
        GO(Tph[b + j * Nv], VACij[i + k * No]);
        GO(Tph[c + k * Nv], VABij[i + j * No]);

#undef GO
      } // for loop j
}

// instantiate
template __MAYBE_GLOBAL__ void singles_contribution<double>(size_t No,
                                                            size_t Nv,
                                                            size_t a,
                                                            size_t b,
                                                            size_t c,
                                                            double *const Tph,
                                                            double *const VABij,
                                                            double *const VACij,
                                                            double *const VBCij,
                                                            double *Zijk);

template __MAYBE_GLOBAL__ void
singles_contribution<Complex>(size_t No,
                              size_t Nv,
                              size_t a,
                              size_t b,
                              size_t c,
                              DataFieldType<Complex> *const Tph,
                              DataFieldType<Complex> *const VABij,
                              DataFieldType<Complex> *const VACij,
                              DataFieldType<Complex> *const VBCij,
                              DataFieldType<Complex> *Zijk);
// Singles contribution:2 ends here

// [[file:~/cuda/atrip/atrip.org::*Doubles%20contribution][Doubles
// contribution:2]]
template <typename F>
void doubles_contribution(size_t const No,
                          size_t const Nv
                          // -- VABCI
                          ,
                          DataPtr<F> const VABph,
                          DataPtr<F> const VACph,
                          DataPtr<F> const VBCph,
                          DataPtr<F> const VBAph,
                          DataPtr<F> const VCAph,
                          DataPtr<F> const VCBph
                          // -- VHHHA
                          ,
                          DataPtr<F> const VhhhA,
                          DataPtr<F> const VhhhB,
                          DataPtr<F> const VhhhC
                          // -- TA
                          ,
                          DataPtr<F> const TAphh,
                          DataPtr<F> const TBphh,
                          DataPtr<F> const TCphh
                          // -- TABIJ
                          ,
                          DataPtr<F> const TABhh,
                          DataPtr<F> const TAChh,
                          DataPtr<F> const TBChh
                          // -- TIJK
                          // , DataPtr<F> Tijk_
                          ,
                          DataFieldType<F> *Tijk_,
                          // -- tmp buffers
                          DataFieldType<F> *_t_buffer,
                          DataFieldType<F> *_vhhh) {
  const size_t NoNo = No * No;

  DataFieldType<F> *Tijk = (DataFieldType<F> *)Tijk_;

#if defined(ATRIP_USE_DGEMM)
#  if defined(HAVE_CUDA)
#    define REORDER(__II, __JJ, __KK)                                          \
      reorder<<<1, No>>>(reorder_proxy<DataFieldType<F>, __II##__JJ##__KK>{},  \
                         No,                                                   \
                         Tijk,                                                 \
                         _t_buffer)
#    define DGEMM_PARTICLES(__A, __B)                                          \
      atrip::xgemm<F>("T",                                                     \
                      "N",                                                     \
                      (int const *)&NoNo,                                      \
                      (int const *)&No,                                        \
                      (int const *)&Nv,                                        \
                      &one,                                                    \
                      (DataFieldType<F> *)__A,                                 \
                      (int const *)&Nv,                                        \
                      (DataFieldType<F> *)__B,                                 \
                      (int const *)&Nv,                                        \
                      &zero,                                                   \
                      _t_buffer,                                               \
                      (int const *)&NoNo)
#    define DGEMM_HOLES(__A, __B, __TRANSB)                                    \
      atrip::xgemm<F>("N",                                                     \
                      __TRANSB,                                                \
                      (int const *)&NoNo,                                      \
                      (int const *)&No,                                        \
                      (int const *)&No,                                        \
                      &m_one,                                                  \
                      __A,                                                     \
                      (int const *)&NoNo,                                      \
                      (DataFieldType<F> *)__B,                                 \
                      (int const *)&No,                                        \
                      &zero,                                                   \
                      _t_buffer,                                               \
                      (int const *)&NoNo)
#    define MAYBE_CONJ(_conj, _buffer)                                         \
      do {                                                                     \
        acc::maybe_conjugate<<<1, 1>>>((DataFieldType<F> *)_conj,              \
                                       (DataFieldType<F> *)_buffer,            \
                                       NoNoNo);                                \
      } while (0)

  // END CUDA
  // ////////////////////////////////////////////////////////////////////

#  else

  // NONCUDA
  // /////////////////////////////////////////////////////////////////////

#    define REORDER(__II, __JJ, __KK)                                          \
      reorder(reorder_proxy<DataFieldType<F>, __II##__JJ##__KK>{},             \
              No,                                                              \
              Tijk,                                                            \
              _t_buffer)
#    define DGEMM_PARTICLES(__A, __B)                                          \
      atrip::xgemm<F>("T",                                                     \
                      "N",                                                     \
                      (int const *)&NoNo,                                      \
                      (int const *)&No,                                        \
                      (int const *)&Nv,                                        \
                      &one,                                                    \
                      __A,                                                     \
                      (int const *)&Nv,                                        \
                      __B,                                                     \
                      (int const *)&Nv,                                        \
                      &zero,                                                   \
                      _t_buffer,                                               \
                      (int const *)&NoNo)
#    define DGEMM_HOLES(__A, __B, __TRANSB)                                    \
      atrip::xgemm<F>("N",                                                     \
                      __TRANSB,                                                \
                      (int const *)&NoNo,                                      \
                      (int const *)&No,                                        \
                      (int const *)&No,                                        \
                      &m_one,                                                  \
                      __A,                                                     \
                      (int const *)&NoNo,                                      \
                      __B,                                                     \
                      (int const *)&No,                                        \
                      &zero,                                                   \
                      _t_buffer,                                               \
                      (int const *)&NoNo)
#    define MAYBE_CONJ(_conj, _buffer)                                         \
      acc::maybe_conjugate((DataFieldType<F> *)_conj,                          \
                           (DataFieldType<F> *)_buffer,                        \
                           NoNoNo);
#  endif

  F one{1.0}, m_one{-1.0}, zero{0.0};
  const size_t NoNoNo = No * NoNo;

// !!!! Zeroing vectors Tijk, _t_buffer and _vhhh
#  if defined(HAVE_CUDA)

#    if !defined(ATRIP_ONLY_DGEMM)
  {
    const size_t elements = NoNoNo * sizeof(DataFieldType<F>) / 4;
    WITH_CHRONO(
        "double:zeroing",
        _CHECK_CUDA_SUCCESS("Zeroing Tijk",
                            cuMemsetD32_v2((CUdeviceptr)Tijk, 0x00, elements));
        _CHECK_CUDA_SUCCESS(
            "Zeroing t buffer",
            cuMemsetD32_v2((CUdeviceptr)_t_buffer, 0x00, elements));
        _CHECK_CUDA_SUCCESS(
            "Zeroing vhhh buffer",
            cuMemsetD32_v2((CUdeviceptr)_vhhh, 0x00, elements));)
  }
#    endif

#  else
  std::memset((void *)_t_buffer, 0x00, NoNoNo * sizeof(DataFieldType<F>));
  std::memset((void *)_vhhh, 0x00, NoNoNo * sizeof(DataFieldType<F>));
  std::memset((void *)Tijk, 0x00, NoNoNo * sizeof(DataFieldType<F>));
#  endif /* HAVE_CUDA */

#  if defined(ATRIP_ONLY_DGEMM)
#    undef MAYBE_CONJ
#    undef REORDER
#    define MAYBE_CONJ(a, b)                                                   \
      do {                                                                     \
      } while (0)
#    define REORDER(i, j, k)                                                   \
      do {                                                                     \
      } while (0)
#  endif /* defined(ATRIP_ONLY_DGEMM) */

  // HOLES
  WITH_CHRONO("doubles:holes", {
    // VhhhC[i + k*No + L*NoNo] * TABhh[L + j*No]; H1
    MAYBE_CONJ(_vhhh, VhhhC);
    WITH_CHRONO("doubles:holes:1", /**/
                DGEMM_HOLES(_vhhh, TABhh, "N");
                REORDER(I, K, J);)
    // VhhhC[j + k*No + L*NoNo] * TABhh[i + L*No]; H0
    WITH_CHRONO("doubles:holes:2", /**/
                DGEMM_HOLES(_vhhh, TABhh, "T");
                REORDER(J, K, I);)

    // VhhhB[i + j*No + L*NoNo] * TAChh[L + k*No]; H5
    MAYBE_CONJ(_vhhh, VhhhB);
    WITH_CHRONO("doubles:holes:3", /**/
                DGEMM_HOLES(_vhhh, TAChh, "N");
                REORDER(I, J, K);)
    // VhhhB[k + j*No + L*NoNo] * TAChh[i + L*No]; H3
    WITH_CHRONO("doubles:holes:4", /**/
                DGEMM_HOLES(_vhhh, TAChh, "T");
                REORDER(K, J, I);)

    // VhhhA[j + i*No + L*NoNo] * TBChh[L + k*No]; H1
    MAYBE_CONJ(_vhhh, VhhhA);
    WITH_CHRONO("doubles:holes:5", /**/
                DGEMM_HOLES(_vhhh, TBChh, "N");
                REORDER(J, I, K);)
    // VhhhA[k + i*No + L*NoNo] * TBChh[j + L*No]; H4
    WITH_CHRONO("doubles:holes:6", /**/
                DGEMM_HOLES(_vhhh, TBChh, "T");
                REORDER(K, I, J);)
  })
#  undef MAYBE_CONJ

  // PARTICLES
  WITH_CHRONO("doubles:particles", {
    // TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv]; P0
    WITH_CHRONO("doubles:particles:1", /**/
                DGEMM_PARTICLES(TAphh, VBCph);
                REORDER(I, J, K);)
    // TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv]; P3
    WITH_CHRONO("doubles:particles:2", /**/
                DGEMM_PARTICLES(TAphh, VCBph);
                REORDER(I, K, J);)
    // TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv]; P5
    WITH_CHRONO("doubles:particles:3", /**/
                DGEMM_PARTICLES(TCphh, VABph);
                REORDER(K, I, J);)
    // TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv]; P2
    WITH_CHRONO("doubles:particles:4", /**/
                DGEMM_PARTICLES(TCphh, VBAph);
                REORDER(K, J, I);)
    // TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv]; P1
    WITH_CHRONO("doubles:particles:5", /**/
                DGEMM_PARTICLES(TBphh, VACph);
                REORDER(J, I, K);)
    // TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv]; P4
    WITH_CHRONO("doubles:particles:6", /**/
                DGEMM_PARTICLES(TBphh, VCAph);
                REORDER(J, K, I);)
  })

#  undef REORDER
#  undef DGEMM_HOLES
#  undef DGEMM_PARTICLES
#else
  const size_t NoNv = No * Nv;
  for (size_t k = 0; k < No; k++)
    for (size_t j = 0; j < No; j++)
      for (size_t i = 0; i < No; i++) {
        const size_t ijk = i + j * No + k * NoNo, jk = j + k * No;
        Tijk[ijk] = 0.0; // :important
        // HOLE DIAGRAMS: TABHH and VHHHA
        for (size_t L = 0; L < No; L++) {
          // t[abLj] * V[Lcik]        H1
          // t[baLi] * V[Lcjk]        H0      TODO: conjugate T for complex
          Tijk[ijk] -= TABhh[L + j * No] * VhhhC[i + k * No + L * NoNo];
          Tijk[ijk] -= TABhh[i + L * No] * VhhhC[j + k * No + L * NoNo];

          // t[acLk] * V[Lbij]        H5
          // t[caLi] * V[Lbkj]        H3
          Tijk[ijk] -= TAChh[L + k * No] * VhhhB[i + j * No + L * NoNo];
          Tijk[ijk] -= TAChh[i + L * No] * VhhhB[k + j * No + L * NoNo];

          // t[bcLk] * V[Laji]        H2
          // t[cbLj] * V[Laki]        H4
          Tijk[ijk] -= TBChh[L + k * No] * VhhhA[j + i * No + L * NoNo];
          Tijk[ijk] -= TBChh[j + L * No] * VhhhA[k + i * No + L * NoNo];
        }
        // PARTILCE DIAGRAMS: TAPHH and VABPH
        for (size_t E = 0; E < Nv; E++) {
          // t[aEij] * V[bcEk]        P0
          // t[aEik] * V[cbEj]        P3 // TODO: CHECK THIS ONE, I DONT KNOW
          Tijk[ijk] += TAphh[E + i * Nv + j * NoNv] * VBCph[E + k * Nv];
          Tijk[ijk] += TAphh[E + i * Nv + k * NoNv] * VCBph[E + j * Nv];

          // t[cEki] * V[abEj]        P5
          // t[cEkj] * V[baEi]        P2
          Tijk[ijk] += TCphh[E + k * Nv + i * NoNv] * VABph[E + j * Nv];
          Tijk[ijk] += TCphh[E + k * Nv + j * NoNv] * VBAph[E + i * Nv];

          // t[bEji] * V[acEk]        P1
          // t[bEjk] * V[caEi]        P4
          Tijk[ijk] += TBphh[E + j * Nv + i * NoNv] * VACph[E + k * Nv];
          Tijk[ijk] += TBphh[E + j * Nv + k * NoNv] * VCAph[E + i * Nv];
        }
      }
#endif /* defined(ATRIP_USE_DGEMM) */
}

// instantiate templates
template void doubles_contribution<double>(size_t const No,
                                           size_t const Nv
                                           // -- VABCI
                                           ,
                                           DataPtr<double> const VABph,
                                           DataPtr<double> const VACph,
                                           DataPtr<double> const VBCph,
                                           DataPtr<double> const VBAph,
                                           DataPtr<double> const VCAph,
                                           DataPtr<double> const VCBph
                                           // -- VHHHA
                                           ,
                                           DataPtr<double> const VhhhA,
                                           DataPtr<double> const VhhhB,
                                           DataPtr<double> const VhhhC
                                           // -- TA
                                           ,
                                           DataPtr<double> const TAphh,
                                           DataPtr<double> const TBphh,
                                           DataPtr<double> const TCphh
                                           // -- TABIJ
                                           ,
                                           DataPtr<double> const TABhh,
                                           DataPtr<double> const TAChh,
                                           DataPtr<double> const TBChh
                                           // -- TIJK
                                           ,
                                           DataFieldType<double> *Tijk,
                                           // -- tmp buffers
                                           DataFieldType<double> *_t_buffer,
                                           DataFieldType<double> *_vhhh);

template void doubles_contribution<Complex>(size_t const No,
                                            size_t const Nv
                                            // -- VABCI
                                            ,
                                            DataPtr<Complex> const VABph,
                                            DataPtr<Complex> const VACph,
                                            DataPtr<Complex> const VBCph,
                                            DataPtr<Complex> const VBAph,
                                            DataPtr<Complex> const VCAph,
                                            DataPtr<Complex> const VCBph
                                            // -- VHHHA
                                            ,
                                            DataPtr<Complex> const VhhhA,
                                            DataPtr<Complex> const VhhhB,
                                            DataPtr<Complex> const VhhhC
                                            // -- TA
                                            ,
                                            DataPtr<Complex> const TAphh,
                                            DataPtr<Complex> const TBphh,
                                            DataPtr<Complex> const TCphh
                                            // -- TABIJ
                                            ,
                                            DataPtr<Complex> const TABhh,
                                            DataPtr<Complex> const TAChh,
                                            DataPtr<Complex> const TBChh
                                            // -- TIJK
                                            ,
                                            DataFieldType<Complex> *Tijk,
                                            // -- tmp buffers
                                            DataFieldType<Complex> *_t_buffer,
                                            DataFieldType<Complex> *_vhhh);
// Doubles contribution:2 ends here

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:2]]
} // namespace atrip
// Epilog:2 ends here
