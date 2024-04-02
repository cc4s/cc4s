#pragma once

#include <atrip/Utils.hpp>

#if defined(HAVE_CUDA)
#  include <cuda.h>
#endif

#include <sstream>

#if defined(HAVE_CUDA) && defined(__CUDACC__)
#  define __MAYBE_GLOBAL__ __global__
#  define __MAYBE_DEVICE__ __device__
#  define __MAYBE_HOST__ __host__
#  define __INLINE__ __inline__
#else
#  define __MAYBE_GLOBAL__
#  define __MAYBE_DEVICE__
#  define __MAYBE_HOST__
#  define __INLINE__ inline
#endif

#if defined(HAVE_CUDA)
#  define ACC_FUNCALL(fname, i, j, ...) fname<<<(i), (j)>>>(__VA_ARGS__)
#else
#  define ACC_FUNCALL(fname, i, j, ...) fname(__VA_ARGS__)
#endif /*  defined(HAVE_CUDA) */

#define _CHECK_CUDA_SUCCESS(message, ...)                                      \
  do {                                                                         \
    CUresult result = __VA_ARGS__;                                             \
    if (result != CUDA_SUCCESS) {                                              \
      auto msg = _FORMAT("\t!!CUDA_ERROR(%d): %s:%d\v%s",                      \
                         result,                                               \
                         __FILE__,                                             \
                         __LINE__,                                             \
                         message);                                             \
      std::cerr << msg;                                                        \
      throw msg;                                                               \
    }                                                                          \
  } while (0)

#define _CHECK_CUBLAS_SUCCESS(message, ...)                                    \
  do {                                                                         \
    cublasStatus_t result = __VA_ARGS__;                                       \
    if (result != 0) {                                                         \
      auto msg = _FORMAT("\t!!CUBLAS_ERROR(%d): %s:%d\v%s",                    \
                         result,                                               \
                         __FILE__,                                             \
                         __LINE__,                                             \
                         message);                                             \
      std::cerr << msg;                                                        \
      throw msg;                                                               \
    }                                                                          \
  } while (0)

#define _CUDA_MALLOC(msg, ptr, __size)                                         \
  do {                                                                         \
    WITH_CHRONO("malloc:cuda",                                                 \
                const CUresult error = cuMemAlloc((ptr), (__size));)           \
    if (*(ptr) == 0UL) { throw("UNSUFICCIENT GPU MEMORY for " msg); }          \
    if (error != CUDA_SUCCESS) {                                               \
      std::stringstream s;                                                     \
      s << "CUDA Error allocating memory for " << (msg) << "code " << error    \
        << "\n";                                                               \
      throw s.str();                                                           \
    }                                                                          \
  } while (0)

#define _CUDA_FREE(msg, ptr) _CHECK_CUDA_SUCCESS(msg, cuMemFree(ptr))
