#pragma once

#include <atrip/Utils.hpp>

#if defined(HAVE_CUDA)
#  include <cuda.h>
#  define CUBLASAPI
#  include <cublas_api.h>
#  include <cublas_v2.h>
#  include <nvToolsExt.h>
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

#define _CUDA_MALLOC(msg, ptr, __size)                                         \
  do {                                                                         \
    WITH_CHRONO("malloc:device",                                               \
                const CUresult error = cuMemAlloc((ptr), (__size));)           \
    if (*(ptr) == 0UL) { throw("UNSUFICCIENT GPU MEMORY for " msg); }          \
    if (error != CUDA_SUCCESS) {                                               \
      std::stringstream s;                                                     \
      s << "CUDA Error allocating memory for " << (msg) << "code " << error    \
        << "\n";                                                               \
      throw s.str();                                                           \
    }                                                                          \
  } while (0)

#define _CUDA_FREE(msg, ptr) ACC_CHECK_SUCCESS(msg, cuMemFree(ptr))
