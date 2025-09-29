#ifndef HIP_HPP_
#define HIP_HPP_

#ifndef __HIP_PLATFORM_AMD__
#  define __HIP_PLATFORM_AMD__
#endif

#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>

#include <omnitrace/categories.h>
#include <omnitrace/types.h>
#include <omnitrace/user.h>

#define HIPBLAS_V2
#define ROCM_MATHLIBS_API_USE_HIP_COMPLEX
#include <hipblas.h>

#if defined(HAVE_HIP) && defined(__HIPCC__)
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

#define _HIP_MALLOC(msg, ptr, __size)                                          \
  do {                                                                         \
    WITH_CHRONO("malloc:device",                                               \
                const hipError_t error = hipMalloc((ptr), (__size));)          \
    if (*(ptr) == 0UL) { throw("UNSUFICCIENT GPU MEMORY for " msg); }          \
    if (error != hipSuccess) {                                                 \
      std::stringstream s;                                                     \
      s << "HIP Error allocating memory for " << (msg) << "code " << error     \
        << "\n";                                                               \
      throw s.str();                                                           \
    }                                                                          \
  } while (0)

#define _HIP_FREE(msg, ptr) ACC_CHECK_SUCCESS(msg, hipFree(ptr))

#endif
