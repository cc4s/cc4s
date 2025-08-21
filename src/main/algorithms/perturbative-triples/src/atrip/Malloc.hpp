#include <type_traits>

#if defined(HAVE_CUDA)
#  include <atrip/CUDA.hpp>
#elif defined(HAVE_HIP)
#  include <atrip/HIP.hpp>
#endif

// This can be used also in ATRIP_DRY
#if defined(HAVE_CUDA)
#  define FREE_DATA_PTR(msg, ptr) _CUDA_FREE(msg, (CUdeviceptr)ptr)
#else
#  define FREE_DATA_PTR(msg, ptr) std::free(ptr)
#endif /* defined(HAVE_CUDA) */

#if defined(ATRIP_DRY)

#  define MALLOC_HOST_DATA(msg, ptr, size) *(void **)ptr = (void *)malloc(16)

#  if defined(HAVE_CUDA)
#    define MALLOC_DATA_PTR(msg, ptr, _size)                                   \
      _CUDA_MALLOC(msg, (CUdeviceptr *)ptr, 16)
#  elif defined(HAVE_HIP)
#    define MALLOC_DATA_PTR(msg, ptr, _size) _HIP_MALLOC(msg, ptr, 16)
#  else
#    define MALLOC_DATA_PTR(msg, ptr, _size) MALLOC_HOST_DATA(msg, ptr, 16)
#  endif /* defined(HAVE_CUDA) */

#else

#  define MALLOC_HOST_DATA(msg, ptr, size) *(void **)ptr = (void *)malloc(size)

#  if defined(HAVE_CUDA)
#    define MALLOC_DATA_PTR(msg, ptr, _size)                                   \
      _CUDA_MALLOC(msg, (CUdeviceptr *)ptr, _size)
#  elif defined(HAVE_HIP)
#    define MALLOC_DATA_PTR(msg, ptr, _size) _HIP_MALLOC(msg, ptr, _size)
#  else
#    define MALLOC_DATA_PTR(msg, ptr, _size) MALLOC_HOST_DATA(msg, ptr, _size)
#  endif /* defined(HAVE_CUDA) */

#endif /* defined(ATRIP_DRY) */
