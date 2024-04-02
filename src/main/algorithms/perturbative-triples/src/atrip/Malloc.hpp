#include <atrip/CUDA.hpp>
#include <type_traits>

// This can be used also in ATRIP_DRY
#if defined(HAVE_CUDA)
#  define FREE_DATA_PTR(msg, ptr) _CUDA_FREE(msg, ptr)
#else
#  define FREE_DATA_PTR(msg, ptr) std::free(ptr)
#endif /* defined(HAVE_CUDA) */

#if defined(ATRIP_DRY)

#  define MALLOC_HOST_DATA(msg, ptr, size) *(void **)ptr = (void *)malloc(16)

#  if defined(HAVE_CUDA)
#    define MALLOC_DATA_PTR(msg, ptr, _size) _CUDA_MALLOC(msg, ptr, 16)
#  else
#    define MALLOC_DATA_PTR(msg, ptr, _size) MALLOC_HOST_DATA(msg, ptr, 16)
#  endif /* defined(HAVE_CUDA) */

#else

#  define MALLOC_HOST_DATA(msg, ptr, size) *(void **)ptr = (void *)malloc(size)

#  if defined(HAVE_CUDA)
#    define MALLOC_DATA_PTR(msg, ptr, _size) _CUDA_MALLOC(msg, ptr, _size)
#  else
#    define MALLOC_DATA_PTR(msg, ptr, _size) MALLOC_HOST_DATA(msg, ptr, _size)
#  endif /* defined(HAVE_CUDA) */

#endif /* defined(ATRIP_DRY) */
