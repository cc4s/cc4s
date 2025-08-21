#ifndef ACC_HPP_
#define ACC_HPP_

#include "config.h"
#include <atrip/Utils.hpp>

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// -- C U D A
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

#if defined(HAVE_CUDA)

#  include <atrip/CUDA.hpp>
#  define ACC_FUNCALL(fname, i, j, ...) fname<<<(i), (j)>>>(__VA_ARGS__)
#  define ACC_CHECK_SUCCESS(message, ...)                                      \
    do {                                                                       \
      CUresult result = __VA_ARGS__;                                           \
      if (result != CUDA_SUCCESS) {                                            \
        auto msg = _FORMAT("\t!!CUDA_ERROR(%d): %s:%d\v%s",                    \
                           result,                                             \
                           __FILE__,                                           \
                           __LINE__,                                           \
                           message);                                           \
        std::cerr << msg;                                                      \
        throw msg;                                                             \
      }                                                                        \
    } while (0)
#  define ACC_CHECK_BLAS(message, ...)                                         \
    do {                                                                       \
      cublasStatus_t result = __VA_ARGS__;                                     \
      if (result != 0) {                                                       \
        auto msg = _FORMAT("\t!!CUBLAS_ERROR(%d): %s:%d\v%s",                  \
                           result,                                             \
                           __FILE__,                                           \
                           __LINE__,                                           \
                           message);                                           \
        std::cerr << msg;                                                      \
        throw msg;                                                             \
      }                                                                        \
    } while (0)
#  define ACC_INIT(arg) cuInit(arg)
#  define ACC_DEVICE_GET_COUNT(ptr) cuDeviceGetCount(ptr)
#  define ACC_DEVICE CUdevice
#  define ACC_CONTEXT CUcontext
#  define ACC_DEVICE_GET(ptr, i) cuDeviceGet(ptr, i)
#  define ACC_CONTEXT_CREATE(ptr, i, dev) cuCtxCreate(ptr, i, dev)
#  define ACC_CONTEXT_SET_CURRENT(ctx) cuCtxSetCurrent(ctx)
#  define ACC_DEVICE_GET_PROPERTIES(ptr, dev) cuDeviceGetProperties(ptr, dev)
#  define ACC_BLAS_CREATE(ptr) cublasCreate(ptr)
#  define ACC_BLAS_STATUS cublasStatus_t
#  define ACC_BLAS_HANDLE cublasHandle_t
#  define ACC_BLAS_SGEMM cublasSgemm
#  define ACC_BLAS_DGEMM cublasDgemm
#  define ACC_BLAS_ZGEMM cublasZgemm
#  define ACC_BLAS_SCOPY cublasScopy
#  define ACC_BLAS_DCOPY cublasDcopy
#  define ACC_BLAS_ZCOPY cublasZcopy
#  define ACC_BLAS_OP cublasOperation_t
#  define ACC_BLAS_OP_C CUBLAS_OP_C
#  define ACC_BLAS_OP_N CUBLAS_OP_N
#  define ACC_BLAS_OP_T CUBLAS_OP_T
#  define ACC_DEVICE_SYNCHRONIZE() cuCtxSynchronize()
#  define ACC_MEMCPY_HOST_TO_DEV(devptr, hostptr, size)                        \
    cuMemcpyHtoD((devptr), (hostptr), (size))
#  define ACC_MEMCPY_DEV_TO_HOST(hostptr, devptr, size)                        \
    cuMemcpyDtoH((hostptr), (devptr), size)
// types

#  define ACC_DOUBLE_COMPLEX cuDoubleComplex
#  define ACC_BLAS_DOUBLE_COMPLEX cuDoubleComplex
#  define ACC_DEVICE_PTR CUdeviceptr
// complex
#  define ACC_COMPLEX_MUL cuCmul
#  define ACC_COMPLEX_DIV cuCdiv
#  define ACC_COMPLEX_REAL cuCreal
#  define ACC_COMPLEX_SUB cuCsub
#  define ACC_COMPLEX_ADD cuCadd
// memset
#  define ACC_MEM_SET_D32 cuMemsetD32_v2
#  define ACC_FREE cuMemFree
#  define ACC_MEMCPY cuMemcpy
#  define ACC_MEM_GET_INFO cuMemGetInfo
#  define ACC_GET_DEVICE_NAME cuDeviceGetName
#  define ACC_DEVICE_TOTAL_MEM cuDeviceTotalMem
#  define ACC_MEM_ALLOC cuMemAlloc

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// -- H I P
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

#elif defined(HAVE_HIP)

#  include <atrip/HIP.hpp>
#  define ACC_FUNCALL(fname, i, j, ...) fname<<<(i), (j)>>>(__VA_ARGS__)
#  define ACC_CHECK_SUCCESS(message, ...)                                      \
    do {                                                                       \
      hipError_t result = __VA_ARGS__;                                         \
      if (result != hipSuccess) {                                              \
        auto msg = _FORMAT("\t!!HIP_ERROR(%d): %s:%d\v%s",                     \
                           result,                                             \
                           __FILE__,                                           \
                           __LINE__,                                           \
                           message);                                           \
        std::cerr << msg;                                                      \
        throw msg;                                                             \
      }                                                                        \
    } while (0)
#  define ACC_CHECK_BLAS(message, ...)                                         \
    do {                                                                       \
      hipblasStatus_t result = __VA_ARGS__;                                    \
      if (result != 0) {                                                       \
        auto msg = _FORMAT("\t!!HIPBLAS_ERROR(%d): %s:%d\v%s",                 \
                           result,                                             \
                           __FILE__,                                           \
                           __LINE__,                                           \
                           message);                                           \
        std::cerr << msg;                                                      \
        throw msg;                                                             \
      }                                                                        \
    } while (0)
#  define ACC_INIT(arg) hipInit(arg)
#  define ACC_DEVICE_GET_COUNT(ptr) hipGetDeviceCount(ptr)
#  define ACC_DEVICE hipDevice_t
#  define ACC_CONTEXT hipCtx_t
#  define ACC_DEVICE_GET(ptr, i) hipDeviceGet(ptr, i)
#  define ACC_CONTEXT_CREATE(ptr, i, dev) hipCtxCreate(ptr, i, dev)
#  define ACC_CONTEXT_SET_CURRENT(ctx) hipCtxSetCurrent(ctx)
#  define ACC_DEVICE_GET_PROPERTIES(ptr, dev) hipDeviceGetProperties(ptr, dev)
#  define ACC_BLAS_CREATE(ptr) hipblasCreate(ptr)
#  define ACC_BLAS_STATUS hipblasStatus_t
#  define ACC_BLAS_HANDLE hipblasHandle_t
#  define ACC_BLAS_SGEMM hipblasSgemm
#  define ACC_BLAS_DGEMM hipblasDgemm
#  define ACC_BLAS_ZGEMM hipblasZgemm
#  define ACC_BLAS_SCOPY hipblasScopy
#  define ACC_BLAS_DCOPY hipblasDcopy
#  define ACC_BLAS_ZCOPY hipblasZcopy
#  define ACC_BLAS_OP hipblasOperation_t
#  define ACC_BLAS_OP_C HIPBLAS_OP_C
#  define ACC_BLAS_OP_N HIPBLAS_OP_N
#  define ACC_BLAS_OP_T HIPBLAS_OP_T
#  define ACC_DEVICE_SYNCHRONIZE() hipDeviceSynchronize()
#  define ACC_MEMCPY_HOST_TO_DEV(devptr, hostptr, size)                        \
    hipMemcpyHtoD((devptr), (hostptr), (size))
#  define ACC_MEMCPY_DEV_TO_HOST(hostptr, devptr, size)                        \
    hipMemcpyDtoH((hostptr), (devptr), (size))
// types
#  define ACC_DOUBLE_COMPLEX hipDoubleComplex
#  define ACC_BLAS_DOUBLE_COMPLEX hipDoubleComplex
#  define ACC_DEVICE_PTR hipDeviceptr_t
// complex
#  define ACC_COMPLEX_MUL hipCmul
#  define ACC_COMPLEX_DIV hipCdiv
#  define ACC_COMPLEX_REAL hipCreal
#  define ACC_COMPLEX_SUB hipCsub
#  define ACC_COMPLEX_ADD hipCadd
// memset
#  define ACC_MEM_SET_D32 hipMemsetD32
#  define ACC_FREE hipFree
#  define ACC_MEMCPY hipMemcpyDtoD
#  define ACC_MEM_GET_INFO hipMemGetInfo
#  define ACC_GET_DEVICE_NAME hipDeviceGetName
#  define ACC_DEVICE_TOTAL_MEM hipDeviceTotalMem
#  define ACC_MEM_ALLOC hipMalloc

#else

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// -- NONE
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

#  define ACC_FUNCALL(fname, i, j, ...) fname(__VA_ARGS__)
#  define ACC_DEVICE_SYNCHRONIZE()                                             \
    do {                                                                       \
    } while (0)
#  define ACC_FREE std::free

#  define ACC_BLAS_OP (char *)
#  define ACC_BLAS_OP_C "C"
#  define ACC_BLAS_OP_N "N"
#  define ACC_BLAS_OP_T "T"

#  define __MAYBE_GLOBAL__
#  define __MAYBE_DEVICE__
#  define __MAYBE_HOST__
#  define __INLINE__ inline

#endif

#endif
