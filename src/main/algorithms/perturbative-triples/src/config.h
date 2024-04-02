/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Enable the allocation of additional free buffers at runtime if the free
   buffer pool runs empty. */
/* #undef ATRIP_ALLOCATE_ADDITIONAL_FREE_POINTERS */

/* When using MPI, assume support for CUDA aware mpi by the given MPI
   implementation. */
/* #undef ATRIP_CUDA_AWARE_MPI */

/* Atrip debug level */
#define ATRIP_DEBUG 1

/* Disable the step of slicing tensors for CTF, this is useful for example for
   benchmarking or testing. */
/* #undef ATRIP_DONT_SLICE */

/* Enable mpi staging buffers int the GPU code. Every source will first be
   copied to this buffer and then subsequently sent. You should probably
   enable ATRIP_ALLOCATE_ADDITIONAL_FREE_POINTERS since the additional buffers
   are free pointers and they might run out, so dynamic allocation of free
   buffers is necessary. */
/* #undef ATRIP_MPI_STAGING_BUFFERS */

/* Run slow but correct code for the mapping of (iteration, rank) to tuple of
   the naive tuple distribution. */
/* #undef ATRIP_NAIVE_SLOW */

/* Run only the parts of atrip that involve dgemm calls, this is useful for
   benchmarking and testing the code, it is intended for developers of Atrip.
   */
/* #undef ATRIP_ONLY_DGEMM */

/* Print the tuple energies for each (abc) tuple at then end of the run. */
/* #undef ATRIP_PRINT_TUPLES */

/* When using CUDA, activate storing all sources (slices of the input tensors)
   in the GPU. This means that a lot of GPUs will be needed. */
/* #undef ATRIP_SOURCES_IN_GPU */

/* Use dgemm for the doubles equations */
#define ATRIP_USE_DGEMM 1

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Code with accelerator support */
/* #undef HAVE_ACC */

/* Define if you have a BLAS library. */
/* #undef HAVE_BLAS */

/* Use CTF library */
/* #undef HAVE_CTF */

/* AMD CUDA support flag */
/* #undef HAVE_CUDA */

/* define if the compiler supports basic C++11 syntax */
#define HAVE_CXX11 1

/* Wether we are using HIP */
/* #undef HAVE_HIP */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
/* #undef HAVE_LAPACK */

/* Define to 1 if you have the `omnitrace-user' library (-lomnitrace-user). */
/* #undef HAVE_LIBOMNITRACE_USER */

/* Define if you have the MPI library. */
#define HAVE_MPI 1

/* Define if you have a SCALAPACK library. */
/* #undef HAVE_SCALAPACK */

/* define if the compiler supports static_cast<> */
#define HAVE_STATIC_CAST /**/

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* define if the compiler supports basic templates */
#define HAVE_TEMPLATES /**/

/* define if the compiler recognizes typename */
#define HAVE_TYPENAME /**/

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Wether we are using HIP */
/* #undef HOP_TARGET_HIP */

/* Name of package */
#define PACKAGE "atrip"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://gihtub.com/alejandrogallo/atrip"

/* Define to the full name of this package. */
#define PACKAGE_NAME "atrip"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "atrip 0.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "atrip"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.0.0"

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "0.0.0"

/* Use profiler for HIP or CUDA */
/* #undef WITH_PROFILER */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
