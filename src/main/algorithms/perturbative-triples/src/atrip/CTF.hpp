#ifndef ATRIP_CTF_HPP_
#define ATRIP_CTF_HPP_
#if defined(HAVE_CTF)

#  define WITH_CTF(...) __VA_ARGS__

#  if defined(__NVCC__)
#    pragma nv_diagnostic_push
#    if defined __NVCC_DIAG_PRAGMA_SUPPORT__
// http://www.ssl.berkeley.edu/~jimm/grizzly_docs/SSL/opt/intel/cc/9.0/lib/locale/en_US/mcpcom.msg
#      pragma nv_diag_suppress partial_override
#    else
#      pragma diag_suppress partial_override
#    endif
#    include <ctf.hpp>
#    pragma nv_diagnostic_pop
#  else
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wvla"
#    pragma GCC diagnostic ignored "-Wnonnull"
#    pragma GCC diagnostic ignored "-Wall"
#    pragma GCC diagnostic ignored "-Wint-in-bool-context"
#    pragma GCC diagnostic ignored "-Wunused-parameter"
#    pragma GCC diagnostic ignored "-Wdeprecated-copy"
#    pragma GCC diagnostic ignored "-Woverloaded-virtual"
#    include <ctf.hpp>
#    pragma GCC diagnostic pop
#  endif /* defined(__NVCC__) */
#else
#  define WITH_CTF(...)
#endif /* defined(HAVE_CTF) */
#endif
