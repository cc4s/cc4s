#include <extern/ScaLapack.hpp>

void cc4s::pgemm(
  const char *opA, const char *opB,
  const int *m, const int *n, const int *k,
  const double *alpha,
  const double *a, const int *ia, const int *ja, const int *desca,
  const double *b, const int *ib, const int *jb, const int *descb,
  const double *beta,
  double *C, const int *ic, const int *jc, const int *descc
) {
  pdgemm_(
    opA, opB,
    m, n, k,
    alpha,
    a, ia, ja, desca,
    b, ib, jb, descb,
    beta,
    C, ic, jc, descc
  );
}

void cc4s::pgemm(
  const char *opA, const char *opB,
  const int *m, const int *n, const int *k,
  const complex *alpha,
  const complex *a, const int *ia, const int *ja, const int *desca,
  const complex *b, const int *ib, const int *jb, const int *descb,
  const complex *beta,
  complex *C, const int *ic, const int *jc, const int *descc
) {
  pzgemm_(
    opA, opB,
    m, n, k,
    alpha,
    a, ia, ja, desca,
    b, ib, jb, descb,
    beta,
    C, ic, jc, descc
  );
}

