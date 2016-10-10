#include <extern/ScaLapack.hpp>

void cc4s::pheev(
  const char *jobz, const char *upperLower,
  const int *m,
  const double *a, const int *ia, const int *ja, const int *desca,
  double *lambda,
  double *z, const int *iz, const int *jz, const int *descz,
  int *info
) {
  double optimalWork;
  int workCount(-1);
  pdsyev_(
    jobz, upperLower,
    m,
    a, ia, ja, desca,
    lambda,
    z, iz, jz, descz,
    &optimalWork, &workCount, info
  );
  workCount = static_cast<int>(optimalWork+0.5);
  double *work(new double[workCount]);
  pdsyev_(
    jobz, upperLower,
    m,
    a, ia, ja, desca,
    lambda,
    z, iz, jz, descz,
    work, &workCount, info
  );
  delete[] work;
}

void cc4s::pheev(
  const char *jobz, const char *upperLower,
  const int *m,
  const complex *a, const int *ia, const int *ja, const int *desca,
  double *lambda,
  complex *z, const int *iz, const int *jz, const int *descz,
  int *info
) {
  complex optimalWork;
  double optimalRealWork;
  int workCount(-1), realWorkCount(-1);
  pzheev_(
    jobz, upperLower,
    m,
    a, ia, ja, desca,
    lambda,
    z, iz, jz, descz,
    &optimalWork, &workCount, &optimalRealWork, &realWorkCount,
    info
  );
  workCount = static_cast<int>(std::real(optimalWork)+0.5);
  complex *work(new complex[workCount]);
  double *realWork(new double[static_cast<int64_t>(optimalRealWork+0.5)]);
  pzheev_(
    jobz, upperLower,
    m,
    a, ia, ja, desca,
    lambda,
    z, iz, jz, descz,
    work, &workCount, realWork, &realWorkCount,
    info
  );
  delete[] work; delete[] realWork;
}


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

