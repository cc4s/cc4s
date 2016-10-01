/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCALAPACK_DEFINED
#define SCALAPACK_DEFINED

#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
  int numroc_(
    int *n, int *nb, int *iproc, int *isrcproc, int *nprocs
  );
  void pdgesvd_(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const double *a, const int *ia, const int *ja, const int *desca,
    double *s, double *u, const int *iu, const int *ju, const int *descu,
    double *vt, const int *ivt, const int *jvt, const int *descvt,
    double *work, const int *lwork, double *rwork, int *info
  );
  void pzgesvd_(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const complex *a, const int *ia, const int *ja, const int *desca,
    double *s, complex *u, const int *iu, const int *ju, const int *descu,
    complex *vt, const int *ivt, const int *jvt, const int *descvt,
    complex *work, const int *lwork, double *rwork, int *info
  );
};

#endif

