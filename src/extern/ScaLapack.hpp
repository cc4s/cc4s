/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_DEFINED
#define SCA_LAPACK_DEFINED

#include <math/Complex.hpp>

// TODO: use define for name mangling: underscore or not

extern "C" {
  void descinit_(
    int *descriptor,
    const int *globalRows,  const int *globalColumns,
    const int *blockRows, const int *blockColumns,
    const int *rowOffset, const int *columnOffset,
    const int *context,
    const int *leadingDimensionSize,
    int *info
  );
  int numroc_(
    const int *globalSize, const int *blockSize,
    const int *iproc, const int *isrcproc, const int *nprocs
  );
  int indxl2g_(
    const int *localIndex, const int *blockSize,
    const int *iproc, const int *isrcproc, const int *nprocs
  );
  int indxg2l_(
    const int *globalIndex, const int *blockSize,
    const int *iproc, const int *isrcproc, const int *nprocs
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

  void pdgemm_(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const double *alpha,
    const double *a, const int *ia, const int *ja, const int *desca,
    const double *b, const int *ib, const int *jb, const int *descb,
    const double *beta,
    double *C, const int *ic, const int *jc, const int *descc
  );
  void pzgemm_(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const complex *alpha,
    const complex *a, const int *ia, const int *ja, const int *desca,
    const complex *b, const int *ib, const int *jb, const int *descb,
    const complex *beta,
    complex *C, const int *ic, const int *jc, const int *descc
  );
};

// overloaded functions wrappers for invocation from template methods
namespace cc4s {
  void pgesvd(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const double *a, const int *ia, const int *ja, const int *desca,
    double *s, double *u, const int *iu, const int *ju, const int *descu,
    double *vt, const int *ivt, const int *jvt, const int *descvt,
    double *work, const int *lwork, double *rwork, int *info
  );
  void pgesvd(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const complex *a, const int *ia, const int *ja, const int *desca,
    double *s, complex *u, const int *iu, const int *ju, const int *descu,
    complex *vt, const int *ivt, const int *jvt, const int *descvt,
    complex *work, const int *lwork, double *rwork, int *info
  );

  void pgemm(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const double *alpha,
    const double *a, const int *ia, const int *ja, const int *desca,
    const double *b, const int *ib, const int *jb, const int *descb,
    const double *beta,
    double *C, const int *ic, const int *jc, const int *descc
  );
  void pgemm(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const complex *alpha,
    const complex *a, const int *ia, const int *ja, const int *desca,
    const complex *b, const int *ib, const int *jb, const int *descb,
    const complex *beta,
    complex *C, const int *ic, const int *jc, const int *descc
  );
}

#endif

