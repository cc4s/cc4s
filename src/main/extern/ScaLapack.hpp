/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SCA_LAPACK_DEFINED
#define SCA_LAPACK_DEFINED

#include <math/Real.hpp>
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

  void pdsyev_(
    const char *jobz, const char *upperLower,
    const int *m,
    const cc4s::Real<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *lambda,
    cc4s::Real<64> *z, const int *iz, const int *jz, const int *descz,
    cc4s::Real<64> *work, const int *lwork,
    int *info
  );
  void pzheev_(
    const char *jobz, const char *upperLower,
    const int *m,
    const cc4s::Complex<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *lambda,
    cc4s::Complex<64> *z, const int *iz, const int *jz, const int *descz,
    cc4s::Complex<64> *work, const int *lwork,
    cc4s::Real<64> *realWork, const int *lRealCount,
    int *info
  );

  void pdsyevd_(
    const char *jobz, const char *upperLower,
    const int *m,
    const cc4s::Real<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *lambda,
    cc4s::Real<64> *z, const int *iz, const int *jz, const int *descz,
    cc4s::Real<64> *work, const int *lwork,
    int *iwork, const int *liwork,
    int *info
  );
  void pzheevd_(
    const char *jobz, const char *upperLower,
    const int *m,
    const cc4s::Complex<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *lambda,
    cc4s::Complex<64> *z, const int *iz, const int *jz, const int *descz,
    cc4s::Complex<64> *work, const int *lwork,
    cc4s::Real<64> *realWork, const int *lRealCount,
    int *iwork, const int *liwork,
    int *info
  );

  void pdgesvd_(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const cc4s::Real<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *s, cc4s::Real<64> *u,
    const int *iu, const int *ju, const int *descu,
    cc4s::Real<64> *vt, const int *ivt, const int *jvt, const int *descvt,
    cc4s::Real<64> *work, const int *lwork,
    int *info
  );
  void pzgesvd_(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const cc4s::Complex<64> *a, const int *ia, const int *ja, const int *desca,
    cc4s::Real<64> *s,
    cc4s::Complex<64> *u, const int *iu, const int *ju, const int *descu,
    cc4s::Complex<64> *vt, const int *ivt, const int *jvt, const int *descvt,
    cc4s::Complex<64> *work, const int *lwork, cc4s::Real<64> *rwork,
    int *info
  );

  void pdgemm_(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const cc4s::Real<64> *alpha,
    const cc4s::Real<64> *a, const int *ia, const int *ja, const int *desca,
    const cc4s::Real<64> *b, const int *ib, const int *jb, const int *descb,
    const cc4s::Real<64> *beta,
    cc4s::Real<64> *C, const int *ic, const int *jc, const int *descc
  );
  void pzgemm_(
    const char *opA, const char *opB,
    const int *m, const int *n, const int *k,
    const cc4s::Complex<64> *alpha,
    const cc4s::Complex<64> *a, const int *ia, const int *ja, const int *desca,
    const cc4s::Complex<64> *b, const int *ib, const int *jb, const int *descb,
    const cc4s::Complex<64> *beta,
    cc4s::Complex<64> *C, const int *ic, const int *jc, const int *descc
  );
}

#endif

