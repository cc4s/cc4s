/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_MATRIX_DEFINED
#define SCA_LAPACK_MATRIX_DEFINED

#include <util/BlacsWorld.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ScaLapackDescriptor {
  public:
    ScaLapackDescriptor(
      BlacsWorld *blacsWorld,
      int lens_[2]
    );
    int dataType;     // 1 for dense matrices
    int blacsContext;
    int lens[2];
    int blockSize[2];
    int offset[2];
    int localLens[2]; // NOTE: assumes column major format
    // NOTE: the last element localLens[1]=number of local columns
    // is actually not part of the scalapack descriptor
  };
  
  template <typename F>
  class ScaLapackMatrix: public ScaLapackDescriptor {
  public:
    // copy other ScaLapack matrix
    ScaLapackMatrix(ScaLapackMatrix<F> &A);
    // construct from CTF::Matrix
    ScaLapackMatrix(CTF::Matrix<F> &A, BlacsWorld *blacsWorld);
    ~ScaLapackMatrix();
    const int *getDescriptor() {
      return &dataType;
    }

    BlacsWorld *blacsWorld;
    int64_t *localIndices;
    F *localValues;
  };
}

#endif

