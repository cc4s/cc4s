/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_MATRIX_DEFINED
#define SCA_LAPACK_MATRIX_DEFINED

#include <util/BlacsWorld.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ScaLapackDescriptor {
  public:
    ScaLapackDescriptor(BlacsWorld *blacsWorld, int lens_[2], int blockSize);
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
    /**
     * \brief Copies the content ScaLapack matrix
     * inheriting its distribution and its BlacsWorld.
     */
    ScaLapackMatrix(ScaLapackMatrix<F> &A);
    /**
     * \brief Construct a ScaLapack matrix from the content of a
     * CTF::Matrix in the given BlacsWorld.
     */
    ScaLapackMatrix(
      CTF::Matrix<F> &A, BlacsWorld *blacsWorld, int blockSize = 64
    );
    /**
     * \brief Frees all resources associated with the ScaLapack matrix on
     * all processes.
     */
    ~ScaLapackMatrix();

    /**
     * \brief Returns a pointer to a ScaLapack matrix descriptor, i.e.
     * the first 9 integers of the class ScaLapackDescriptor,
     * which are in the format ScaLapack routines expect it.
     */
    const int *getDescriptor() const {
      return &dataType;
    }

    /**
     * \brief Returns the pointer to the block data owned by this process.
     */
    const F *getLocalValues() const {
      return localValues;
    }

    /**
     * \brief Returns the pointer to the mutable block data owned by this
     * process.
     */
    F *getLocalValues() {
      return localValues;
    }

    /**
     * \brief Writes the content of the ScaLapackMatrix to the given
     * CTF::Matrix. The CTF::Matrix must be allocated in the correct shape.
     */
    void write(CTF::Matrix<F> &A);

  protected:
    /**
     * \brief Retrieves the global row or column index within the full matrix
     * from the given row or column index within the local block.
     */
    int getGlobalIndex(int localIndex, int dimension);

    /**
     * \brief The BlacsWorld specifying the processor grid this matrix
     * is distributed over.
     */
    BlacsWorld *blacsWorld;
    /**
     * \brief The global indices within the local matrix elements owned
     * by this process. The global index of an element in row i and
     * in column j of the full matrix is (i + lens[0]*j).
     */
    int64_t *localIndices;

    /**
     * \brief The block data owned by this process.
     */
    F *localValues;
  };
}

#endif

