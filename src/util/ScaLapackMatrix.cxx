#include <util/ScaLapackMatrix.hpp>

#include <extern/ScaLapack.hpp>
#include <math/Complex.hpp>
#include <util/Log.hpp>

using namespace cc4s;

ScaLapackDescriptor::ScaLapackDescriptor(
  BlacsWorld *blacsWorld,
  int lens_[2],
  int blockSize_
):
  dataType(1),
  blacsContext(blacsWorld->context)
{
  for (int d(0); d < 2; ++d) {
    lens[d] = lens_[d];
    blockSize[d] = std::min(blockSize_, lens_[d]);
    offset[d] = 0;
    localLens[d] = numroc_(
      &lens[d], &blockSize[d],
      &blacsWorld->firstElement[d], &offset[d], &blacsWorld->lens[d]
    );
/*
    std::cout << "rank=" << blacsWorld->rank << ", lens[" << d << "]=" << lens[d] << std::endl;
    std::cout << "rank=" << blacsWorld->rank << ", blockSize[" << d << "]=" << blockSize[d] << std::endl;
    std::cout << "rank=" << blacsWorld->rank << ", offset[" << d << "]=" << offset[d] << std::endl;
    std::cout << "rank=" << blacsWorld->rank << ", localLens[" << d << "]=" << localLens[d] << std::endl;
    std::cout << "rank=" << blacsWorld->rank << ", firstelement[" << d << "]=" << blacsWorld->firstElement[d] << std::endl;
*/
  }
}

template <typename F>
ScaLapackMatrix<F>::ScaLapackMatrix(
  ScaLapackMatrix<F> &A
):
  ScaLapackDescriptor(A),
  blacsWorld(A.blacsWorld),
  localIndices(new int64_t[localLens[0]*localLens[1]]),
  localValues(new F[localLens[0]*localLens[1]])
{
  for (int64_t i(0); i < static_cast<int64_t>(localLens[0])*localLens[1]; ++i) {
    localIndices[i] = A.localIndices[i];
    localValues[i] = A.localValues[i];
  }
}

// instantiate
template
ScaLapackMatrix<cc4s::Float64>::ScaLapackMatrix(ScaLapackMatrix<cc4s::Float64> &A);
template
ScaLapackMatrix<cc4s::Complex64>::ScaLapackMatrix(ScaLapackMatrix<cc4s::Complex64> &A);


template <typename F>
ScaLapackMatrix<F>::ScaLapackMatrix(
  CTF::Matrix<F> &A, BlacsWorld *blacsWorld_, int blockSize
):
  ScaLapackDescriptor(blacsWorld_, A.lens, blockSize),
  blacsWorld(blacsWorld_)
{
  // allocate local data
  localValues = new F[localLens[0]*localLens[1]];
  localIndices = new int64_t[localLens[0]*localLens[1]];
  // determine global indices of local data in block cyclic distribution scheme
  int64_t index(0);
  for (int localJ(0); localJ < localLens[1]; ++localJ) {
    int64_t j(getGlobalIndex(localJ, 1));
    for (int localI(0); localI < localLens[0]; ++localI) {
      int64_t i(getGlobalIndex(localI, 0));
      localIndices[index++] = i + lens[0]*j;
    }
  }
  A.read(localLens[0]*localLens[1], localIndices, localValues);
}

// instantiate
template
ScaLapackMatrix<cc4s::Float64>::ScaLapackMatrix(
  CTF::Matrix<cc4s::Float64> &A, BlacsWorld *blacsWorld, int blockSize
);
template
ScaLapackMatrix<cc4s::Complex64>::ScaLapackMatrix(
  CTF::Matrix<Complex64> &A, BlacsWorld *blacsWorld, int blockSize
);


template <typename F>
ScaLapackMatrix<F>::ScaLapackMatrix(
  CTF::Tensor<F> &A, int lens_[2], BlacsWorld *blacsWorld_, int blockSize
):
  ScaLapackDescriptor(blacsWorld_, lens_, blockSize),
  blacsWorld(blacsWorld_)
{
  // allocate local data
  localValues = new F[localLens[0]*localLens[1]];
  localIndices = new int64_t[localLens[0]*localLens[1]];
  // determine global indices of local data in block cyclic distribution scheme
  int64_t index(0);
  for (int localJ(0); localJ < localLens[1]; ++localJ) {
    int64_t j(getGlobalIndex(localJ, 1));
    for (int localI(0); localI < localLens[0]; ++localI) {
      int64_t i(getGlobalIndex(localI, 0));
      localIndices[index++] = i + lens[0]*j;
    }
  }
  A.read(index, localIndices, localValues);
}

// instantiate
template
ScaLapackMatrix<Float64>::ScaLapackMatrix(
  CTF::Tensor<Float64> &A, int lens_[2], BlacsWorld *blacsWorld, int blockSize
);
template
ScaLapackMatrix<Complex64>::ScaLapackMatrix(
  CTF::Tensor<Complex64> &A, int lens_[2], BlacsWorld *blacsWorld, int blockSize
);


template <typename F>
ScaLapackMatrix<F>::~ScaLapackMatrix() {
  delete[] localValues;
  delete[] localIndices;
}

// instantiate
template
ScaLapackMatrix<Float64>::~ScaLapackMatrix();
template
ScaLapackMatrix<Complex64>::~ScaLapackMatrix();


template <typename F>
void ScaLapackMatrix<F>::write(CTF::Matrix<F> &A) {
  // wait for all processes to finish pending operations
  blacsWorld->barrier();
  A.write(localLens[0]*localLens[1], localIndices, localValues);
}

template <typename F>
void ScaLapackMatrix<F>::write(CTF::Tensor<F> &A) {
  // wait for all processes to finish pending operations
  blacsWorld->barrier();
  A.write(localLens[0]*localLens[1], localIndices, localValues);
}

// instantiate
template
void ScaLapackMatrix<Float64>::write(CTF::Matrix<Float64> &A);
template
void ScaLapackMatrix<Complex64>::write(CTF::Matrix<Complex64> &A);

template
void ScaLapackMatrix<Float64>::write(CTF::Tensor<Float64> &A);
template
void ScaLapackMatrix<Complex64>::write(CTF::Tensor<Complex64> &A);


template <typename F>
int ScaLapackMatrix<F>::getGlobalIndex(int localIndex, int d) {
  // convert to and from Fortran indices
  ++localIndex;
  return indxl2g_(
    &localIndex, &blockSize[d],
    &blacsWorld->firstElement[d], &offset[d], &blacsWorld->lens[d]
  ) - 1;
}

