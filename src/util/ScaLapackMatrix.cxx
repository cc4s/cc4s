#include <util/ScaLapackMatrix.hpp>

#include <extern/ScaLapack.hpp>
#include <math/Complex.hpp>
#include <util/Log.hpp>

using namespace cc4s;

ScaLapackDescriptor::ScaLapackDescriptor(
  BlacsWorld *blacsWorld,
  int lens_[2]
):
  dataType(1),
  blacsContext(blacsWorld->context)
{
  for (int i(0); i < 2; ++i) {
    lens[i] = lens_[i];
    blockSize[i] = std::max(lens[i] / blacsWorld->lens[i], 1);
    offset[i] = 0;
    localLens[i] = numroc_(
      &lens[i], &blockSize[i],
      &blacsWorld->firstElement[i], &offset[i], &blacsWorld->lens[i]
    );
/*
    LOG_RANK(1, "ScaLapackMatrix") << "rank=" << blacsWorld->rank << ", lens[" << i << "]=" << lens[i] << std::endl;
    LOG_RANK(1, "ScaLapackMatrix") << "rank=" << blacsWorld->rank << ", blockSize[" << i << "]=" << blockSize[i] << std::endl;
    LOG_RANK(1, "ScaLapackMatrix") << "rank=" << blacsWorld->rank << ", offset[" << i << "]=" << offset[i] << std::endl;
    LOG_RANK(1, "ScaLapackMatrix") << "rank=" << blacsWorld->rank << ", localLens[" << i << "]=" << localLens[i] << std::endl;
    LOG_RANK(1, "ScaLapackMatrix") << "rank=" << blacsWorld->rank << ", firstelement[" << i << "]=" << blacsWorld->firstElement[i] << std::endl;
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
ScaLapackMatrix<double>::ScaLapackMatrix(ScaLapackMatrix<double> &A);
template
ScaLapackMatrix<complex>::ScaLapackMatrix(ScaLapackMatrix<complex> &A);


template <typename F>
ScaLapackMatrix<F>::ScaLapackMatrix(
  CTF::Matrix<F> &A, BlacsWorld *blacsWorld_
):
  ScaLapackDescriptor(blacsWorld_, A.lens),
  blacsWorld(blacsWorld_)
{
  // allocate local data
  localValues = new F[localLens[0]*localLens[1]];
  localIndices = new int64_t[localLens[0]*localLens[1]];
  // determine global indices of local data in block cyclic distribution scheme
  int64_t index(0);
  for (
    int64_t j(blacsWorld->firstElement[1]);
    j < lens[1];
    j += blacsWorld->lens[1]
  ) {
    for (
      int64_t i(blacsWorld->firstElement[0]);
      i < lens[0];
      i += blacsWorld->lens[0]
    ) {
      localIndices[index++] = i + lens[0]*j;
    }
  }
  A.read(localLens[0]*localLens[1], localIndices, localValues);
}

// instantiate
template
ScaLapackMatrix<double>::ScaLapackMatrix(
  CTF::Matrix<double> &A, BlacsWorld *blacsWorld_
);
template
ScaLapackMatrix<complex>::ScaLapackMatrix(
  CTF::Matrix<complex> &A, BlacsWorld *blacsWorld_
);


template <typename F>
ScaLapackMatrix<F>::~ScaLapackMatrix() {
  delete[] localValues;
  delete[] localIndices;
}

// instantiate
template
ScaLapackMatrix<double>::~ScaLapackMatrix();
template
ScaLapackMatrix<complex>::~ScaLapackMatrix();

