/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <math/Complex.hpp>
#include <tcc/DryTensor.hpp>

#include <array>

using namespace cc4s;

int64_t DryMemory::currentTotalSize = 0, DryMemory::maxTotalSize = 0;
std::vector<DryMemory::ExtendingResource> DryMemory::extendingResources;


template <typename F>
DryMatrix<F>::DryMatrix(
  int rows, int cols, int sym, SourceLocation const &location
):
  DryTensor<F>(
    2, std::array<int,2>{{rows, cols}}.data(),
    std::array<int,2>{{sym, 0}}.data(),
    location
  )
{
}
// instantiate
template
DryMatrix<Float64>::DryMatrix(
  int rows, int cols, int sym, SourceLocation const &location
);
template
DryMatrix<Complex64>::DryMatrix(
  int rows, int cols, int sym, SourceLocation const &location
);

template <typename F>
DryVector<F>::DryVector(
  int elements, SourceLocation const &location
):
  DryTensor<F>(
    1, std::array<int,1>{{elements}}.data(),
    std::array<int,1>{{0}}.data(),
    location
  )
{
}
// instantiate
template
DryVector<Float64>::DryVector(int elements, SourceLocation const &location);
template
DryVector<Complex64>::DryVector(int elements, SourceLocation const &location);


template <typename F>
DryScalar<F>::DryScalar(
  SourceLocation const &location
):
  DryTensor<F>(
    0, std::array<int,0>{{}}.data(),
    std::array<int,0>{{}}.data(),
    location
  )
{
}
// instantiate
template
DryScalar<Float64>::DryScalar(SourceLocation const &location);
template
DryScalar<Complex64>::DryScalar(SourceLocation const &location);

template <typename F>
DryScalar<F>::DryScalar(
  F const value, // will be discarded
  SourceLocation const &location
):
  DryTensor<F>(
    0, std::array<int,0>{{}}.data(),
    std::array<int,0>{{}}.data(),
    location
  )
{
}
// instantiate
template
DryScalar<Float64>::DryScalar(
  const Float64 value, SourceLocation const &location
);
template
DryScalar<Complex64>::DryScalar(
  const Complex64 value, SourceLocation const &location
);
template
DryScalar<Float128>::DryScalar(
  const Float128 value, SourceLocation const &location
);
template
DryScalar<Complex128>::DryScalar(
  const Complex128 value, SourceLocation const &location
);
