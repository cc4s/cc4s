/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <util/DryTensor.hpp>

#include <math/Complex.hpp>
#include <ctf.hpp>
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
    std::array<int,2>{{sym, NS}}.data(),
    location
  )
{
}
// instantiate
template
DryMatrix<double>::DryMatrix(
  int rows, int cols, int sym, SourceLocation const &location
);
template
DryMatrix<complex>::DryMatrix(
  int rows, int cols, int sym, SourceLocation const &location
);

template <typename F>
DryVector<F>::DryVector(
  int elements, SourceLocation const &location
):
  DryTensor<F>(
    1, std::array<int,1>{{elements}}.data(),
    std::array<int,1>{{NS}}.data(),
    location
  )
{
}
// instantiate
template
DryVector<double>::DryVector(int elements, SourceLocation const &location);
template
DryVector<complex>::DryVector(int elements, SourceLocation const &location);


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
DryScalar<double>::DryScalar(SourceLocation const &location);
template
DryScalar<complex>::DryScalar(SourceLocation const &location);

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
DryScalar<double>::DryScalar(
  double const value, SourceLocation const &location
);
template
DryScalar<complex>::DryScalar(
  complex const value, SourceLocation const &location
);

