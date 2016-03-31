/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <util/DryTensor.hpp>

#include <math/Complex.hpp>
#include <ctf.hpp>
#include <array>

using namespace cc4s;

int64_t DryMemory::currentTotalSize = 0, DryMemory::maxTotalSize = 0;

template <typename F>
DryMatrix<F>::DryMatrix(
  int rows, int cols, int sym
):
  DryTensor<F>(
    2, std::array<int,2>{{rows, cols}}.data(),
    std::array<int,2>{{sym, NS}}.data()
  )
{
}
// instantiate
template
DryMatrix<double>::DryMatrix(int rows, int cols, int sym);
template
DryMatrix<complex>::DryMatrix(int rows, int cols, int sym);

template <typename F>
DryVector<F>::DryVector(
  int elements
):
  DryTensor<F>(
    1, std::array<int,1>{{elements}}.data(),
    std::array<int,1>{{NS}}.data()
  )
{
}
// instantiate
template
DryVector<double>::DryVector(int elements);
template
DryVector<complex>::DryVector(int elements);


template <typename F>
DryScalar<F>::DryScalar(
):
  DryTensor<F>(
    0, std::array<int,0>{{}}.data(),
    std::array<int,0>{{}}.data()
  )
{
}
// instantiate
template
DryScalar<double>::DryScalar();
template
DryScalar<complex>::DryScalar();

template <typename F>
DryScalar<F>::DryScalar(
  F const value // will be discarded
):
  DryTensor<F>(
    0, std::array<int,0>{{}}.data(),
    std::array<int,0>{{}}.data()
  )
{
}
// instantiate
template
DryScalar<double>::DryScalar(double const value);
template
DryScalar<complex>::DryScalar(complex const value);

