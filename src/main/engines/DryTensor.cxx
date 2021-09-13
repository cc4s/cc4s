#include <engines/DryTensor.hpp>

#include <math/Complex.hpp>
#include <array>

using namespace cc4s;

size_t DryMemory::currentTotalSize = 0, DryMemory::maxTotalSize = 0;
std::vector<DryMemory::ExtendingResource> DryMemory::extendingResources;


template <typename F>
DryMatrix<F>::DryMatrix(
  size_t rows, size_t cols, int sym, SourceLocation const &location
):
  DryTensor<F>(
    2, std::array<size_t,2>{{rows, cols}}.data(),
    std::array<int,2>{{sym, 0}}.data(),
    location
  )
{
}
// instantiate
template
DryMatrix<Real<64>>::DryMatrix(
  size_t rows, size_t cols, int sym, SourceLocation const &location
);
template
DryMatrix<Complex<64>>::DryMatrix(
  size_t rows, size_t cols, int sym, SourceLocation const &location
);

template <typename F>
DryVector<F>::DryVector(
  size_t elements, SourceLocation const &location
):
  DryTensor<F>(
    1, std::array<size_t,1>{{elements}}.data(),
    std::array<int,1>{{0}}.data(),
    location
  )
{
}
// instantiate
template
DryVector<Real<64>>::DryVector(size_t elements, SourceLocation const &location);
template
DryVector<Complex<64>>::DryVector(size_t elements, SourceLocation const &location);


template <typename F>
DryScalar<F>::DryScalar(
  SourceLocation const &location
):
  DryTensor<F>(
    0, std::array<size_t,0>{{}}.data(),
    std::array<int,0>{{}}.data(),
    location
  )
{
}
// instantiate
template
DryScalar<Real<64>>::DryScalar(SourceLocation const &location);
template
DryScalar<Complex<64>>::DryScalar(SourceLocation const &location);

template <typename F>
DryScalar<F>::DryScalar(
  F const value, // will be discarded
  SourceLocation const &location
):
  DryTensor<F>(
    0, std::array<size_t,0>{{}}.data(),
    std::array<int,0>{{}}.data(),
    location
  )
{
}
// instantiate
template
DryScalar<Real<64>>::DryScalar(
  const Real<64> value, SourceLocation const &location
);
template
DryScalar<Complex<64>>::DryScalar(
  const Complex<64> value, SourceLocation const &location
);
template
DryScalar<Real<128>>::DryScalar(
  const Real<128> value, SourceLocation const &location
);
template
DryScalar<Complex<128>>::DryScalar(
  const Complex<128> value, SourceLocation const &location
);
