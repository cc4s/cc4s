#include <mixers/DiisMixer.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>

#include <math/IterativePseudoInverse.hpp>

#include <array>

using namespace CTF;
using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(DiisMixer);

template <typename F>
DiisMixer<F>::DiisMixer(
  Algorithm *algorithm
):
  Mixer<F>(algorithm), next(nullptr), nextResiduum(nullptr)
{
  int N( algorithm->getRealArgument("maxResidua", 4) );
  LOG(1,"DiisMixer") << "maxResidua=" << N << std::endl;

  amplitudes.resize(N);
  residua.resize(N);
  nextIndex = 0;
  count = 0;

  // generate initial overlap matrix
  std::array<int,2> lens{{N+1, N+1}};
  std::array<int,2> syms{{NS, NS}};
  B = NEW(CTF::Tensor<F>,
    lens.size(), lens.data(), syms.data(), *Cc4s::world, "B"
  );
  CTF::Tensor<F> one(false, *B);
  one["ij"] += 1.0;
  std::array<int,2> upperLeftBegin{{0,0}};
  std::array<int,2> upperRightEnd{{1,N+1}};
  std::array<int,2> lowerLeftEnd{{N+1,1}};
  B->slice(
    upperLeftBegin.data(), upperRightEnd.data(), 1.0,
    one,
    upperLeftBegin.data(), upperRightEnd.data(), +1.0
  );
  B->slice(
    upperLeftBegin.data(), lowerLeftEnd.data(), 1.0,
    one,
    upperLeftBegin.data(), lowerLeftEnd.data(), -1.0
  );
}

template <typename F>
DiisMixer<F>::~DiisMixer() {
}

template <typename F>
void DiisMixer<F>::append(
  const PTR(FockVector<F>) &A, const PTR(FockVector<F>) &R
) {
  // replace amplidue and residuum at nextIndex
  amplitudes[nextIndex] = A;
  residua[nextIndex] = R;

  // generate 1x1 matrix to enter new overlap elements
  const int N(residua.size());
  std::array<int,2> lens{{1, 1}};
  std::array<int,2> syms{{NS, NS}};
  CTF::Tensor<F> one(lens.size(), lens.data(), syms.data(), *Cc4s::world, "1");
  one["ij"] += 1.0;
  for (int i(0); i < N; ++i) {
    if (residua[i]) {
      std::array<int,2> colBegin{{nextIndex+1,i+1}};
      std::array<int,2> colEnd{{nextIndex+2,i+2}};
      std::array<int,2> rowBegin{{i+1,nextIndex+1}};
      std::array<int,2> rowEnd{{i+2,nextIndex+2}};
      std::array<int,2> oneBegin{{0,0}};
      std::array<int,2> oneEnd{{1,1}};
      F overlap( 2.0*std::real(residua[i]->dot(*R)) );
      B->slice(
        colBegin.data(), colEnd.data(), 0.0,
        one,
        oneBegin.data(), oneEnd.data(), overlap
      );
      if (i == nextIndex) continue;
      B->slice(
        rowBegin.data(), rowEnd.data(), 0.0,
        one,
        oneBegin.data(), oneEnd.data(), overlap
      );
    }
  }
  if (count < N) ++count;

  // now, pseudo-invert upper left corner of B and read out its first column
  std::array<int,2> upperLeftBegin{{0, 0}};
  std::array<int,2> lowerRightEnd{{count+1, count+1}};
  std::array<int,2> firstColEnd{{count+1,1}};
  std::vector<F> column(count+1);
  IterativePseudoInverse<F>(
    B->slice(upperLeftBegin.data(), lowerRightEnd.data()), 1e-14
  ).get().slice(
    upperLeftBegin.data(), firstColEnd.data()
  ).read_all(column.data());
  LOG(1, "DiisMixer") << "lambda" << "=" << column[0] << std::endl;

  next = NEW(FockVector<F>, *A);
  *next *= F(0);
  nextResiduum = NEW(FockVector<F>, *R);
  *nextResiduum *= F(0);
  for (int j(0); j < count; ++j) {
    int i( (nextIndex+N-j) % N );
    LOG(1, "DiisMixer") << "w^(-" << (j+1) << ")=" << column[i+1] << std::endl;
    *next += column[i+1] * *amplitudes[i];
    *nextResiduum += column[i+1] * *residua[i];
  }

  nextIndex = (nextIndex+1) % N;
}

template <typename F>
PTR(const FockVector<F>) DiisMixer<F>::get() {
  return next;
}

template <typename F>
PTR(const FockVector<F>) DiisMixer<F>::getResiduum() {
  return nextResiduum;
}

// instantiate
template class DiisMixer<cc4s::Float64>;
template class DiisMixer<cc4s::Complex64>;

