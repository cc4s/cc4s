#include <mixers/DiisMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <algorithms/Algorithm.hpp>
#include <Data.hpp>
#include <extern/Lapack.hpp>

using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(DiisMixer);

template <typename F, typename TE>
DiisMixer<F,TE>::DiisMixer(
  const Ptr<MapNode> &arguments
):
  Mixer<F,TE>(arguments), next(nullptr), nextResiduum(nullptr)
{
  N = arguments->getValue<int>("maxResidua", 4);
  LOG(1,"DiisMixer") << "maxResidua=" << N << std::endl;

  amplitudes.resize(N);
  residua.resize(N);
  nextIndex = 0;
  count = 0;
  int M(N+1);
  // set up overlap matrix
  B.resize(M*M,0.0);
  for ( int i(1); i<M; i++){ B[i] = -1;}
  for ( int i(M); i<M*M; i=i+M){ B[i] = -1;}
}

template <typename F, typename TE>
DiisMixer<F,TE>::~DiisMixer() {
}


template <typename F, typename TE>
void DiisMixer<F,TE>::append(
  const Ptr<FockVector<F,TE>> &A,
  const Ptr<FockVector<F,TE>> &R
) {

  // replace amplidue and residuum at nextIndex
  amplitudes[nextIndex] = A;
  residua[nextIndex] = R;

  // write the overlap matrix for the new residua
  for (int i(0); i < N; ++i) {
    if (residua[i]) {
      F overlap( 2.0*std::real(residua[i]->dot(*R)) );
      int j((i+1)*(N+1)+nextIndex+1);
      B[j] = overlap;
      j = (nextIndex+1)*(N+1)+i+1;
      B[j] = overlap;
    }
  }

  // now, pseudo-invert upper left corner of B and read out its first column
  if (count < N) ++count;
  int dim(count+1);
  std::vector<F> column(count+1);
  std::vector<F> matrix(dim*dim);
  for ( int m(0); m < dim; m++)
  for ( int n(0); n < dim; n++)
    matrix[n+dim*m] = B[n+(N+1)*m];

  column = inverse(matrix, dim);

  next = New<FockVector<F,TE>>(*A);
  *next *= F(0);
  for (int j(0); j < count; ++j) {
    int i( (nextIndex+N-j) % N );
    LOG(1, "DiisMixer") << "w^(-" << (j+1) << ")=" << column[i+1] << std::endl;
    *next += column[i+1] * *amplitudes[i];
  }
  nextIndex = (nextIndex+1) % N;
}

template <typename F, typename TE>
Ptr<const FockVector<F,TE>> DiisMixer<F,TE>::get() {
  return next;
}

template <typename F, typename TE>
Ptr<const FockVector<F,TE>> DiisMixer<F,TE>::getResiduum() {
  return nextResiduum;
}

template <typename F, typename TE>
std::vector<double> DiisMixer<F,TE>::inverse(
  std::vector<double> matrix, int N
){
  std::vector<double> column(N,0);
  std::vector<int> ipiv(N);
  std::vector<double> work(N);
  column[0] = -1.0;

  int one(1); int info;
  dsysv_("U", &N, &one, matrix.data(), &N, ipiv.data(), column.data(), &N, work.data(), &N, &info);
  if ( info != 0) throw "problem diagonalization\n";
  return column;
}

template <typename F, typename TE>
std::vector<complex<double>> DiisMixer<F,TE>::inverse(
  std::vector<complex<double>> matrix, int N
){
  std::vector<complex<double>> column(N,0);
  //TODO
  return column;
}





// instantiate
template class DiisMixer<Real<64>, DryTensorEngine>;
template class DiisMixer<Complex<64>, DryTensorEngine>;
template class DiisMixer<Real<64>, DefaultTensorEngine>;
template class DiisMixer<Complex<64>, DefaultTensorEngine>;

