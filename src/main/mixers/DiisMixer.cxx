/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <mixers/DiisMixer.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <algorithms/Algorithm.hpp>
#include <Node.hpp>
#include <extern/Lapack.hpp>

using namespace cc4s;

MIXER_REGISTRAR_DEFINITION(DiisMixer)

template <typename F, typename TE>
DiisMixer<F,TE>::DiisMixer(
  const Ptr<MapNode> &arguments
):
  Mixer<F,TE>(arguments), next(nullptr), nextResiduum(nullptr)
{
  N = arguments->getValue<size_t>("maxResidua", 4);
  LOG() << "maxResidua=" << N << std::endl;

  amplitudes.resize(N);
  residua.resize(N);
  nextIndex = 0;
  count = 0;
  size_t M(N+1);
  // set up overlap matrix
  B.resize(M*M,0);
  for ( size_t i(1); i<M; i++){ B[i] = -1;}
  for ( size_t i(M); i<M*M; i=i+M){ B[i] = -1;}
}

template <typename F, typename TE>
DiisMixer<F,TE>::~DiisMixer() {
}


template <typename F, typename TE>
void DiisMixer<F,TE>::append(
  const Ptr<TensorUnion<F,TE>> &A,
  const Ptr<TensorUnion<F,TE>> &R
) {

  // replace amplidue and residuum at nextIndex
  amplitudes[nextIndex] = A;
  residua[nextIndex] = R;

  // write the overlap matrix for the new residua
  for (size_t i(0); i < N; ++i) {
    if (residua[i]) {
      F overlap( 2.0*residua[i]->dot(*R) );
      size_t j((i+1)*(N+1)+nextIndex+1);
      B[j] = overlap;
      j = (nextIndex+1)*(N+1)+i+1;
      B[j] = overlap;
    }
  }

  // now, pseudo-invert upper left corner of B and read out its first column
  if (count < N) ++count;
  size_t dim(count+1);
  std::vector<F> column(count+1);
  std::vector<F> matrix(dim*dim);
  for ( size_t m(0); m < dim; m++)
  for ( size_t n(0); n < dim; n++)
    matrix[n+dim*m] = B[n+(N+1)*m];

  if (!Cc4s::options->dryRun) {
    column = inverse(matrix, dim);
  }

  next = New<TensorUnion<F,TE>>(*A);
  *next *= F(0);
  nextResiduum = New<TensorUnion<F,TE>>(*R);
  *nextResiduum *= F(0);
//  OUT() << "\tDiis: ";
  for (size_t j(0); j < count; ++j) {
    size_t i( (nextIndex+N-j) % N );
    LOG() << "w^(-" << (j+1) << ")=" << column[i+1] << std::endl;
//    OUT() << "w^(-" << (j+1) << ")=" << column[i+1] << ", ";
    *next += column[i+1] * *amplitudes[i];
    *nextResiduum += column[i+1] * *residua[i];
  }
//  OUT() << std::endl;
  nextIndex = (nextIndex+1) % N;
  residuumNorm = sqrt(real(nextResiduum->dot(*nextResiduum)));
}

template <typename F, typename TE>
Ptr<TensorUnion<F,TE>> DiisMixer<F,TE>::get() {
  return next;
}

template <typename F, typename TE>
Real<> DiisMixer<F,TE>::getResiduumNorm() {
  return residuumNorm;
}

template <typename F, typename TE>
std::vector<Real<>> DiisMixer<F,TE>::inverse(
  std::vector<Real<>> matrix, size_t N
){
  std::vector<Real<>> column(N,0);
  std::vector<int> ipiv(N);
  std::vector<Real<>> work(N);
  column[0] = -1.0;

  int one(1); int info; int n(N);
  dsysv_("U", &n, &one, matrix.data(), &n, ipiv.data(), column.data(), &n, work.data(), &n, &info);
  if ( info != 0) THROW("Diagonalization failed");
  return column;
}

template <typename F, typename TE>
std::vector<Complex<>> DiisMixer<F,TE>::inverse(
  std::vector<Complex<>> matrix, size_t N
){
  std::vector<Complex<>> column(N,0);
  std::vector<int> ipiv(N);
  std::vector<Complex<>> work(N);
  int one(1); int info; int n(N);
  column[0] = -1.0;

  zsysv_("U", &n, &one, matrix.data(), &n, ipiv.data(), column.data(), &n, work.data(), &n, &info);
  if ( info != 0) THROW("Diagonalization failed");

  return column;
}





// instantiate
template class cc4s::DiisMixer<Real<64>, DefaultDryTensorEngine>;
template class cc4s::DiisMixer<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::DiisMixer<Real<64>, DefaultTensorEngine>;
template class cc4s::DiisMixer<Complex<64>, DefaultTensorEngine>;

