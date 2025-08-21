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

#ifndef LAPACK_WRAPPER_DEFINED
#define LAPACK_WRAPPER_DEFINED

#include <Complex.hpp>

#include <vector>
#include <sstream>
#include <ctf.hpp>

namespace cc4s {

//TODO: replace double with Real<> stuff

  //  matmul(A, B, C, M, N, K);
  //  matmul(A, B, C, M, N, K, false, false, 2.0, 0.5):  2 * A * B + 0.5 * C
  void matmul(const std::vector<Real<>>& A, const std::vector<Real<>>& B,
              std::vector<Real<>>& C, int M, int N, int K,
              bool transA = false, bool transB = false,
              Real<> alpha = 1.0, Real<> beta = 0.0)
  {
    char tA = transA ? 'T' : 'N';
    char tB = transB ? 'T' : 'N';

    int lda = transA ? K : M;
    int ldb = transB ? N : K;
    int ldc = M;

    dgemm_(&tA, &tB,
           &M, &N, &K,
           &alpha,
           A.data(), &lda,
           B.data(), &ldb,
           &beta,
           C.data(), &ldc);
  }

  //compute_eigen(A, N, eigenvalues, true);
  void compute_eigen(std::vector<Real<>>& A,
                     int N,
                     std::vector<Real<>>& eigenvalues,
                     bool compute_vectors = true) {
    if (A.size() != N * N) {
        throw std::invalid_argument("Matrix A must be square (N x N).");
    }

    eigenvalues.resize(N);
    int lda = N;
    char jobz = compute_vectors ? 'V' : 'N';  // 'V' for eigenvectors, 'N' for only eigenvalues
    char uplo = 'U';  // Assume upper triangular part of A is used (symmetric matrix)

    int lwork = -1;  // Query optimal workspace size
    Real<> work_query;
    int info;

    // First call to get optimal work size
    dsyev_(&jobz, &uplo, &N, A.data(), &lda, eigenvalues.data(), &work_query, &lwork, &info);

    if (info != 0) {
        throw std::runtime_error("dsyev_ failed in workspace query.");
    }

    // Allocate the optimal workspace
    lwork = static_cast<int>(work_query);
    std::vector<Real<>> work(lwork);

    // Compute eigenvalues and (optionally) eigenvectors
    dsyev_(&jobz, &uplo, &N, A.data(), &lda, eigenvalues.data(), work.data(), &lwork, &info);

    if (info != 0) {
        throw std::runtime_error("dsyev_ failed during eigen decomposition.");
    }
  }

}

#endif

