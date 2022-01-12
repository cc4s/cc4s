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

#ifndef LAPACK_INVERSE_DEFINED
#define LAPACK_INVERSE_DEFINED

#include <extern/Lapack.hpp>
#include <Complex.hpp>
#include <LapackMatrix.hpp>
#include <Exception.hpp>
#include <Log.hpp>

#include <vector>

namespace cc4s {
  // base template
  template <typename F=real>
  class LapackInverse;


  // specialization for complex
  template <>
  class LapackInverse<Complex64> {
  public:
    LapackInverse(
      const LapackMatrix<Complex64> &A_
    ): invA(A_) {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("Inverse requries a square matrix");
      }
      int rows(A_.getRows());
      std::vector<Complex64> work(rows*rows);
      int workSize(work.size());
      std::vector<int> rowPermutation(rows);
      int info;

      zgetrf_(
        &rows, &rows,
        invA.getValues(), &rows,
        rowPermutation.data(),
        &info
      );
      zgetri_(
        &rows,
        invA.getValues(), &rows,
        rowPermutation.data(),
        work.data(), &workSize,
        &info
      );
      if (info < 0) {
        std::stringstream stream;
        stream << "Argument " << -info << " of ZGETRI is illegal";
        throw EXCEPTION(stream.str());
      }
      if (info > 0) {
        throw EXCEPTION("Singular matrix cannot be inverted");
      }
    }

    ~LapackInverse() {
    }

    const LapackMatrix<Complex64> &get() const {
      return invA;
    }
  protected:
    LapackMatrix<Complex64> invA;
  };
}

#endif

