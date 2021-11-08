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

#ifndef TCC_COSTS_DEFINED
#define TCC_COSTS_DEFINED

#include <limits>

namespace cc4s {
  class Costs {
  public:
    Costs(
      int64_t const elementsCount_
    ):
      maxElementsCount(elementsCount_),
      elementsCount(elementsCount_),
      multiplicationsCount(0),
      additionsCount(0)
    {
    }

    Costs(
      int64_t const maxElementsCount_,
      int64_t const elementsCount_,
      int64_t const multiplicationsCount_,
      int64_t const additionsCount_
    ):
      maxElementsCount(maxElementsCount_),
      elementsCount(elementsCount_),
      multiplicationsCount(multiplicationsCount_),
      additionsCount(additionsCount_)
    {
    }

    Costs(
      Costs const &a
    ):
      maxElementsCount(a.maxElementsCount),
      elementsCount(a.elementsCount),
      multiplicationsCount(a.multiplicationsCount),
      additionsCount(a.additionsCount)
    {
    }

    Costs &operator +=(Costs const &a) {
      maxElementsCount = std::max(maxElementsCount, a.maxElementsCount);
      elementsCount += a.elementsCount;
      multiplicationsCount += a.multiplicationsCount;
      additionsCount += a.additionsCount;
      return *this;
    }

    /**
     * \brief Maximum number of tensor elements of storage required
     * during the evaluation .
     **/
    int64_t maxElementsCount;
    /**
     * \brief Number of tensor elements of storage required by a result.
     **/
    int64_t elementsCount;
    /**
     * \brief Number of tensor element multiplication required for
     * the evaluation of this operation.
     **/
    int64_t multiplicationsCount;
    /**
     * \brief Number of tensor elements additions required for
     * the evaluation of this operation.
     **/
    int64_t additionsCount;
  };

  inline Costs operator +(Costs const &a, Costs const &b) {
    Costs result(a);
    result += b;
    return result;
  }
}

#endif

