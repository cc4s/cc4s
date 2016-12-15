/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_COSTS_DEFINED
#define TCC_COSTS_DEFINED

#include <limits>

namespace tcc {
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

  // TODO: use cost model for comparison
  inline bool operator <(Costs const &a, Costs const &b) {
    return (a.multiplicationsCount + a.additionsCount) <
      (b.multiplicationsCount + b.additionsCount);
  }
}

#endif

