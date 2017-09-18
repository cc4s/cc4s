#ifndef FOCK_VECTOR_DEFINED
#define FOCK_VECTOR_DEFINED

#include <math/ComplexTensor.hpp>
#include <math/Vector.hpp>

#include <vector>
#include <string>
#include <initializer_list>

#include <ctf.hpp>

namespace cc4s {
  template <typename F=double>
  class FockVector {
  public:
    typedef F FieldType;

    std::vector<CTF::Tensor<F>> componentTensors;
    std::vector<std::string> componentIndices;

    FockVector() {
    }

    FockVector(
      const FockVector<F> &a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices),
      indexEnds(a.componentTensors.size())
    {
      initializeIndexTranslation();
    }

    FockVector(
      std::initializer_list<CTF::Tensor<F>> componentTensors_,
      std::initializer_list<std::string> componentIndices_
    ):
      componentTensors(componentTensors_),
      componentIndices(componentIndices_),
      indexEnds(componentTensors_.size())
    {
      initializeIndexTranslation();
    }

    FockVector<F> *operator =(const FockVector<F> &a) {
      // asuume shape of rhs
      componentTensors = a.componentTensors;
      componentIndices = a.componentIndices;
      return *this;
    }

    FockVector<F> &operator += (FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator -= (FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += F(-1) * a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator *= (const F s) {
      CTF::Scalar<F> scalar(s);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        componentTensors[i][indices] *= scalar[""];
      }
      return *this;
    }

    F dot(FockVector<F> &a) {
      checkCompatibilityTo(a);
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        // We need also the indices for a in general, since we might have
        // {{T['ijab']}}.dot(Q['abij']), which is still valid, although
        // the indices are not equal
        const char *aIndices(a.componentIndices[i].c_str());
        CTF::Bivar_Function<F> fDot(&cc4s::dot<F>);
        // add to result
        result.contract(
          1.0, componentTensors[i], indices, a.componentTensors[i], aIndices,
          1.0, "", fDot
        );
      }
      return result.get_val();
    }

    int64_t getDimension() const {
      return indexEnds.back();
    }

    /**
     * \brief Reads out all locally stored values together with their
     * respective indices.
     **/
    std::vector<std::pair<int64_t,F>> readLocal() {
      int64_t elementsCount(0);
      std::vector<std::pair<int64_t,F>> elements;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        int64_t componentValuesCount;
        int64_t *componentIndices;
        F *componentValues;
        componentTensors[i].read_local(
          &componentValuesCount, &componentIndices, &componentValues
        );

        elements.resize(elementsCount+componentValuesCount);
        for (int64_t k(0); k < componentValuesCount; ++k) {
          // translate index within component tensor to Fock vector index
          elements[elementsCount+k].first = getTotalIndex(
            i, componentIndices[k]
          );
          elements[elementsCount+k].second = componentValues[k];
        }
        elementsCount += componentValuesCount;
        free(componentIndices);
        free(componentValues);
      }
      return elements;
    }

    void write(const std::vector<std::pair<int64_t,F>> &elements) {
      // vectors to contain indices and values for each component tensor
      std::vector<std::vector<int64_t>> tensorIndices(componentTensors.size());
      std::vector<std::vector<F>> tensorValues(componentTensors.size());
      std::vector<int64_t> indexEnds(componentTensors.size());

      for (uint64_t k(0); k < elements.size(); ++k) {
        int component;
        int64_t componentIndex;
        fromTotalIndex(elements[k].first, component, componentIndex);
        // write to respective component tensor
        tensorIndices[component].push_back(componentIndex);
        tensorValues[component].push_back(elements[k].second);
      }

      // write data of each tensor
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        componentTensors[i].write(
          tensorIndices[i].size(),
          tensorIndices[i].data(),
          tensorValues[i].data()
        );
      }
    }
  protected:
    std::vector<int64_t> indexEnds;

    /**
     * \Brief Initializes the index ends vector needed for the
     * index translation methods getTotalIndex and fromTotalIndex.
     **/
    void initializeIndexTranslation() {
      int64_t indexBase(0);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        int64_t tensorIndexSize(1);
        for (int d(0); d < componentTensors[i].order; ++d) {
          tensorIndexSize *= componentTensors[i].lens[d];
        }
        indexEnds[i] = indexBase += tensorIndexSize;
      }
    }

    /**
     * \Brief Translates the given component and the index into
     * its element into a total index between 0 and getDimension()-1.
     **/
    int64_t getTotalIndex(
      const int component, const int64_t componentIndex
    ) const {
      int64_t base(component > 0 ? indexEnds[component-1] : 0);
      return base + componentIndex;
    }

    /**
     * \Brief Translates the given total index between 0 and getDimension()-1
     * into a component number and an index into the corresponding component
     * tensor.
     **/
    void fromTotalIndex(
      const int64_t totalIndex, int &component, int64_t &componentIndex
    ) const {
      component = 0;
      int64_t base(0);
      while (totalIndex >= indexEnds[component]) {
        base = indexEnds[component++];
      }
      componentIndex = totalIndex - base;
    }

    void checkCompatibilityTo(const FockVector<F> &a) const {
      if (
        componentTensors.size() != a.componentTensors.size() ||
        componentIndices.size() != a.componentIndices.size()
      ) {
        throw EXCEPTION("Number of component tensors does no match");
      }
      // TODO: check shapes.
    }
  };

  template <typename F>
  FockVector<F> inline operator +(
    FockVector<F> &a, FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result += b;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator -(
    FockVector<F> &a, FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result -= b;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator *(FockVector<F> &a, const F &s) {
    FockVector<F> result(a);
    result *= s;
    return result;
  }

  template <typename F>
  FockVector<F> inline operator *(const F &s, FockVector<F> &a) {
    FockVector<F> result(a);
    result *= s;
    return result;
  }
}

#endif

/*
  example usage:

  CTF::Scalar<> nulliesA;
  CTF::Tensor<> singlesA(2, vo, ...);
  CTF::Tensor<> doublesA(4, vvoo, ...);
  
  FockVector<double> fockA({{nulliesA, singlesA, doublesA}}, {{"", "ai", "abij"}});
  FockVector<double> fockB({{nulliesB, singlesB, doublesB}}, {{"", "ai", "abij"}});

  fockA += fockB;
  fockA.dot(fockB);
*/
