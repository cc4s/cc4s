#ifndef FOCK_VECTOR_DEFINED
#define FOCK_VECTOR_DEFINED

#include <math/ComplexTensor.hpp>

#include <vector>
#include <string>
#include <initializer_list>

#include <ctf.hpp>

namespace cc4s {
  template <typename F=double>
  class FockVector {
  public:
    std::vector<CTF::Tensor<F>> componentTensors;
    std::vector<std::string> componentIndices;

    FockVector(
      FockVector<F> const &a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices)
    {
    }
  
    FockVector(
      std::initializer_list<CTF::Tensor<F>> componentTensors_
      std::initializer_list<std::string> componentIndices_
    ):
      componentTensors(componentTensors_),
      componentIndices(componentIndices_)
    {
      initializeIndexTranslation();
    }

    FockVector<F> *operator =(const FockVector<F> &a) {
      // asuume shape of rhs
      componentTensors = a.componentTensors;
      componentIndices = a.componentIndices;
      initializeIndexTranslation();
      return *this;
    }

    FockVector<F> &operator += (FockVector<F> &a) {
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        // TODO: check whether given vector a is compatible w.r.t.
        // number and shape components
        char const *indices(componentIndices[i].c_str());
        componentTensors[i][indices] += a.componentTensors[i][indices];
      }
      return *this;
    }

    FockVector<F> &operator *= (F const s) {
      // TODO: implement scalar multiplication
    }

    F dot(FockVector<F> &a) {
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        // TODO: check whether given vector a is compatible w.r.t.
        // number and shape components
        char const *indices(componentIndices[i].c_str());
        CTF::Tensor<F> conjThis(componentTensors[i]);
        cc4s::conjugate(conjThis);
        // add to result
        result[""] += conjThis[indices] * a.componentTensors[i][indices];
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
        tensorIndices[i].reserve(tensorIndices[i].size()+1);
        tensorValues[i].reserve(tensorIndices[i].size()+1);
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
      while (component < static_cast<int>(indexEnds.size())) {
        if (totalIndex < indexEnds[component]) break;
        base = indexEnds[component++];
      }
      if (component >= static_cast<int>(indexEnds.size())) {
        throw new EXCEPTION("Index out bounds");
      }
      componentIndex = totalIndex - base;
    }

    /**
     * \brief Check if two fock vectors are dual conjugated of each other.
     * e.g.:
     *  Suppose we have Aabij and Bijab,
     *    \f[ A^{ab}_{ij} = c^{+}_a c^{+}_b c_j c_i \f]
     *  and also
     *    \f[ (A^{ab}_{ij})^+ = c^{+}_i c^{+}_j c_b c_a \f] = A^{ij}_{ab}
     *  therefore it should check that the length
     *  of the index a in Aabij is the same as the length of
     *  a in Bijab, the len of b in Aabij the same as b in Bijab, etc...
     *
     *  TODO: Improve speed?
     */
    void checkDualCompatibility(const FockVector<F> &a) const {
      checkCompatibilityTo(a);
      for (int i ; i < componentTensors.size() ; i++) {
        CTF::Tensor<F> *A(a.componentTensors[i]), *B(componentTensors[i]);
        std::string aIndices(a.componentTensors[i].c_str()),
          bIndices(componentTensors[i].c_str());
        int indexLens(A->order());
        for (int j ; j < indexLens ; j++) {
          std::size_t bIndexPos(
            bIndices.find(aIndices[j])
          );
          if (bIndexPos == std::string::npos) {
            throw EXCEPTION("Index of fock vectors do not match");
          }
          if (A->lens[j] != B->lens[bIndexPos]) {
            throw EXCEPTION("Number of component tensors does no match");
          }
        }
      }
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
  FockVector<F> inline operator +(FockVector<F> &a, FockVector<F> &b) {
    FockVector<F> result(a);
    result += b;
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
