#ifndef FOCK_VECTOR_DEFINED
#define FOCK_VECTOR_DEFINED

#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <math/MathFunctions.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <ostream>

#include <ctf.hpp>

namespace cc4s {
  template <typename F=double>
  class FockVector {
  public:
    typedef F FieldType;

    std::vector<PTR(CTF::Tensor<F>)> componentTensors;
    std::vector<std::string> componentIndices;

    /**
     * \brief Default constructor for an empty Fock vector without elements.
     **/
    FockVector() {
    }

    /**
     * \brief Move constructor taking possession of the tensors owned by a.
     **/
    FockVector(
      FockVector<F> &&a
    ):
      componentTensors(a.componentTensors),
      componentIndices(a.componentIndices),
      indexEnds(a.componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Copy constructor copying the tensors owned by a.
     **/
    FockVector(
      const FockVector<F> &a
    ):
      componentTensors(a.componentTensors.size()),
      componentIndices(a.componentIndices),
      indexEnds(a.componentTensors.size())
    {
      copyComponents(a.componentTensors);
      buildIndexTranslation();
    }

    /**
     * \brief Move constructor taking possession of the tensors given.
     **/
    FockVector(
      const std::vector<PTR(CTF::Tensor<F>)> &tensors,
      const std::vector<std::string> &indices
    ):
      componentTensors(tensors),
      componentIndices(indices),
      indexEnds(componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Move constructor taking possession of the tensors given
     * by the iterators.
     **/
    template <typename TensorsIterator, typename IndicesIterator>
    FockVector(
      TensorsIterator tensorsBegin, TensorsIterator tensorsEnd,
      IndicesIterator indicesBegin, IndicesIterator indicesEnd
    ):
      componentTensors(tensorsBegin, tensorsEnd),
      componentIndices(indicesBegin, indicesEnd),
      indexEnds(componentTensors.size())
    {
      buildIndexTranslation();
    }

    /**
     * \brief Retrieves the i-th component tensor. Note that
     * the CTF::Tensor is not const since rearrangement may be
     * required also in non-modifying tensor operations.
     **/
    const PTR(CTF::Tensor<F>) &get(const int i) const {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component tensor.
     **/
    PTR(CTF::Tensor<F>) &get(const int i) {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component indices.
     **/
    const std::string &getIndices(const int i) const {
      return componentIndices[i];
    }

    /**
     * \brief Retrieves the i-th component indices as modifiable string.
     **/
    std::string &getIndices(const int i) {
      return componentIndices[i];
    }

    /**
     * \brief Move assignment operator taking possession of the tensors
     * owned by a.
     **/
    FockVector<F> &operator =(const FockVector<F> &&a) {
      componentTensors = a.componentTensors;
      componentIndices = a.componentIndices;
      buildIndexTranslation();
      return *this;
    }

    /**
     * \brief Copy assignment operator copying the tensors owned by a.
     **/
    FockVector<F> &operator =(const FockVector<F> &a) {
      componentIndices = a.componentIndices;
      copyComponents(a.componentTensors);
      buildIndexTranslation();
      return *this;
    }

    /**
     * \brief Add-to assignment operator.
     **/
    FockVector<F> &operator += (const FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        get(i)->sum(+1.0, *a.get(i), indices, 1.0, indices);
      }
      return *this;
    }

    FockVector<F> &operator -= (const FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        get(i)->sum(-1.0, *a.get(i), indices, 1.0, indices);
      }
      return *this;
    }

    FockVector<F> &operator *= (const F s) {
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        get(i)->sum(s, *get(i), indices, 0.0, indices);
      }
      return *this;
    }

    FockVector<F> conjugateTranspose() const {
      FockVector<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        int order(getIndices(i).length() / 2);
        std::vector<int> transposedLens(get(i)->lens, get(i)->lens + 2*order);
        std::rotate(
          transposedLens.begin(),
          transposedLens.begin() + order,
          transposedLens.begin() + 2*order
        );
        result.componentTensors.push_back(
          NEW(CTF::Tensor<F>,
            transposedLens.size(), transposedLens.data(),
            get(i)->sym, *get(i)->wrld,
            (std::string(get(i)->get_name()) + "*").c_str()
          )
        );
        result.componentIndices.push_back(
          getIndices(i).substr(order, 2*order) + getIndices(i).substr(0, order)
        );
        CTF::Univar_Function<F> fConj(cc4s::conj<F>);
        result.get(i)->sum(
          1.0, *get(i), getIndices(i).c_str(),
          0.0, result.getIndices(i).c_str(), fConj
        );
      }
      return std::move(result);
    }

    F braket(const FockVector<F> &a) const {
      checkDualCompatibility(a);
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        const char *aIndices(a.getIndices(i).c_str());
        // add to result
        result[""] += (*get(i))[indices] * (*a.get(i))[aIndices];
      }
      return result.get_val();
    }

    F dot(const FockVector<F> &a) const {
      checkCompatibilityTo(a);
      CTF::Scalar<F> result;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        // We need also the indices for a in general, since we might have
        // {{T['ijab']}}.dot(Q['abij']), which is still valid, although
        // the indices are not equal
        // FIXME: shouldn't that be a braket?
        const char *aIndices(a.getIndices(i).c_str());
        CTF::Bivar_Function<F> fDot(&cc4s::dot<F>);
        // add to result
        result.contract(
          1.0, *get(i), indices, *a.get(i), aIndices,
          1.0, "", fDot
        );
      }
      return result.get_val();
    }

    /**
     * \brief Get the dimension of the FockVector.
     */
    int64_t getFockDimension() const {
      return componentTensors.size();
    }

    int64_t getDimension() const {
      return indexEnds.back();
    }

    /**
     * \brief Reads out all locally stored values together with their
     * respective indices.
     **/
    std::vector<std::pair<int64_t,F>> readLocal() const {
      int64_t elementsCount(0);
      std::vector<std::pair<int64_t,F>> elements;
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        int64_t componentValuesCount;
        int64_t *componentIndices;
        F *componentValues;
        get(i)->read_local(
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
        get(i)->write(
          tensorIndices[i].size(),
          tensorIndices[i].data(),
          tensorValues[i].data()
        );
      }
    }
  protected:
    std::vector<int64_t> indexEnds;

    void copyComponents(const std::vector<PTR(CTF::Tensor<F>)> &components) {
      componentTensors.resize(components.size());
      for (size_t i(0); i < components.size(); ++i) {
        componentTensors[i] = NEW(CTF::Tensor<F>, *components[i]);
      }
    }

    /**
     * \Brief Builds the index ends vector needed for the
     * index translation methods getTotalIndex and fromTotalIndex.
     **/
    void buildIndexTranslation() {
      indexEnds.resize(componentTensors.size());
      int64_t indexBase(0);
      for (unsigned int i(0); i < componentTensors.size(); ++i) {
        int64_t tensorIndexSize(1);
        for (int d(0); d < get(i)->order; ++d) {
          tensorIndexSize *= get(i)->lens[d];
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
      for (int i(0); i < componentTensors.size() ; i++) {
        int indexLens(a.get(i)->order());
        for (int j(0); j < indexLens; j++) {
          int indexPos( get(i).find(a.getIndicies(i)[j]) );
          if (indexPos == std::string::npos) {
            throw EXCEPTION("Indices of fock vectors do not match");
          }
          if (a.get(i)->lens[j] != get(i)->lens[indexPos]) {
            throw EXCEPTION("Shapes of component tensors does not match");
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

  // sum of vectors, copy version
  template <typename F>
  inline FockVector<F> operator +(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result += b;
    return std::move(result);
  }
  // sum of vectors, left operand movable
  template <typename F>
  inline FockVector<F> &&operator +(
    FockVector<F> &&a, const FockVector<F> &b
  ) {
    a += b;
    return std::move(a);
  }
  // sum of vectors, right operand movable
  template <typename F>
  inline FockVector<F> &&operator +(
    FockVector<F> &a, const FockVector<F> &&b
  ) {
    b += a;
    return std::move(b);
  }

  // difference of vectors, copy version
  template <typename F>
  inline FockVector<F> operator -(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result -= b;
    return std::move(result);
  }
  // difference of vectors, left operand movable
  template <typename F>
  inline FockVector<F> &&operator -(
    FockVector<F> &&a, const FockVector<F> &b
  ) {
    a -= b;
    return std::move(a);
  }
  // difference of vectors, right operand movable
  template <typename F>
  inline FockVector<F> &&operator -(
    const FockVector<F> &a, FockVector<F> &&b
  ) {
    b -= a;
    // TODO: directly invoke sum to prevent extra multiplication
    b *= F(-1);
    return std::move(b);
  }

  // scalar multiplication from right, copy version
  template <typename F>
  inline FockVector<F> operator *(const FockVector<F> &a, const F &s) {
    FockVector<F> result(a);
    result *= s;
    return std::move(result);
  }
  // scalar multiplication from right, vector operand movable
  template <typename F>
  inline FockVector<F> &&operator *(FockVector<F> &&a, const F &s) {
    a *= s;
    return std::move(a);
  }

  // scalar multiplication from left, copy version
  template <typename F>
  inline FockVector<F> operator *(const F &s, const FockVector<F> &a) {
    FockVector<F> result(a);
    result *= s;
    return std::move(result);
  }
  // scalar multiplication from left, vector operand movable
  template <typename F>
  inline FockVector<F> &&operator *(const F &s, FockVector<F> &&a) {
    a *= s;
    return std::move(a);
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const FockVector<F> &a
  ) {
    stream << "( ";
    stream << a.get(0) << "[" << a.getIndices(0) << "]";
    for (size_t i(1); i < a.componentTensors.size(); ++i) {
      stream << ", " << a.get(i) << "[" << a.getIndices(i) << "]";
    }
    return stream << " )";
  }
}

#endif

