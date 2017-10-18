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
  /**
   * \brief Represents the direct sum of CTF::Tensors and provides the
   * vector space operations of addition, scalar multiplication, inner product,
   * complex conjugation to get dual vectors and matrix multiplication
   * between vectors and duals, which yields a scalar.
   **/
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
    const PTR(CTF::Tensor<F>) &get(const size_t i) const {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component tensor.
     **/
    PTR(CTF::Tensor<F>) &get(const size_t i) {
      return componentTensors[i];
    }

    /**
     * \brief Retrieves the i-th component indices.
     **/
    const std::string &getIndices(const size_t i) const {
      return componentIndices[i];
    }

    /**
     * \brief Retrieves the i-th component indices as modifiable string.
     **/
    std::string &getIndices(const size_t i) {
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
     * \brief Add-to assignment operator adding each component of a
     * to the respective component of this FockVector.
     **/
    FockVector<F> &operator += (const FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size(); ++i) {
        const char *indices(componentIndices[i].c_str());
        get(i)->sum(+1.0, *a.get(i), indices, 1.0, indices);
      }
      return *this;
    }

    /**
     * \brief Subtract-from assignment operator subtracting each component of a
     * from the respective component of this FockVector.
     **/
    FockVector<F> &operator -= (const FockVector<F> &a) {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        get(i)->sum(-1.0, *a.get(i), indices, 1.0, indices);
      }
      return *this;
    }

    /**
     * \brief Multiply-by assignment operator scalar multiplying each component
     * each component of this FockVector by the given scalar.
     **/
    FockVector<F> &operator *= (const F s) {
      for (size_t i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        get(i)->sum(s, *get(i), indices, 0.0, indices);
      }
      return *this;
    }

    /**
     * \brief Creates and returns the conjugate transpose of this FockVector.
     * The first and the second half of the inidices in each component are
     * swapped for the transposition. For real types F the conjugation
     * does nothing.
     **/
    FockVector<F> conjugateTranspose() const {
      FockVector<F> result;
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t order(getIndices(i).length() / 2);
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

    /**
     * \brief Returns the matrix product of this bra-FockVector with the
     * given dual ket-FockVector ket.
     **/
    F braket(const FockVector<F> &ket) const {
      checkDualCompatibility(ket);
      CTF::Scalar<F> result;
      for (size_t i(0); i < componentTensors.size(); ++i) {
        const char *indices(getIndices(i).c_str());
        const char *ketIndices(ket.getIndices(i).c_str());
        // add to result
        result[""] += (*get(i))[indices] * (*ket.get(i))[ketIndices];
      }
      return result.get_val();
    }

    /**
     * \brief Returns the inner product of this ket-FockVector with the
     * given ket-FockVector a. The elements of this FockVector are conjugated
     * in the inner product, i.e. this->dot(a) yields the same results as
     * this->conjugateTranspose().braket(a).
     **/
    F dot(const FockVector<F> &a) const {
      checkCompatibilityTo(a);
      CTF::Scalar<F> result;
      for (size_t i(0); i < componentTensors.size(); ++i) {
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
     * \brief Get the number of component tensors of this FockVector.
     */
    size_t getComponentsCount() const {
      return componentTensors.size();
    }

    /**
     * \brief Get the total number of degrees of freedom represented by this
     * FockVector, i.e. the total number of field values contained in all
     * component tensors. The indices used by read and write are between
     * 0 and getDimension()-1.
     */
    size_t getDimension() const {
      return indexEnds.back();
    }

    /**
     * \Brief Translates the given component and component index into
     * its element into an index between 0 and getDimension()-1.
     **/
    size_t getIndex(const size_t component, const size_t componentIndex) const {
      size_t base(component > 0 ? indexEnds[component-1] : 0);
      return base + componentIndex;
    }

    /**
     * \Brief Translates the given index between 0 and getDimension()-1
     * into a component number and component index into the corresponding
     * component tensor.
     **/
    void fromIndex(
      const size_t index, size_t &component, size_t &componentIndex
    ) const {
      component = 0;
      size_t base(0);
      while (component < indexEnds.size()) {
        if (index < indexEnds[component]) break;
        base = indexEnds[component++];
      }
      if (component >= indexEnds.size()) {
        throw new EXCEPTION("Index out bounds");
      }
      componentIndex = index - base;
    }

    /**
     * \brief Reads out all locally stored values together with their
     * respective indices. The indices are between 0 and getDimension()-1.
     **/
    std::vector<std::pair<size_t,F>> readLocal() const {
      size_t elementsCount(0);
      std::vector<std::pair<size_t,F>> elements;
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t componentValuesCount;
        size_t *componentIndices;
        F *componentValues;
        get(i)->read_local(
          reinterpret_cast<int64_t *>(&componentValuesCount),
          reinterpret_cast<int64_t **>(&componentIndices),
          &componentValues
        );

        elements.resize(elementsCount+componentValuesCount);
        for (size_t k(0); k < componentValuesCount; ++k) {
          // translate index within component tensor to FockVector index
          elements[elementsCount+k].first = getIndex(i, componentIndices[k]);
          elements[elementsCount+k].second = componentValues[k];
        }
        elementsCount += componentValuesCount;
        free(componentIndices);
        free(componentValues);
      }
      return elements;
    }

    /**
     * \brief Writes the given values together with their
     * respective indices. The indices are between 0 and getDimension()-1.
     **/
    void write(const std::vector<std::pair<size_t,F>> &elements) {
      // vectors to contain indices and values for each component tensor
      std::vector<std::vector<size_t>> tensorIndices(componentTensors.size());
      std::vector<std::vector<F>> tensorValues(componentTensors.size());

      for (size_t k(0); k < elements.size(); ++k) {
        size_t component;
        size_t componentIndex;
        fromIndex(elements[k].first, component, componentIndex);
        // write to respective component tensor
        tensorIndices[component].push_back(componentIndex);
        tensorValues[component].push_back(elements[k].second);
      }

      // write data of each tensor
      for (size_t i(0); i < componentTensors.size(); ++i) {
        tensorIndices[i].reserve(tensorIndices[i].size()+1);
        tensorValues[i].reserve(tensorIndices[i].size()+1);
        get(i)->write(
          tensorIndices[i].size(),
          reinterpret_cast<int64_t *>(tensorIndices[i].data()),
          tensorValues[i].data()
        );
      }
    }

  protected:
    /**
     * \brief The end of the FockVector index range for each component.
     * This vector is used for translating component number and indices
     * into FockVector indicies.
     **/
    std::vector<size_t> indexEnds;

    /**
     * \Brief Builds the index ends vector needed for the
     * index translation methods getIndex and fromIndex.
     **/
    void buildIndexTranslation() {
      indexEnds.resize(componentTensors.size());
      size_t indexBase(0);
      for (size_t i(0); i < componentTensors.size(); ++i) {
        size_t tensorIndexSize(1);
        for (int d(0); d < get(i)->order; ++d) {
          tensorIndexSize *= get(i)->lens[d];
        }
        indexEnds[i] = indexBase += tensorIndexSize;
      }
    }

    /**
     * \brief Sets this FockVector's component tensors by copying the given
     * component tensors. Called by copy constructors and copy assignments.
     **/
    void copyComponents(const std::vector<PTR(CTF::Tensor<F>)> &components) {
      componentTensors.resize(components.size());
      for (size_t i(0); i < components.size(); ++i) {
        componentTensors[i] = NEW(CTF::Tensor<F>, *components[i]);
      }
    }

    /**
     * \brief Check if two FockVectors are dual conjugated of each other.
     * e.g.:
     *  Suppose we have Aabij and Bijab,
     *    \f[ A^{ab}_{ij} = c^{+}_a c^{+}_b c_j c_i \f]
     *  and also
     *    \f[ (A^{ab}_{ij})^+ = c^{+}_i c^{+}_j c_b c_a \f] = A^{ij}_{ab}
     *  therefore it should check that the length
     *  of the index a in Aabij is the same as the length of
     *  a in Bijab, the len of b in Aabij the same as b in Bijab, etc...
     **/
    // TODO: Improve speed?
    void checkDualCompatibility(const FockVector<F> &a) const {
      checkCompatibilityTo(a);
      for (size_t i(0); i < componentTensors.size() ; i++) {
        size_t indexLens(a.get(i)->order());
        for (size_t j(0); j < indexLens; j++) {
          size_t indexPos( get(i).find(a.getIndicies(i)[j]) );
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

  /**
   * \brief Returns the sum of two FockVectors a and b, where
   * neither a nor b are modified.
   **/
  template <typename F>
  inline FockVector<F> operator +(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result += b;
    return std::move(result);
  }
  /**
   * \brief Returns the sum of two FockVectors a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator +(
    FockVector<F> &&a, const FockVector<F> &b
  ) {
    a += b;
    return std::move(a);
  }
  /**
   * \brief Returns the sum of two FockVectors a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator +(
    FockVector<F> &a, const FockVector<F> &&b
  ) {
    b += a;
    return std::move(b);
  }

  /**
   * \brief Returns the difference between two FockVectors a and b, where
   * neither a nor b are modified.
   **/
  template <typename F>
  inline FockVector<F> operator -(
    const FockVector<F> &a, const FockVector<F> &b
  ) {
    FockVector<F> result(a);
    result -= b;
    return std::move(result);
  }
  /**
   * \brief Returns the difference between two FockVectors a and b, where
   * a is movable and will be used for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator -(
    FockVector<F> &&a, const FockVector<F> &b
  ) {
    a -= b;
    return std::move(a);
  }
  /**
   * \brief Returns the difference between two FockVectors a and b, where
   * b is movable and will be used for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator -(
    const FockVector<F> &a, FockVector<F> &&b
  ) {
    b -= a;
    // TODO: directly invoke sum to prevent extra multiplication by -1
    b *= F(-1);
    return std::move(b);
  }

  /**
   * \brief Returns the scalar multiple of the FockVector a
   * right-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F>
  inline FockVector<F> operator *(const FockVector<F> &a, const F s) {
    FockVector<F> result(a);
    result *= s;
    return std::move(result);
  }
  /**
   * \brief Returns the scalar multiple of the FockVector a
   * right-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator *(FockVector<F> &&a, const F s) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Returns the scalar multiple of the FockVector a
   * left-multiplied with the scalar s, where a is not modified.
   **/
  template <typename F>
  inline FockVector<F> operator *(const F s, const FockVector<F> &a) {
    FockVector<F> result(a);
    result *= s;
    return std::move(result);
  }
  /**
   * \brief Returns the scalar multiple of the FockVector a
   * left-multiplied with the scalar s, where a movable and will be used
   * for the result.
   **/
  template <typename F>
  inline FockVector<F> &&operator *(const F s, FockVector<F> &&a) {
    a *= s;
    return std::move(a);
  }

  /**
   * \brief Writes the FockVector a to the given stream and returns it
   * for further stream operations.
   **/
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

