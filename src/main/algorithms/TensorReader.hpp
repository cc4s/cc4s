#ifndef TENSOR_READER_DEFINED 
#define TENSOR_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  using Natural = uint64_t;

  class TensorDimensionProperty {
  public:
    // property name for reference in tensor expressions,
    // e.g. "spin" for the notation "i.spin"
    std::string name;
    // map: index -> property
    std::map<Natural,Natural> propertyOfIndex;
    // index for reverse-lookup of the map index -> property
    std::map<Natural,std::set<Natural>> indicesOfProperty;
  };

  class TensorDimension {
  public:
    // dimension name for reference in tensor expressions,
    // e.g. "State" for indicating that a dimension referes to a state p
    std::string name;
    std::map<std::string, Ptr<TensorDimensionProperty>> properties;
  };

  class TensorDimensionPropertyReference {
  public:
    Natural dimension;
    Ptr<TensorDimensionProperty> property;
  };

  class TensorNonZeroCondition {
  public:
    std::vector<TensorDimensionPropertyReference> dimensionPropertyReferences;
    // list tuples of dimension properties
    std::vector<std::vector<Natural>> tuples;
  };

  class TensorNonZeroConditions {
  public:
    // NOTE: currently, only "all" is supported, i.e. all of the listed
    // conditions must be met for a tensor entry to be non-zero.
    std::vector<Ptr<TensorNonZeroCondition>> all;
  };


  class TensorNonZeroBlockElementIterator {
  public:
    std::vector<Natural> lens;
    std::vector<std::vector<Natural>> indices;
    std::vector<Natural> indexPositions;

    TensorNonZeroBlockElementIterator(
      const std::vector<Natural> &lens_,
      std::vector<std::vector<Natural>> &&indices_
    ): lens(lens_), indices(indices_), indexPositions(indices_.size(), 0) {
    }

    TensorNonZeroBlockElementIterator &operator ++() {
      if (!atEnd()) {
        ++(indexPositions[0]);
        Natural d(0);
        while (
          d+1 < indexPositions.size()  &&
          indexPositions[d] >= indices[d].size()
        ) {
          // carry?
          indexPositions[d] = 0;
          ++d;
          ++(indexPositions[d]);
        }
      }
      return *this;
    }

    bool atEnd() const {
      if (indexPositions.size() == 0) return true;
      const Natural lastIndex(indexPositions.size()-1);
      return indexPositions[lastIndex] >= indices[lastIndex].size();
    }

    Natural getGlobalIndex() const {
      Natural I(0), D(indexPositions.size()-1);
      for (Natural d(0); d < indexPositions.size(); ++d) {
        I *= lens[D-d];
        I += indices[D-d][indexPositions[D-d]];
      }
      return I;
    }

    Natural getCount() const {
      Natural count(1);
      for (Natural d(0); d < indices.size(); ++d) {
        count *= indices[d].size();
      }
      return count;
    }

    std::string print() {
      std::stringstream s;
      s << "(";
      for (Natural d(0); d < indexPositions.size(); ++d) {
        s << " " << indices[d][indexPositions[d]] << " ";
      }
      s << ")";
      return s.str();
    }
  };

  class TensorNonZeroBlockIterator {
  public:
    // the non-zero conditions over which this iterator iterates
    Ptr<TensorNonZeroConditions> nonZeroConditions;
    // references the position within the tuple list in each non-zero condition
    std::vector<Natural> tuplePositions;
    bool ended;

    TensorNonZeroBlockIterator(
      const Ptr<TensorNonZeroConditions> &nonZeroConditions_
    ):
      nonZeroConditions(nonZeroConditions_),
      tuplePositions(nonZeroConditions_->all.size(), 0),
      ended(false)
    {
    }

    TensorNonZeroBlockIterator(
      const Ptr<TensorNonZeroConditions> &nonZeroConditions_,
      std::vector<Natural> &&tuplePositions_
    ):
      nonZeroConditions(nonZeroConditions_),
      tuplePositions(tuplePositions_),
      ended(false)
    {
    }

    TensorNonZeroBlockIterator &operator ++() {
      if (ended) return *this;
      if (tuplePositions.size() > 0) {
        ++(tuplePositions[0]);
        Natural d(0);
        while (
          d+1 < tuplePositions.size()  &&
          tuplePositions[d] >= nonZeroConditions->all[d]->tuples.size()
        ) {
          // carry?
          tuplePositions[d] = 0;
          ++d;
          ++(tuplePositions[d]);
        }
        if (tuplePositions[d] >= nonZeroConditions->all[d]->tuples.size()) {
          ended = true;
        }
      } else {
        ended = true;
      }
      return *this;
    }

    bool atEnd() const {
      return ended;
    }

    TensorNonZeroBlockElementIterator getElementIterator(
      const std::vector<Natural> &lens
    ) {
      // build index map
      std::vector<Natural> conditionsCountOnDimension(lens.size());
      std::vector<std::vector<Natural>> matchCountOnDimensionIndex(lens.size());
      for (Natural d(0); d < lens.size(); ++d) {
        matchCountOnDimensionIndex[d].resize(lens[d]);
      }
      // go through all conditions
      Natural c(0);
      for (auto tp: tuplePositions) {
        auto nonZeroCondition(nonZeroConditions->all[c]);
        for (Natural dp(0); dp < nonZeroCondition->dimensionPropertyReferences.size(); ++dp) {
          auto dimensionPropertyReference(
            nonZeroCondition->dimensionPropertyReferences[dp]
          );
          auto d(dimensionPropertyReference.dimension);
          auto property(dimensionPropertyReference.property);
          // one more condition on dimension d
          ++conditionsCountOnDimension[d];
          // go through all selected indices of this dimension and mark them
          auto selectedProperty(nonZeroCondition->tuples[tp][dp]);
          auto indices(property->indicesOfProperty[selectedProperty]);
          for (auto index: indices) {
            ++matchCountOnDimensionIndex[d][index];
          }
        }
        ++c;
      }
      // finally, build list of indices, mathcing all conditions
      std::vector<std::vector<Natural>> indices(lens.size());
      for (Natural d(0); d < lens.size(); ++d) {
        for (Natural i(0); i < lens[d]; ++i) {
          if (matchCountOnDimensionIndex[d][i] == conditionsCountOnDimension[d]) {
            indices[d].push_back(i);
          }
        }
      }
      return TensorNonZeroBlockElementIterator(lens, std::move(indices));
    }

    std::string print() const {
      std::stringstream s;
      Natural d(0);
      for (auto tp: tuplePositions) {
        auto nonZeroCondition(nonZeroConditions->all[d]);
        for (Natural i(0); i < nonZeroCondition->dimensionPropertyReferences.size(); ++i) {
          auto dimensionPropertyReference(
            nonZeroCondition->dimensionPropertyReferences[i]
          );
          s << " " <<
            dimensionPropertyReference.dimension << "." <<
            dimensionPropertyReference.property->name << "=" <<
            nonZeroCondition->tuples[tp][i];
        }
        s << ";";
        ++d;
      }
      return s.str();
    }
  };

  class TensorReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorReader)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    void readData(
      const Ptr<MapNode> &tensor,
      const std::string &scalarType,
      const bool binary
    );
    template <typename F, typename TE>
    Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const std::vector<Ptr<TensorDimension>> &dimensions,
      const Ptr<TensorNonZeroConditions> &nonZeroConditions,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readBinary(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const std::vector<Ptr<TensorDimension>> &dimensions,
      const Ptr<TensorNonZeroConditions> &nonZeroConditions,
      const SourceLocation &sourceLocation
    );
  };
}

#endif

