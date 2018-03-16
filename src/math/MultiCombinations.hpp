#ifndef MULTI_COMBINATIONS_DEFINED
#define MULTI_COMBINATIONS_DEFINED

#include <vector>

namespace cc4s {
  /**
   * \brief Enumerates all multi combinations of k balls in n slots.
   * The balls are indistinguishable but the slots are not and
   * each slot can contain any number of balls.
   * There are binomial(n+k-1, k) multi combinations.
   * The enumeration is implemented by an iterator returning
   * the slot numbers for each ball for each possibility.
   **/
  class MultiCombinations {
  public:
    MultiCombinations(
      const unsigned int n_, const unsigned int k_
    ): n(n_), k(k_) {
    }

    class Iterator:
      public std::iterator<std::input_iterator_tag, std::vector<unsigned int>>
    {
    public:
      Iterator(const Iterator &rhs): n(rhs.n), slots(rhs.slots) {
      }
      Iterator(Iterator &&rhs): n(rhs.n), slots(std::move(rhs.slots)) {
      }
      Iterator &operator =(const Iterator &rhs) {
        n = rhs.n;
        slots = rhs.slots;
        return *this;
      }
      Iterator &operator =(Iterator &&rhs) {
        n = rhs.n;
        slots = std::move(rhs.slots);
        return *this;
      }
      bool operator ==(const Iterator &rhs) const {
        return slots == rhs.slots;
      }
      bool operator !=(const Iterator &lhs) const {
        return !(*this == lhs);
      }
      /**
       * \brief Pre-increment advancing the iterator to the next multi
       * combination and returning the iterator in the advanced state.
       **/
      Iterator &operator ++() {
        // start moving the first ball to next slot
        advanceBall(0);
        return *this;
      }
      /**
       * \brief Post-increment advancing the iterator to the next multi
       * combination and returning an iterator before advancement.
       **/
      Iterator operator ++(int) {
        Iterator before(*this);
        ++*this;
        return before;
      }
      /**
       * \brief Returns the slot numbers for each of the k balls
       **/
      const std::vector<unsigned int> &operator *() const {
        return slots;
      }
      /**
       * \brief Returns the slot numbers for each of the k balls
       **/
      const std::vector<unsigned int> &operator ->() const {
        return slots;
      }
      friend void swap(Iterator &lhs, Iterator &rhs);
    protected:
      Iterator(
        unsigned int n_, unsigned int k, bool begin = true
      ): n(n_), slots(k) {
        // start with putting all k balls in the first slot
        // end with putting all k balls in the n+1st slot
        for (unsigned int i(0); i < k; ++i) {
          slots[i] = begin ? 0 : n_;
        }
      }
      /**
       * \brief Advance the ith ball and possibly all subsequent ones.
       / Return which slot the ith ball ends up.
       **/
      unsigned int advanceBall(unsigned int i) {
        if (++slots[i] >= n) {
          // ball reaches end
          if (i+1 < slots.size()) {
            // not the last: advance next ball and put this ball in same slot
            slots[i] = advanceBall(i+1);
          } else {
            // last ball: iterator ends, put them all in (unexisting) slot n
            slots[i] = n;
          }
        }
        return slots[i];
      }
      unsigned int n;
      std::vector<unsigned int> slots;
      friend class MultiCombinations;
    };

  public:
    /**
     * \brief Returns the iterator to the first multi combination.
     **/
    Iterator begin() const {
      return Iterator(n,k,true);
    }
    /**
     * \brief Returns the end iterator of multi combination.
     * Although it can be dereferenced without failure the contained slot
     * numbers exceed the range of valid slots.
     **/
    Iterator end() const {
      return Iterator(n,k,false);
    }

    unsigned int n, k;
  };
}

#endif

