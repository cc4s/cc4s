#ifndef TCC_SCOPE_DEFINED
#define TCC_SCOPE_DEFINED

namespace cc4s {
  /**
   * \brief Contains frequency counts for all indices within tht
   * contraction including the indices on the left hand side.
   * These counts are used to efficiently determine which indices
   * are outer indices and which are contracted indices for each
   * individual contraction.
   **/
  class Scope {
  public:
    Scope(
      const std::string &file_ = "<anonymous>", const size_t line_ = 0
    ): triedPossibilitiesCount(0), file(file_), line(line_) {
      for (unsigned int i(0); i < INDICES_COUNT; ++i) counts[i] = 0;
    }

    void add(const std::string &indices, const int step = +1) {
      for (unsigned int d(0); d < indices.length(); ++d) {
        counts[static_cast<uint8_t>(indices[d])] += step;
      }
    }

    void add(const char *indices, const int step = +1) {
      for (unsigned int d(0); indices[d] != 0; ++d) {
        counts[static_cast<uint8_t>(indices[d])] += step;
      }
    }

    int operator [](char index) const {
      return counts[static_cast<uint8_t>(index)];
    }

    static constexpr unsigned int INDICES_COUNT = UINT8_MAX;
    unsigned int counts[INDICES_COUNT];

    // used during contraction compilation
    size_t triedPossibilitiesCount;

    // FIXME: user SourceLocation instead of (file,line) pair
    /**
     * \brief Source file of this scope.
     **/
    std::string file;
    /**
     * \brief Source file linenumber of this scope.
     **/
    size_t line;
  };
}

#endif

