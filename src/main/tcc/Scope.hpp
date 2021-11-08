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

