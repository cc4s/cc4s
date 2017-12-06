/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all rights reserved.*/
#ifndef MP2_EOM_DEFINED
#define MP2_EOM_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

namespace cc4s {
  class UccsdAmplitudesFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UccsdAmplitudesFromCoulombIntegrals);
    UccsdAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~UccsdAmplitudesFromCoulombIntegrals();

    // Mask, one and two body parts
    PTR(CTF::Tensor<>) Mai, Mabij;

    bool maskIsGiven();
    virtual void run();
    void createMask();
    virtual std::string getAbbreviation() { return "Uccsd"; }

  protected:
    /**
     * \brief Implements the iterate method with the DRCCD iteration.
     * \param[in] i Iteration number
     */
    virtual PTR(FockVector<double>) getResiduum(
      const int iteration, const PTR(const FockVector<double>) &amplitudes
    );

    virtual PTR(FockVector<complex>) getResiduum(
      const int iteration, const PTR(const FockVector<complex>) &amplitudes
    );

  };

  /**
   * \brief Class to parse a string of comma separated range delimiters.
   * Note: It only works for unsigned integers.
   */
  class RangeParser {
    public:
      /**
       * \brief Constructor from a string, it will parse the string
       * automatically.
       */
      RangeParser(std::string rawRange_) : rawRange(rawRange_) {
        parse();
      }
      /**
       * \brief Get the range in the form of a vector of integers
       */
      std::vector<int> getRange() const {
        return parsedRange;
      }

      /**
       * \brief Parse individual range atoms.
       * For example this method will map "2" into std::vector<int>{2}
       * and "2-3" into std::vector<int>{2,3}.
       */
      std::vector<int> atomicRangeToRange(const std::string atomicRange) {
        std::vector<int> range;
        int low(0), high(0);
        int dash_pos(atomicRange.find("-"));
        if (dash_pos == -1) {
          low = std::stoi(atomicRange);
          high = low;
        } else {
          low = std::stoi(atomicRange.substr(0, dash_pos));
          high = std::stoi(atomicRange.substr(
            dash_pos+1, atomicRange.length() - dash_pos - 1
          ));
        }
        #ifdef DEBUG
        LOG(1, "Parser") << "Low " << low << std::endl;
        LOG(1, "Parser") << "high " << high << std::endl;
        #endif
        for (int i(std::min(low, high)) ;
            i <= std::max(low, high) ; i++) {
          range.push_back(i);
        }
        return range;
      }

      /**
       * \brief Main parser method, it will parse a comma separated string
       * of range delimiters.
       *
       * E.g. it will turn "2  ,3  , 4, 8-10" into {2,3,4,8,9,10}.
       * Note that it is supposed to be whitespace safe.
       */
      void parse() {
        std::string buff("");
        char comma = ',', space = ' ';
        for (unsigned int i(0) ; i < rawRange.length() ; i++) {
          if (rawRange[i] == comma) {
          } else if (rawRange[i] == space) {
            continue;
          } else {
            buff.append(rawRange, i, 1);
            if (i != rawRange.length() - 1){
              continue;
            }
          }
          // Parse the numbers here
          auto tempRange(atomicRangeToRange(buff));
          parsedRange.insert(
            parsedRange.end(), tempRange.begin(), tempRange.end()
          );
          buff.clear();
        }
      }

    int get_max() const {
      auto max_it(std::max_element(parsedRange.begin(), parsedRange.end()));
      return *max_it;
    }

    protected:
      const std::string rawRange;
      std::vector<int> parsedRange;
  };

  inline std::ostream &operator<<(
    std::ostream &stream, const RangeParser &parser
  ) {
    for (auto i : parser.getRange()) {
      stream << i << " ";
    }
    return stream << std::endl;
  }

}

#endif

