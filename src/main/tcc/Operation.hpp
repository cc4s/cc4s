#ifndef TCC_OPERATION_DEFINED
#define TCC_OPERATION_DEFINED

#include <tcc/Costs.hpp>

namespace cc4s {
  template <typename TE>
  class Operation {
  public:
    /**
     * \brief Total executed floating point operations.
     **/
    static size_t flops;

    virtual ~Operation() {
    }

    virtual void execute() = 0;

    virtual size_t getLatestSourceVersion() = 0;

    virtual operator std::string () const = 0;

    /**
     * \brief Costs to evaluate this operation in time and memory
     **/
    Costs costs;

    std::string file;
    size_t line;

  protected:
    Operation(
      const Costs &costs_, const std::string &file_, const size_t line_
    ): costs(costs_), file(file_), line(line_) {
    }
    class ProtectedToken {
    };
  };

  template<typename TE>
  size_t cc4s::Operation<TE>::flops = 0;
}

#endif

