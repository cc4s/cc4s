/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_OPERATION_DEFINED
#define TCC_OPERATION_DEFINED

#include <tcc/Costs.hpp>

namespace cc4s {
  template <typename TE>
  class Operation {
  public:
    virtual ~Operation() {
    }

    virtual void execute() = 0;

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
}

#endif

