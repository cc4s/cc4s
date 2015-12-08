/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <Data.hpp>
#include <string>
#include <ctf.hpp>

namespace cc4s {
  class Argument {
    public:
      Argument(
        std::string const &name_
      ): name(name_), data(name_) {
      }
      Argument(
        std::string const &name_, std::string const &data_
      ): name(name_), data(data_) {
      }
      std::string const &getName() const { return name; }
      std::string const &getData() const { return data; }
    protected:
      std::string name, data;
  };

  class Algorithm {
    public:
      Algorithm(std::vector<Argument const *> const &argumentList);
      virtual ~Algorithm();

      // retrieving input arguments
      std::string getTextArgument(std::string const &argumentName);
      int64_t getIntegerArgument(std::string const &argumentName);
      double getRealArgument(std::string const &argumentName);
      template <typename F=double>
      CTF::Tensor<F> *getTensorArgument(std::string const &argumentName);

      // typing, allocating and setting output arguments
      template <typename F=double>
      void allocatedTensorArgument(
        std::string const &argumentName, CTF::Tensor<F> *tensor
      );
      void setRealArgument(std::string const &argumentName, double const value);

    protected:
      Data *getArgumentData(std::string const &argumentName);
      std::map<std::string, std::string> arguments;
  };
}

  #endif

