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
    Algorithm(std::vector<Argument> const &argumentList);
    virtual ~Algorithm();
    virtual std::string getName() = 0;
    virtual void run() = 0;

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

  class AlgorithmFactory {
  public:
    typedef std::map<
      std::string,
      std::function<Algorithm *(std::vector<Argument> const &)>
    > AlgorithmMap;

    /**
     * \brief Creates an algorithm object of the algorithm type specified
     * by the given name. The given arguments are passed to the algorithm
     * constructor.
     * The instantiated algorithm must be registered using the
     * AlgorithmRegistrar class.
     */
    static Algorithm *create(
      std::string const &name, std::vector<Argument> const &arguments
    ) {
      auto iterator(getAlgorithmMap()->find(name));
      return iterator != getAlgorithmMap()->end() ?
        iterator->second(arguments) : nullptr;
    }
  protected:
    static AlgorithmMap *getAlgorithmMap() {
      return algorithmMap ? algorithmMap : (algorithmMap = new AlgorithmMap);
    }
    static AlgorithmMap *algorithmMap;
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename AlgorithmType>
  Algorithm *createAlgorithm(std::vector<Argument> const &arguments) {
    return new AlgorithmType(arguments);
  }

  /**
   * \brief Class to be statically instantiated by an algorithm to register
   * it in the AlgorithmFactory. Registered algorithms can be instantiated
   * from the cc4s control language.
   */
  template <typename AlgorithmType>
  class AlgorithmRegistrar: protected AlgorithmFactory {
  public:
    /**
     * \brief Constructs the registrating instance. The algorithm type
     * must be given as template argument, the algorithm name as
     * method argument.
     */
    AlgorithmRegistrar(std::string const &name) {
      (*getAlgorithmMap())[name] = &createAlgorithm<AlgorithmType>;
    }
  };

  /**
   * \brief Auxiliary macro declaring the algorithm registrar for
   * the algorithm type of the given name. This macro is to be
   * used in the algorith declaration within the .hpp file.
   * Note that name is a symbol name not a string.
   */
  #define ALGORITHM_REGISTRAR_DECLARATION(NAME) \
    virtual std::string getName() { return #NAME; } \
    static AlgorithmRegistrar<NAME> registrar_;
  /**
   * \brief Auxiliary macro defining the algorithm registrar for
   * the algorithm type of the given name. This macro is to be
   * used in the algorithm definition within the .cxx file.
   * Note that name is a symbol name not a string.
   */
  #define ALGORITHM_REGISTRAR_DEFINITION(NAME) \
    AlgorithmRegistrar<NAME> NAME::registrar_(#NAME);
}

#endif

