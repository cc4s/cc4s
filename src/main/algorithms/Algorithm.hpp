/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ALGORITHM_DEFINED
#define ALGORITHM_DEFINED

#include <Data.hpp>

// tensor engine selection
#include <engines/DryTensorEngine.hpp>
#include <engines/CtfTensorEngine.hpp>
namespace cc4s {
  typedef cc4s::CtfTensorEngine DefaultTensorEngine;
}

namespace cc4s {
  class Algorithm {
  public:
    Algorithm();
    virtual ~Algorithm();
    virtual std::string getName() = 0;
    virtual Ptr<MapNode> run(const Ptr<MapNode> &arguments) = 0;
    virtual Ptr<MapNode> dryRun(const Ptr<MapNode> &arguments);
  };

  class AlgorithmFactory {
  public:
    typedef std::map<std::string,std::function<Ptr<Algorithm>()>> AlgorithmMap;

    /**
     * \brief Creates an algorithm object of the algorithm type specified
     * by the given name. The given arguments are passed to the algorithm
     * constructor.
     * The instantiated algorithm must be registered using the
     * AlgorithmRegistrar class.
     */
    static Ptr<Algorithm> create(const std::string &name) {
      auto iterator(getAlgorithmMap()->find(name));
      return iterator != getAlgorithmMap()->end() ?
        iterator->second() : nullptr;
    }
  protected:
    static Ptr<AlgorithmMap> getAlgorithmMap() {
      return algorithmMap ? algorithmMap : (algorithmMap = New<AlgorithmMap>());
    }
    static Ptr<AlgorithmMap> algorithmMap;
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename AlgorithmType>
  Ptr<Algorithm> createAlgorithm() {
    return New<AlgorithmType>();
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
    AlgorithmRegistrar(const std::string &name) {
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

