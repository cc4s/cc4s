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
      const std::string &name_
    ): name(name_), data(name_) {
    }
    Argument(
      const std::string &name_, const std::string &data_
    ): name(name_), data(data_) {
    }
    const std::string &getName() const { return name; }
    const std::string &getData() const { return data; }
  protected:
    std::string name, data;
  };

  class Algorithm {
  public:
    Algorithm(const std::vector<Argument> &argumentList);
    virtual ~Algorithm();
    virtual std::string getName() = 0;
    virtual void run() = 0;
    virtual void dryRun();

    bool isArgumentGiven(const std::string &argumentName);
    // retrieving input arguments
    const std::string &getTextArgument(const std::string &argumentName);
    const std::string &getTextArgument(
      const std::string &argumentName, const std::string &defaultValue
    );
    bool getBooleanArgument(const std::string &name);
    bool getBooleanArgument(
      const std::string &name, const bool defaultValue
    );
    int64_t getIntegerArgument(const std::string &argumentName);
    int64_t getIntegerArgument(
      const std::string &argumentName, const int64_t defaultValue
    );
    Real<> getRealArgument(const std::string &argumentName);
    Real<> getRealArgument(
      const std::string &argumentName, const Real<> defaultValue
    );
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    PTR(ESC(tcc::Tensor<F,TE>)) getTensorArgument(
      const std::string &argumentName
    );

    // typing, allocating and setting output arguments
    /**
     * \brief Specifies the location of an output tensor data.
     * \param[in] argumentName The argument name as specified in the cc4s file
     * \param[in] tensor The reference of the tensor data allocated by the
     * caller and later freed by the system if not needed any further.
     * \note
     * often explicit instantiation may be necessary, e.g.
     * \code{.cxx}
     * setTensorArgument<complex>(complexTensor);
     * \endcode
     */
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    void setTensorArgument(
      const std::string &argumentName, const PTR(ESC(tcc::Tensor<F,TE>)) &tensor
    );
    void setRealArgument(const std::string &argumentName, const Real<> value);
    void setIntegerArgument(
      const std::string &argumentName, const int64_t value
    );

  protected:
    // type promotions:
    Real<> getRealArgumentFromInteger(const PTR(IntegerData) &data);
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    F getRealArgumentFromTensor(
      const PTR(ESC(const TensorData<F,TE>)) &tensorData
    );
    template <typename F=Real<>, typename TE=DefaultTensorEngine>
    PTR(ESC(tcc::Tensor<F,TE>)) getTensorArgumentFromReal(
      const PTR(RealData) &realData
    );

    PTR(Data) getArgumentData(const std::string &argumentName);
    std::map<std::string, std::string> arguments;
  };

  class AlgorithmFactory {
  public:
    typedef std::map<
      std::string,
      std::function<PTR(Algorithm)(const std::vector<Argument> &)>
    > AlgorithmMap;

    /**
     * \brief Creates an algorithm object of the algorithm type specified
     * by the given name. The given arguments are passed to the algorithm
     * constructor.
     * The instantiated algorithm must be registered using the
     * AlgorithmRegistrar class.
     */
    static PTR(Algorithm) create(
      const std::string &name, const std::vector<Argument> &arguments
    ) {
      auto iterator(getAlgorithmMap()->find(name));
      return iterator != getAlgorithmMap()->end() ?
        iterator->second(arguments) : nullptr;
    }
  protected:
    static PTR(AlgorithmMap) getAlgorithmMap() {
      return algorithmMap ? algorithmMap : (algorithmMap = NEW(AlgorithmMap));
    }
    static PTR(AlgorithmMap) algorithmMap;
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename AlgorithmType>
  PTR(Algorithm) createAlgorithm(const std::vector<Argument> &arguments) {
    return NEW(AlgorithmType, arguments);
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

