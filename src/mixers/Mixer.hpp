/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef MIXER_DEFINED
#define MIXER_DEFINED

#include <math/Complex.hpp>
#include <string>
#include <ctf.hpp>

namespace cc4s {
  template <typename F>
  class Mixer {
  public:
    Mixer();
    virtual ~Mixer();
    virtual std::string getName() = 0;

    /**
     * \brief Appends a copy of the given tensor to the mixer.
     * The passed tensor will not be modified, albeit it might be
     * redistributed during copying.
     */
    virtual void append(CTF::Tensor<F> &A) = 0;
    virtual CTF::Tensor<F> &getNext() = 0;

    static std::string indices(CTF::Tensor<F> const &A);
  };

  template <typename F>
  class MixerFactory {
  public:
// FIXME: find out why typedef doesn't work in this case
/*
    typedef std::map<
      std::string,
      std::function<Mixer<F> *()>
    > MixerMap;
*/
    /**
     * \brief Creates a mixer object of the mixer type specified
     * by the given name.
     * The instantiated mixer must be registered using the
     * MixerRegistrar class.
     */
    static Mixer<F> *create(std::string const &name) {
      auto iterator(getMixerMap()->find(name));
      return iterator != getMixerMap()->end() ?
        iterator->second() : nullptr;
    }
  protected:
    static std::map<std::string,std::function<Mixer<F> *()>> *getMixerMap() {
      return mixerMap ? mixerMap : (mixerMap = new std::map<std::string,std::function<Mixer<F> *()>>);
    }
    static std::map<std::string,std::function<Mixer<F> *()>> *mixerMap;
/*
    static MixerMap *getMixerMap() {
      return mixerMap ? mixerMap : (mixerMap = new MixerMap);
    }
    static MixerMap *mixerMap;
*/
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename F, typename MixerType>
  Mixer<F> *createMixer() {
    return new MixerType();
  }

  /**
   * \brief Class to be statically instantiated by a mixer to register
   * it in the MixerFactory. Registered mixers can be instantiated
   * from the cc4s control language.
   */
  template <typename F, typename MixerType>
  class MixerRegistrar: protected MixerFactory<F> {
  public:
    /**
     * \brief Constructs the registrating instance. The mixer type
     * must be given as template argument, the mixer name as
     * method argument.
     */
    MixerRegistrar(std::string const &name) {
      (*MixerFactory<F>::getMixerMap())[name] = &createMixer<F, MixerType>;
    }
  };

  /**
   * \brief Auxiliary macro declaring the mixer registrar for
   * the mixer type of the given name. This macro is to be
   * used in the algorith declaration within the .hpp file.
   * Note that name is a symbol name not a string.
   */
  #define MIXER_REGISTRAR_DECLARATION(NAME) \
    virtual std::string getName() { return #NAME; } \
    static MixerRegistrar<F, NAME<F>> registrar_;

  /**
   * \brief Auxiliary macro defining the mixer registrar for
   * the mixer type of the given name. This macro is to be
   * used in the mixer definition within the .cxx file.
   * Note that name is a symbol name not a string.
   */
  #define MIXER_REGISTRAR_DEFINITION(NAME) \
    template <typename F> \
    MixerRegistrar<F, NAME<F>> NAME<F>::registrar_(#NAME);
}

#endif

