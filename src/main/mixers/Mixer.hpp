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

#ifndef MIXER_DEFINED
#define MIXER_DEFINED

#include <Data.hpp>
#include <math/Complex.hpp>
#include <math/TensorUnion.hpp>
#include <util/SharedPointer.hpp>

#include <string>

namespace cc4s {
  template <typename F, typename TE>
  class Mixer {
  public:
    Mixer(const Ptr<MapNode> &arguments);
    virtual ~Mixer();

    /**
     * \brief Returns the name of the implenting mixer.
     **/
    virtual std::string getName() = 0;

    /**
     * \brief Appends the given pair (A,R) of TensorUnions to the mixer,
     * where R is the residuum when using the amplitudes A.
     * The mixer may use the given amplitudes and residua to provide
     * an estimated amplitude with a lower expected residuum.
     * A and R are not expected to change upon return.
     **/
    virtual void append(
      const Ptr<TensorUnion<F,TE>> &A,
      const Ptr<TensorUnion<F,TE>> &R
    ) = 0;

    /**
     * \brief Returns the current best estimate of the amplitudes
     * according to previously given pairs of amplitudes and residua.
     * Requires one or more previous calls to append.
     * The returned TensorUnions must not be changed.
     **/
    virtual Ptr<const TensorUnion<F,TE>> get() = 0;

    /**
     * \brief Returns the estimated residuum of the current best estimate
     * of the amplitdues according to previously given pairs of amplitudes
     * and residua.
     * Requires one or more previous calls to append.
     * The returned TensorUnions must not be changed.
     **/
    virtual Real<> getResiduumNorm() = 0;

    Ptr<MapNode> arguments;
  };

  template <typename F, typename TE>
  class MixerFactory {
  public:
    /**
     * \brief Creates a mixer object of the mixer type specified
     * by the given name.
     * The instantiated mixer must be registered using the
     * MixerRegistrar class.
     */
    static Ptr<Mixer<F,TE>> create(
      std::string const &name, const Ptr<MapNode> &arguments
    ) {
      auto iterator(getMixerMap()->find(name));
      return iterator != getMixerMap()->end() ?
        iterator->second(arguments) : Ptr<Mixer<F,TE>>();
    }
  protected:
    typedef std::map<
      std::string,
      std::function<Ptr<Mixer<F,TE>> (const Ptr<MapNode> &arguments)>
    > MixerMap;

    static Ptr<MixerMap> getMixerMap() {
      return mixerMap ? mixerMap : (mixerMap = New<MixerMap>());
    }
    static Ptr<MixerMap> mixerMap;
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename F, typename TE, typename MixerType>
  Ptr<Mixer<F,TE>> createMixer(const Ptr<MapNode> &arguments) {
    return New<MixerType>(arguments);
  }

  /**
   * \brief Class to be statically instantiated by a mixer to register
   * it in the MixerFactory. Registered mixers can be instantiated
   * from the cc4s control language.
   */
  template <typename F, typename TE, typename MixerType>
  class MixerRegistrar: protected MixerFactory<F,TE> {
  public:
    /**
     * \brief Constructs the registrating instance. The mixer type
     * must be given as template argument, the mixer name as
     * method argument.
     */
    MixerRegistrar(std::string const &name) {
      (*MixerFactory<F,TE>::getMixerMap())[name] = &createMixer<F,TE,MixerType>;
    }
  };

  /**
   * \brief Auxiliary macro declaring the mixer registrar for
   * the mixer type of the given name. This macro is to be
   * used in the mixer declaration within the .hpp file.
   * Note that name is a symbol name not a string.
   */
  #define MIXER_REGISTRAR_DECLARATION(NAME) \
    virtual std::string getName() { return #NAME; } \
    static MixerRegistrar<F,TE,NAME<F,TE>> registrar_;

  /**
   * \brief Auxiliary macro defining the mixer registrar for
   * the mixer type of the given name. This macro is to be
   * used in the mixer definition within the .cxx file.
   * Note that name is a symbol name not a string.
   */
  #define MIXER_REGISTRAR_DEFINITION(NAME) \
    template <typename F, typename TE> \
    MixerRegistrar<F,TE,NAME<F,TE>> NAME<F,TE>::registrar_(#NAME);
}

#endif

