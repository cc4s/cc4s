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

#ifndef COUPLED_CLUSTER_METHOD_DEFINED
#define COUPLED_CLUSTER_METHOD_DEFINED

#include <Node.hpp>
#include <math/Complex.hpp>
#include <math/TensorUnion.hpp>
#include <util/SharedPointer.hpp>

#include <string>

namespace cc4s {
  template <typename F, typename TE>
  class CoupledClusterMethod {
  public:
    CoupledClusterMethod(const Ptr<MapNode> &arguments);
    virtual ~CoupledClusterMethod();

    /**
      * \brief Returns the name of the implementing coupled cluster method.
     **/
    virtual std::string getName() = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<TensorUnion<F,TE>> getResiduum(
      const Ptr<TensorUnion<F,TE>> &amplitudes
    ) = 0;

    Ptr<MapNode> arguments;
  };

  template <typename F, typename TE>
  class CoupledClusterMethodFactory {
  public:
    /**
     * \brief Creates a method object of the coupled cluster method type
     * specified by the given name.
     * The instantiated method must be registered using the
     * CoupledClusterMethodRegistrar class.
     */
    static Ptr<CoupledClusterMethod<F,TE>> create(
      std::string const &name, const Ptr<MapNode> &arguments
    ) {
      auto iterator(getMethodMap()->find(name));
      return iterator != getMethodMap()->end() ?
        iterator->second(arguments) : nullptr;
    }
  protected:
    typedef std::map<
      std::string,
      std::function<
        Ptr<CoupledClusterMethod<F,TE>> (const Ptr<MapNode> &arguments)
      >
    > CoupledClusterMethodMap;

    static Ptr<CoupledClusterMethodMap> getMethodMap() {
      return methodMap ?
        methodMap : (methodMap = New<CoupledClusterMethodMap>());
    }
    static Ptr<CoupledClusterMethodMap> methodMap;
  };

  /**
   * \brief template function creating an instance of the given class.
   */
  template <typename F, typename TE, typename MethodType>
  Ptr<CoupledClusterMethod<F,TE>> createMethod(const Ptr<MapNode> &arguments) {
    return New<MethodType>(arguments);
  }

  /**
   * \brief Class to be statically instantiated by a method to register
   * it in the CoupledClusterMethodFactory.
   * Registered methods can be instantiated
   * from the cc4s control language.
   */
  template <typename F, typename TE, typename MethodType>
  class CoupledClusterMethodRegistrar:
    protected CoupledClusterMethodFactory<F,TE>
  {
  public:
    /**
     * \brief Constructs the registrating instance. The method type
     * must be given as template argument, the method name as
     * method argument.
     */
    CoupledClusterMethodRegistrar(std::string const &name) {
      (*CoupledClusterMethodFactory<F,TE>::getMethodMap())[name] =
        &createMethod<F,TE,MethodType>;
    }
  };

  /**
   * \brief Auxiliary macro declaring the method registrar for
   * the method type of the given name. This macro is to be
   * used in the mixer declaration within the .hpp file.
   * Note that name is a symbol name not a string.
   */
  #define COUPLED_CLUSTER_METHOD_REGISTRAR_DECLARATION(NAME) \
    std::string getName() override { return #NAME; } \
    static CoupledClusterMethodRegistrar<F,TE,NAME<F,TE>> registrar_;

  /**
   * \brief Auxiliary macro defining the method registrar for
   * the method type of the given name. This macro is to be
   * used in the method definition within the .cxx file.
   * Note that name is a symbol name not a string.
   */
  #define COUPLED_CLUSTER_METHOD_REGISTRAR_DEFINITION(NAME) \
    template <typename F, typename TE> \
    CoupledClusterMethodRegistrar<F,TE,NAME<F,TE>> NAME<F,TE>::registrar_(#NAME);
}

#endif

