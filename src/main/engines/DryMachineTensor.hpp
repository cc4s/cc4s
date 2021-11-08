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

#ifndef DRY_MACHINE_TENSOR_DEFINED
#define DRY_MACHINE_TENSOR_DEFINED

#include <engines/DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <string>
#include <mpi.h>
#include <functional>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  template <typename EmulatedTensorEngine> class DryTensorEngine;

  /**
   * \brief MachineTensor adapter for a DryTensor
   **/
  template <typename F, typename EmulatedTensorEngine>
  class DryMachineTensor {
  protected:
    class ProtectedToken {
    };

  public:
    typedef DryTensor<F> T;
    typedef DryTensorEngine<EmulatedTensorEngine> TensorEngine;
    using ETE = EmulatedTensorEngine;

    // constructors called by factory
    DryMachineTensor(
      const std::vector<size_t> &lens,
      const std::string &name,
      const ProtectedToken &
    ):
      tensor(
        lens.size(), lens.data(),
        std::vector<int>(lens.size(), 0).data()
      )
    {
      tensor.set_name(name);
    }

    // copy constructor from DryTensor
    DryMachineTensor(const T &t, const ProtectedToken &): tensor(t) {
    }

    ~DryMachineTensor() {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<DryMachineTensor<F,ETE>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[bIndices] = f(alpha * A[aIndices]) + beta * this[bIndices]
    template <typename Domain>
    void sum(
      Domain alpha,
      const Ptr<DryMachineTensor<Domain,ETE>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const Domain)> &f
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<DryMachineTensor<F,ETE>> &A,
      const std::string &aIndices,
      const Ptr<DryMachineTensor<F,ETE>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    // this[cIndices] = alpha * g(A[aIndices],B[bIndices]) + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<DryMachineTensor<F,ETE>> &A,
      const std::string &aIndices,
      const Ptr<DryMachineTensor<F,ETE>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const std::function<F(const F, const F)> &g
    ) {
      // TODO: allocate intermediates for memory estimation
    }

    void slice(
      F alpha,
      const Ptr<DryMachineTensor<F,ETE>> &A,
      const std::vector<size_t> aBegins,
      const std::vector<size_t> aEnds,
      F beta,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      // allocate tensor for A assuming index reordering
      DryTensor<F> intermediateA(A->tensor, SOURCE_LOCATION);
      // allocate tensor for result assuming index reordering
      DryTensor<F> intermediateResult(this->tensor, SOURCE_LOCATION);
    }

    // TODO: interfaces to be defined: permute, transform

    // read tensor elements to buffer
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
    }

    void readToFile(MPI_File &file, const size_t offset = 0) {
    }

    // write tensor elements from buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
    }

    void writeFromFile(MPI_File &file, const size_t offset = 0) {
    }

    std::vector<size_t> getLens() const {
      return tensor.lens;
    }

    std::string getName() const {
      return tensor.name;
    }

    // adapted DryTensor
    T tensor;

    // TODO: protect
    // create adapater from given DryTensor
    static Ptr<DryMachineTensor<F,ETE>> create(const T &t) {
      return New<DryMachineTensor<F,ETE>>(t, ProtectedToken());
    }

    // create adapter from shape and name
    static Ptr<DryMachineTensor<F,ETE>> create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return New<DryMachineTensor<F,ETE>>(lens, name, ProtectedToken());
    }
  protected:
    friend class Tensor<F,TensorEngine>;
  };
}

#endif

