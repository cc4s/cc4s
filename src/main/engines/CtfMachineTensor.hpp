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

#ifndef CTF_MACHINE_TENSOR_DEFINED
#define CTF_MACHINE_TENSOR_DEFINED

#include <SharedPointer.hpp>

#include <ctf.hpp>
#include <string>
#include <tcc/MemoryTracker.hpp>

// TODO: specify MPI communicator when creating CtfTensorEngine
namespace cc4s {
  template <typename F,typename TE> class Tensor;
  class CtfTensorEngine;

  class MemoryTracker;
  /**
   * \brief MachineTensor adapter for a CTF::Tensor
   **/
  template <typename F>
  class CtfMachineTensor {
  protected:
    class ProtectedToken {
    };

  public:
    typedef CTF::Tensor<F> T;
    typedef CtfTensorEngine TensorEngine;

    // constructors
    CtfMachineTensor(
      const std::vector<size_t> &lens,
      const std::string &name,
      const ProtectedToken &
    ):
      tensor(
        static_cast<int>(lens.size()),
        std::vector<int64_t>(lens.begin(), lens.end()).data(),
        std::vector<int>(0, lens.size()).data(),
        CTF::get_universe(), name.c_str()
      )
    {
    }

    // copy constructor from CTF tensor
    //   CtfMachineTensor(const T &t, const ProtectedToken &): tensor(t) {
    //   }
    CtfMachineTensor(const T &t, const ProtectedToken &) = delete;

    ~CtfMachineTensor() {
    }

    // this[bIndices] = alpha * A[aIndices] + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices
    ) {
      tensor.sum(
        alpha,
        A->tensor, aIndices.c_str(),
        beta,
        bIndices.c_str()
      );
    }

    // this[bIndices] = alpha * f(A[aIndices]) + beta*this[bIndices]
    void sum(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const F)> &f
    ) {
      tensor.sum(
        alpha,
        A->tensor, aIndices.c_str(),
        beta,
        bIndices.c_str(),
        CTF::Univar_Function<F>(f)
      );
    }

    template <typename G>
    void sum(
      G alpha,
      const Ptr<CtfMachineTensor<G>> &A,
      const std::string &aIndices,
      F beta,
      const std::string &bIndices,
      const std::function<F(const G)> &f
    ) {
      CTF::Transform<G,F>(
        std::function<void(const G, F &)>(
          [f,alpha,beta](const G x, F &y) { y = f(alpha*x) + beta*y; }
        )
      ) (
        A->tensor[aIndices.c_str()], tensor[bIndices.c_str()]
      );
    }

    // this[cIndices] = alpha * A[aIndices] * B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      const Ptr<CtfMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices
    ) {
      tensor.contract(
        alpha,
        A->tensor, aIndices.c_str(),
        B->tensor, bIndices.c_str(),
        beta,
        cIndices.c_str()
      );
    }

    // we dont need this anyways
    // instead this will be the low memory variant of provided
    // this[cIndices] = alpha * A[aIndices],B[bIndices] + beta*this[cIndices]
    void contract(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::string &aIndices,
      const Ptr<CtfMachineTensor<F>> &B,
      const std::string &bIndices,
      F beta,
      const std::string &cIndices,
      const double time_memory_weight
//      const std::function<F(const F, const F)> &g
    ) {
#ifdef CTF_LOWMEM_CONTRACTION
      tensor.contract(
        alpha,
        A->tensor, aIndices.c_str(),
        B->tensor, bIndices.c_str(),
        beta,
        cIndices.c_str(),
        time_memory_weight
      );
#endif
    }

    void slice(
      F alpha,
      const Ptr<CtfMachineTensor<F>> &A,
      const std::vector<size_t> aBegins,
      const std::vector<size_t> aEnds,
      F beta,
      const std::vector<size_t> begins,
      const std::vector<size_t> ends
    ) {
      tensor.slice(
        std::vector<int>(begins.begin(), begins.end()).data(),
        std::vector<int>(ends.begin(), ends.end()).data(),
        beta,
        A->tensor,
        std::vector<int>(aBegins.begin(), aBegins.end()).data(),
        std::vector<int>(aEnds.begin(), aEnds.end()).data(),
        alpha
      );
    }

    // TODO: interfaces to be defined: permute, transform

    // read tensor elements to buffer
    void read(
      const size_t elementsCount, const size_t *indexData, F *valueData
    ) {
      tensor.read(
        elementsCount,
        reinterpret_cast<const int64_t *>(indexData),
        valueData
      );
    }

    void readToFile(MPI_File &file, const size_t offset = 0) {
      tensor.write_dense_to_file(file, offset);
    }

    // write tensor elements to buffer
    void write(
      const size_t elementsCount, const size_t *indexData, const F *valueData
    ) {
      tensor.write(
        elementsCount,
        reinterpret_cast<const int64_t *>(indexData),
        valueData
      );
    }

    void writeFromFile(MPI_File &file, const size_t offset = 0) {
      tensor.read_dense_from_file(file, offset);
    }

    std::vector<size_t> getLens() const {
      return std::vector<size_t>(tensor.lens, tensor.lens+tensor.order);
    }

    std::string getName() const {
      return std::string(tensor.get_name());
    }

    /**
     * \brief The adapted CTF tensor
     **/
    T tensor;

  protected:
    static Ptr<CtfMachineTensor<F>> create(const T &t) {
      return NEW(CtfMachineTensor<F>, t, ProtectedToken());
    }

    static Ptr<CtfMachineTensor<F>> create(
      const std::vector<size_t> &lens,
      const std::string &name
    ) {
      return NEW(CtfMachineTensor<F>, lens, name, ProtectedToken());
    }

    friend class Tensor<F,CtfTensorEngine>;
  };

#if defined(CTF_MEMORY_LOGGER)

  class CtfMemoryTracker : public MemoryTracker {
  public:
    CtfMemoryTracker() : is_valid(true) {
      CTF::Logger::instance().set_sink(
        [this](const CTF::MemEvent& ev) { if (is_valid) consume(ev); });
    }

    ~CtfMemoryTracker() override {
      is_valid = false;
      // Make CTF stop calling back into us during late deallocs
      CTF::Logger::instance().set_sink(nullptr);
    }
    // Forward size threshold to the CTF logger
    void set_size_threshold(size_t threshold) {
      CTF::Logger::instance().set_size_threshold(threshold);
    }

    // Local reset: resets only the local snapshot counters
    size_t reset() override {
      local_start_bytes = current_bytes;
      local_peak_bytes = current_bytes;
      local_allocs = 0;
      local_deallocs = 0;
      return current_bytes;
    }

    // Global reset: resets everything
    void global_reset() {
      current_bytes = 0;
      high_watermark = 0;
      global_allocs = 0;
      global_deallocs = 0;
    }

    // Accessors
    size_t current_bytes_allocated() const override { return current_bytes; }
    size_t high_watermark_bytes() const override { return high_watermark; }
    size_t number_of_allocations() const override { return global_allocs; }
    size_t number_of_deallocations() const override { return global_deallocs; }

    // Local high-watermark inside a measured region
    size_t local_high_watermark() const {
      return local_peak_bytes - local_start_bytes;
    }
    size_t local_allocations() const { return local_allocs; }
    size_t local_deallocations() const { return local_deallocs; }

    private:
    void consume(const CTF::MemEvent& ev) {
      // This is extremely verbose
      //LOG() << ev.to_string() << std::endl;
      if (ev.type == CTF::EventType::AlignAlloc || ev.type == CTF::EventType::Alloc) {
        current_bytes += ev.size;
        if (current_bytes > high_watermark)
          high_watermark = current_bytes;
        global_allocs++;
        local_allocs++;
        if (current_bytes > local_peak_bytes)
          local_peak_bytes = current_bytes;

      } else { // Dealloc
         current_bytes -= ev.size;
         global_deallocs++;
         local_deallocs++;
      }
    }

    // Global counters
    size_t current_bytes{0};
    size_t high_watermark{0};
    size_t global_allocs{0};
    size_t global_deallocs{0};

    // Local snapshot counters
    size_t local_start_bytes{0};
    size_t local_peak_bytes{0};
    size_t local_allocs{0};
    size_t local_deallocs{0};

    bool is_valid;
  };


  inline std::unique_ptr<MemoryTracker> createCtfMemoryTracker() {
    return std::unique_ptr<MemoryTracker>(new CtfMemoryTracker());
  }

#endif //CTF_MEMORY_LOGGER
}

#endif

