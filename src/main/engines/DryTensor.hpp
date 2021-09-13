#ifndef DRY_TENSOR_DEFINED
#define DRY_TENSOR_DEFINED

#include <math/Real.hpp>
#include <util/SourceLocation.hpp>
#include <util/Log.hpp>
#include <cstdint>
#include <vector>
#include <string>

namespace cc4s {
  class DryMemory {
  public:
    typedef std::pair<size_t, SourceLocation> ExtendingResource;

    static void allocate(
      size_t size, SourceLocation const &location
    ) {
      currentTotalSize += size;
      if (currentTotalSize > maxTotalSize) {
        maxTotalSize = currentTotalSize;
        extendingResources.push_back(
          ExtendingResource(maxTotalSize, location)
        );
        LOG_LOCATION(location) << "extending memory size to" << size << std::endl;
      } else {
/*
        LOG_LOCATION(location) << "memory size " << size << std::endl;
*/
      }
    }
    static void free(size_t size) {
      currentTotalSize -= size;
    }
    static size_t currentTotalSize, maxTotalSize;
    static std::vector<ExtendingResource> extendingResources;
  };

  template <typename F=Real<>>
  class DryTensor {
  public:
    /**
     * \brief Creates a dry tensor for resource consumption estimation
     * without actually allocating its data.
     * Note that currently only the memory requirement of a DryTensor
     * is estimated. The number of FLOPS should also be estimated in the
     * future.
     * Symmetryies are also ignored at the moment.
     */
    DryTensor(
      const size_t order_, const size_t *lens_, const int *syms_,
      SourceLocation const &location_ = SourceLocation()
    ): order(order_), lens(order_), syms(order_), location(location_) {
      for (size_t i(0); i < order_; ++i) {
        lens[i] = lens_[i];
        syms[i] = syms_[i];
      }
      allocate();
    }
    DryTensor(
      DryTensor const &A, SourceLocation const &location_ = SourceLocation()
    ):
      order(A.order), lens(A.lens), syms(A.syms), location(location_),
      name(A.name)
    {
      allocate();
    }
    virtual ~DryTensor() {
      free();
    }
    virtual void use() {}

    size_t getElementsCount() const {
      size_t elementsCount(1);
      for (auto len: lens) {
        elementsCount *= len;
      }
      return elementsCount;
    }

    void set_name(std::string const &name_) {
      name = name_;
    }
    std::string const &get_name() const {
      return name;
    }

    size_t order;
    std::vector<size_t> lens;
    std::vector<int> syms;
    SourceLocation location;
    std::string name;

  protected:
    void allocate() {
      size = sizeof(F) * getElementsCount();
      DryMemory::allocate(size, location);
    }
    void free() {
      DryMemory::free(size);
    }
    size_t size;
  };

  template <typename F=Real<>>
  class DryMatrix: public DryTensor<F> {
  public:
    DryMatrix(
      size_t rowCount, size_t columnCount, int sym,
      SourceLocation const &location = SourceLocation()
    );
  };

  template <typename F=Real<>>
  class DryVector: public DryTensor<F> {
  public:
    DryVector(
      size_t elementsCount, SourceLocation const &location = SourceLocation()
    );
  };

  template <typename F=Real<>>
  class DryScalar: public DryTensor<F> {
  public:
    DryScalar(SourceLocation const &location = SourceLocation());
    DryScalar(F const value, SourceLocation const &location = SourceLocation());
  };
}

#endif

