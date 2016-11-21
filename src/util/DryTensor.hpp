/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_DEFINED
#define DRY_TENSOR_DEFINED

#include <util/SourceLocation.hpp>
#include <util/Log.hpp>
#include <cstdint>
#include <vector>
#include <string>


namespace cc4s {
  class DryMemory {
  public:
    typedef std::pair<int64_t, SourceLocation> ExtendingResource;

    static void allocate(
      int64_t size, SourceLocation const &location
    ) {
      currentTotalSize += size;
      if (currentTotalSize > maxTotalSize) {
        maxTotalSize = currentTotalSize;
        extendingResources.push_back(
          ExtendingResource(maxTotalSize, location)
        );
        LOG(2, "DryMemory") << "extending size=" << size << " at "
          << location << std::endl;
      } else {
        LOG(2, "DryMemory") << "non-extending size=" << size << " at "
          << location << std::endl;
      }
    }
    static void free(int64_t size) {
      currentTotalSize -= size;
    }
    static int64_t currentTotalSize, maxTotalSize;
    static std::vector<ExtendingResource> extendingResources;
  };

  template <typename F=double>
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
      int order_, int const *lens_, int const *syms_,
      SourceLocation const &location_ = SourceLocation()
    ): order(order_), location(location_) {
      for (int i(0); i < order_; ++i) {
        lens.push_back(lens_[i]);
        syms.push_back(syms_[i]);
      }
      allocate();
    }
    DryTensor(
      DryTensor const &A, SourceLocation const &location_ = SourceLocation()
    ):
      order(A.order), lens(A.lens), syms(A.syms), location(location_)
    {
      allocate();
    }
    virtual ~DryTensor() {
      free();
    }
    virtual void use() {}

    int order;
    std::vector<int> lens, syms;
    SourceLocation location;

  protected:
    void allocate() {
      size = sizeof(F);
      for (int i(0); i < order; ++i) {
        size *= lens[i];
      }
      DryMemory::allocate(size, location);
    }
    void free() {
      DryMemory::free(size);
    }
    int64_t size;
  };

  template <typename F=double>
  class DryMatrix: public DryTensor<F> {
  public:
    DryMatrix(
      int rowCount, int columnCount, int sym,
      SourceLocation const &location = SourceLocation()
    );
  };

  template <typename F=double>
  class DryVector: public DryTensor<F> {
  public:
    DryVector(
      int elementsCount, SourceLocation const &location = SourceLocation()
    );
  };

  template <typename F=double>
  class DryScalar: public DryTensor<F> {
  public:
    DryScalar(SourceLocation const &location = SourceLocation());
    DryScalar(F const value, SourceLocation const &location = SourceLocation());
  };
}

#endif

