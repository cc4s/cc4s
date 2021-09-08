#ifndef SHARED_POINTER_DEFINED
#define SHARED_POINTER_DEFINED

#include <memory>

// provide convenience macros for shared pointers

#define PTR(TYPE) std::shared_ptr<TYPE>
#define WEAK_PTR(TYPE) std::weak_ptr<TYPE>
#define NEW(TYPE, ...) std::make_shared<TYPE>(__VA_ARGS__)

#define DYNAMIC_PTR_CAST(TYPE, VALUE) std::dynamic_pointer_cast<TYPE>(VALUE)

#define THISABLE(TYPE) std::enable_shared_from_this<TYPE>
#define THIS(CLASS) std::dynamic_pointer_cast<CLASS>(this->shared_from_this())

// use this to enclose types containing a comma "," in the template argument
// list
#define ESC(...) __VA_ARGS__

namespace cc4s {
  template <typename T>
  using Ptr = std::shared_ptr<T>;

  template <typename T>
  using WeakPtr = std::weak_ptr<T>;

  template <typename T, typename... Args>
  inline Ptr<T> New(Args&&... args) {
    return std::make_shared<T>(args...);
  }

  template <typename B, typename A>
  inline Ptr<B> dynamicPtrCast(const Ptr<A> &p) {
    return std::dynamic_pointer_cast<B>(p);
  }

  template <typename B>
  class Thisable: public std::enable_shared_from_this<B> {
  public:
    template <typename S>
    inline Ptr<S> toPtr() {
      return std::dynamic_pointer_cast<S>(this->shared_from_this());
    }
  };

}

// support dummy input to Ptr<> from streams
template <typename T>
inline std::istream &operator >>(
  std::istream &stream, const cc4s::Ptr<T> &P
) {
  // do nothing
  return stream;
}


#endif

