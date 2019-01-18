/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SHARED_POINTER_DEFINED
#define SHARED_POINTER_DEFINED

#include <memory>

// provide convenience macros for shared pointers

#define PTR(TYPE) std::shared_ptr<TYPE>
#define WEAK_PTR(TYPE) std::weak_ptr<TYPE>
#define NEW(TYPE, ...) std::make_shared<TYPE>(__VA_ARGS__)

#define DYNAMIC_PTR_CAST(TYPE, VALUE) std::dynamic_pointer_cast<TYPE>(VALUE)

#define THISABLE(TYPE) std::enable_shared_from_this<TYPE>
#define THIS this->shared_from_this()

// use this to enclose types containing a comma "," in the template argument
// list
#define ESC(...) __VA_ARGS__

#endif

