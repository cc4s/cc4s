/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef EXCEPTION_DEFINED
#define EXCEPTION_DEFINED

#include <string>

class Exception {
  public:
    Exception(char const *message_): message(message_) {
    }
    std::string getMessage() {
      return message;
    }
  private:
    std::string message;
};

#endif

