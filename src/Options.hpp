/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

class Options {
  public:
    int no, nv, nG, niter;
    bool profile;

    const int DEFAULT_NO = 4, DEFAULT_NV = 6, DEFAULT_NG = 10, DEFAULT_NITER = 1;
    const bool DEFAULT_PROFILE = false;

    Options(int argumentCount, char **arguments);
};

#endif

