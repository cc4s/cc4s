/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

class Options {
  public:
    /**
     * \brief They are only needed when doing memory tests with random dummy
     * data
     */
    int no, nv, nG;
    int niter;
    bool profile, storeV;

    const int DEFAULT_NO = 4, DEFAULT_NV = 6, DEFAULT_NG = 10, DEFAULT_NITER = 1;
    const bool DEFAULT_PROFILE = false;
    const bool DEFAULT_STORE_V = true;

    Options(int argumentCount, char **arguments);
};

#endif

