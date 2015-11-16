/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

class Options {
  public:
    /**
     * \brief They are only needed when doing memory tests with random dummy
     * data
     */
    int nw;
    int niter;
    bool profile, storeV;
    int logLevel;

    const int DEFAULT_NITER = 8, DEFAULT_NW = 100;
    const bool DEFAULT_PROFILE = false;
    const bool DEFAULT_STORE_V = false;
    const int DEFAULT_LOG_LEVEL = 1;

    Options(int argumentCount, char **arguments);
};

#endif

