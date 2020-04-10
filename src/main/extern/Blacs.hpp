/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef BLACS_DEFINED
#define BLACS_DEFINED

extern "C" {
  void Cblacs_get(int context, int request, int *value);
  int Cblacs_gridinit(int *context, char const *order, int np_row, int np_col);
  void Cblacs_gridinfo(
    int context, int *np_row, int *np_col, int *my_row, int *my_col
  );
  void Cblacs_gridexit(int ictxt);
  void Cblacs_barrier(int ictxt, char const *order);
}

#endif

