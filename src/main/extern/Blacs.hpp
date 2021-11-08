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

