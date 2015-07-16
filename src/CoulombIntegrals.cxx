/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "CoulombIntegrals.hpp"
using namespace CTF;

#include "cc4s.hpp"

/**
 * \brief Allocate all tensors
 */
CoulombIntegrals::CoulombIntegrals() {
  a = new Vector<>(cc4s::nv, *cc4s::world, "Va", cc4s::profile);
  i = new Vector<>(cc4s::no, *cc4s::world, "Vi", cc4s::profile);
  {
    int lens[] = {cc4s::nv, cc4s::no};
    int syms[] = {NS, NS};
    ai = new Tensor<>(2, lens, syms, *cc4s::world, "Vai", cc4s::profile);
  }
  {
    int lens[] = {cc4s::nv, cc4s::nv, cc4s::no, cc4s::no};
    int syms[] = {SY, NS, SY, NS};
    abij = new Tensor<>(4, lens, syms, *cc4s::world, "Vabij",cc4s::profile);
  }
// NOTE: only for testing
  {
    int lens[] = {cc4s::nv, cc4s::nv, cc4s::nv, cc4s::nv};
    int syms[] = {SY, NS, SY, NS};
    abcd = new Tensor<>(4, lens, syms, *cc4s::world, "Vabcd",cc4s::profile);
  }
}

/**
 * \brief Calculates all tensors
 */
void CoulombIntegrals::calculate() {
  int64_t indicesCount, *indices;
  double *values;
  Tensor<> *tensors[] = {
    a, i, ai, abij
// NOTE: only for testing
    , abcd
  };

  for (int i(0); i < (int)sizeof(tensors)/(int)sizeof(*tensors); ++i) {
    tensors[i]->read_local(&indicesCount, &indices, &values);
    for (int64_t j(0); j < indicesCount; ++j) {
      values[j] = ((indices[j]*16+i)%13077)/13077. -.5;
    }
    tensors[i]->write(indicesCount, indices, values);
    free(indices); free(values);
  }
}


/**
 * \brief Calculates a slice xycd of the full tensor abcd
 */
void CoulombIntegrals::calculate_xycd(Tensor<> &xycd, int a, int b) {
  int64_t indicesCount, *indices;
  double *values;
  // lengths in the full Vabcd tensor
  int absoluteLengths[] = {cc4s::nv, cc4s::nv, cc4s::nv, cc4s::nv};
  int positions[4];
  int64_t absoluteIndex;
  xycd.read_local(&indicesCount, &indices, &values);
  for (int j(0); j < indicesCount; ++j) {
    // get position within slice 
    cc4s::from_index(xycd, indices[j], positions);
    // compute absolute position in Vabcd
    positions[0] += a;
    positions[1] += b;
    // get index of absolute position
    absoluteIndex = cc4s::to_index(4, absoluteLengths, positions);
    // use it for calculating the respective element
    // NOTE: Vabcd is the array #4:
    values[j] = ((absoluteIndex*16+4)%13077)/13077. -.5;
  }
  xycd.write(indicesCount, indices, values);
  free(indices), free(values);
}

