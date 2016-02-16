/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <math/CanonicalPolyadicDecomposition.hpp>
#include <math/Complex.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

using namespace CTF;
using namespace cc4s;

template <typename F=double>
void cc4s::composeCanonicalPolyadicDecompositionTensors(
  Tensor<F> &A, Tensor<F> &B, Tensor<F> &C,
  Tensor<F> &T
) {
  Assert(
    A.order == 2 && B.order == 2 && C.order == 2 && T.order == 3 &&
    A.lens[0]==T.lens[0] && B.lens[0]==T.lens[1] && C.lens[0]==T.lens[2] &&
    A.lens[1] == B.lens[1] && B.lens[1] == C.lens[1] && C.lens[1] == A.lens[1],
   "Incompatible tensor shapes for CPD"
  );
  // choose contraction order with minimal memory footprint
  int largestIndex(
    std::max(
      std::max(T.lens[0], T.lens[1]), T.lens[2]
    )
  );
  if (A.lens[0] == largestIndex) {
    LOG(4, "CPD") << "calculating T with largest A..." << std::endl;
    int lens[] = { A.lens[1], B.lens[0], C.lens[0] };
    int syms[] = { NS, NS, NS };
    Tensor<F> BC(3, lens, syms, *T.wrld, "BC", T.profile);
    BC["Rjk"] = B["jR"] * C["kR"];
    T["ijk"] = A["iR"] * BC["Rjk"];
  } else if (B.lens[0] == largestIndex) {
    LOG(4, "CPD") << "calculating T with largest B..." << std::endl;
    int lens[] = { A.lens[1], C.lens[0], A.lens[0] };
    int syms[] = { NS, NS, NS };
    Tensor<F> CA(3, lens, syms, *T.wrld, "CA", T.profile);
    CA["Rki"] = C["kR"] * A["iR"];
    T["ijk"] = B["jR"] * CA["Rki"];
  } else {
    LOG(4, "CPD") << "calculating T with largest C..." << std::endl;
    int lens[] = { A.lens[1], A.lens[0], B.lens[0] };
    int syms[] = { NS, NS, NS };
    Tensor<F> AB(3, lens, syms, *T.wrld, "AB", T.profile);
    AB["Rij"] = A["iR"] * B["jR"];
    T["ijk"] = C["kR"] * AB["Rij"];
  }
}

// instantiate
template
void cc4s::composeCanonicalPolyadicDecompositionTensors(
  Tensor<double> &A, Tensor<double> &B, Tensor<double> &C,
  Tensor<double> &T
);
template
void cc4s::composeCanonicalPolyadicDecompositionTensors(
  Tensor<complex> &A, Tensor<complex> &B, Tensor<complex> &C,
  Tensor<complex> &T
);


template <typename F=double>
void cc4s::contractWithCanonicalPolyadicDecompositionTensors(
  Tensor<F> &T, char const *indicesT,
  Tensor<F> &B, char const idxB,
  Tensor<F> &C, char const idxC,
  Tensor<F> &A, char const idxA
) {
  Assert(
    A.order == 2 && B.order == 2 && C.order == 2 && T.order == 3 &&
    A.lens[1] == B.lens[1] && B.lens[1] == C.lens[1] && C.lens[1] == A.lens[1],
   "Incompatible tensor shapes for CPD"
  );
  char const indicesA[] = { idxA, 'R', 0 };
  char const indicesB[] = { idxB, 'R', 0 };
  char const indicesC[] = { idxC, 'R', 0 };
  // choose contraction order with minimal memory footprint
  int largestIndex(
    std::max(
      std::max(T.lens[0], T.lens[1]), T.lens[2]
    )
  );
  if (A.lens[0] == largestIndex) {
    // A has largest index: contract B and C first
    LOG(4, "CPD") << "applying to T with largest A..." << std::endl;
    const char indicesBC[] = { idxB, idxC, 'R' , 0};
    int lens[] = { B.lens[0], C.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<F> BC(3, lens, syms, *T.wrld, "BC");
    BC[indicesBC] = B[indicesB] * C[indicesC];
    A[indicesA] = T[indicesT] * BC[indicesBC];
  } else if (B.lens[0] == largestIndex) {
    // B has largest index: contract T and B first
    LOG(4, "CPD") << "applying to T with largest B..." << std::endl;
    const char indicesTB[] = { idxA, idxC, 'R' , 0};
    int lens[] = { A.lens[0], C.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<F> TB(3, lens, syms, *T.wrld, "TB");
    TB[indicesTB] = T[indicesT] * B[indicesB];
    A[indicesA] = TB[indicesTB] * C[indicesC];
  } else {
    // C has largest index: contract T and C first
    LOG(4, "CPD") << "applying to T with largest C..." << std::endl;
    const char indicesTC[] = { idxA, idxB, 'R' , 0};
    int lens[] = { A.lens[0], B.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<F> TC(3, lens, syms, *T.wrld, "TC");
    TC[indicesTC] = T[indicesT] * C[indicesC];
    A[indicesA] = TC[indicesTC] * B[indicesB];
  }
}

// instantiate
template
void cc4s::contractWithCanonicalPolyadicDecompositionTensors(
  Tensor<double> &T, char const *indicesT,
  Tensor<double> &conjB, char const idxB,
  Tensor<double> &conjC, char const idxC,
  Tensor<double> &A, char const idxA
);
template
void cc4s::contractWithCanonicalPolyadicDecompositionTensors(
  Tensor<complex> &T, char const *indicesT,
  Tensor<complex> &conjB, char const idxB,
  Tensor<complex> &conjC, char const idxC,
  Tensor<complex> &A, char const idxA
);

