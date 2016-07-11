#include <algorithms/ReduceParticleHoleCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ReduceParticleHoleCoulombVertex);

ReduceParticleHoleCoulombVertex::ReduceParticleHoleCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ReduceParticleHoleCoulombVertex::~ReduceParticleHoleCoulombVertex() {
}

void ReduceParticleHoleCoulombVertex::run() {
  Tensor<complex> *EGH(getTensorArgument<complex>("EnergyMatrix"));

  // TODO: read_all (on all)
  // TODO: on root: call
  /*
        int n = N, lda = LDA, info, lwork;
        complex wkopt;
        complex* work;
        double w[N], rwork[3*N-2];
        complex a[LDA*N] = {
           { 9.14,  0.00}, {-4.37,  9.22}, {-1.98,  1.72}, {-8.96,  9.50},
           { 0.00,  0.00}, {-3.35,  0.00}, { 2.25,  9.51}, { 2.57, -2.40},
           { 0.00,  0.00}, { 0.00,  0.00}, {-4.82,  0.00}, {-3.24, -2.04},
           { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 8.44,  0.00}
        };
        lwork = -1;
        zheev( "Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
        lwork = (int)wkopt.re;
        work = (complex*)malloc( lwork*sizeof(complex) );
        // Solve eigenproblem
        zheev( "Vectors", "Lower", &n, a, &lda, w, work, &lwork, rwork, &info );
        // Check for convergence
        if( info > 0 ) {
                printf( "The algorithm failed to compute eigenvalues.\n" );
                exit( 1 );
        }
        print_rmatrix( "Eigenvalues", 1, n, w, 1 );
        print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
        free( (void*)work );
*/
  // TODO: determine largest eigenvalues in magnitude from negative and positive
  //       end, is there a negative energy scale?
  // TODO: build truncated U
  // TODO: write U to ctf
  // TODO: calculate reduced Gamma

  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  double accuracy(getRealArgument("accuracy", DEFAULT_ACCURACY));
/*
  allocatedTensorArgument<complex>(
    "ReducedParticleHoleCoulombVertex", reducedParticleHoleCoulombVertex
  );
*/
}

void ReduceParticleHoleCoulombVertex::dryRun() {
}

