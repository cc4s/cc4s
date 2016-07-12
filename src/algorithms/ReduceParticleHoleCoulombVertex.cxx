#include <algorithms/ReduceParticleHoleCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Exception.hpp>
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

// FIXME: properly interface with distributed versions of LAPACK
extern "C" {
  void zheev_(
    const char *job, const char *uplo, int *n, complex *a, int *lda,
    double *w, complex *work, int *workCount, double *realWork, int *info
  );
};

void ReduceParticleHoleCoulombVertex::run() {
  Tensor<complex> *EGH(getTensorArgument<complex>("EnergyMatrix"));

  // read all elements on root
  int64_t elementsCount;
  int64_t *indices;
  complex *elements;
  if (EGH->wrld->rank == 0) {
    elementsCount = EGH->lens[0] * EGH->lens[1];
    indices = new int64_t[elementsCount];
    for (int64_t index(0); index < elementsCount; ++index) {
      indices[index] = index;
    }
    elements = new complex[elementsCount];
  } else {
    elementsCount = 0;
    indices = new int64_t[0];
    elements = new complex[0];
  }
  EGH->read(elementsCount, indices, elements);

  LOG(1, "GREDUCE") << "diagonalizing " <<
    EGH->lens[0] << "x" << EGH->lens[0] << " energy matrix..." << std::endl;

  // call LAPACK diagonalization routine
  if (EGH->wrld->rank == 0) {
    int n(EGH->lens[0]), lda(n), info;
    double *eigenValues(new double[n]);
    int workCount;
    complex *work;
    double *realWork(new double[std::max(1, 3*n-2)]);  
    {
      complex optimalWorkCount;
      workCount = -1;
      zheev_(
        "V", "L", &n, elements, &lda,
        eigenValues, &optimalWorkCount, &workCount, realWork, &info
      );
      workCount = static_cast<int>(optimalWorkCount.real());
      work = new complex[workCount];
    }
    // Solve eigenproblem
    zheev_(
      "V", "L", &n, elements, &lda,
      eigenValues, work, &workCount, realWork, &info
    );
    // Check for convergence
    if (info > 0) throw new Exception("Failed to converge eigenvectors of energy matrix");
    delete[] work;
    delete[] realWork;

    double e(0.0);
    for (int i(0); i < n; ++i) {
      LOG(2, "GREDUCE") <<
        "eigenvalue[" << i << "]=" << eigenValues[i] << std::endl;
      e += eigenValues[i];
    }
    LOG(1, "GREDUCE") << "sum of eigenvalues=" << e << std::endl;

    // TODO: determine largest eigenvalues in magnitude from negative and positive
    //       end, is there a negative energy scale?
    // TODO: build truncated U
    // TODO: maybe synchronize with barrier
    delete[] eigenValues;
  }

  // TODO: write U to ctf

  delete[] elements;
  delete[] indices;

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

