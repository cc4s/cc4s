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
  readEnergyMatrix();
  diagonalizeEnergyMatrix();
  truncateUnitaryTransform();
  writeUnitaryTransform();
  reduceVertex();
}

void ReduceParticleHoleCoulombVertex::dryRun() {
}


void ReduceParticleHoleCoulombVertex::readEnergyMatrix() {
  EGH = getTensorArgument<complex>("EnergyMatrix");
  nG = EGH->lens[0];

  // read all elements on root
  if (EGH->wrld->rank == 0) {
    elementsCount = EGH->lens[0] * EGH->lens[1];
    indices = new int64_t[elementsCount];
    for (int64_t index(0); index < elementsCount; ++index) {
      // pack it column major
      indices[index] = index;
    }
    elements = new complex[elementsCount];
  } else {
    elementsCount = 0;
    indices = new int64_t[0];
    elements = new complex[0];
  }
  EGH->read(elementsCount, indices, elements);
}

void ReduceParticleHoleCoulombVertex::diagonalizeEnergyMatrix() {
  LOG(1, "GREDUCE") << "diagonalizing " <<
    nG << "x" << nG << " energy matrix..." << std::endl;

  if (EGH->wrld->rank == 0) {
    // call LAPACK diagonalization routine on root only
    int lda(nG), info;
    eigenValues = new double[nG];
    int workCount;
    complex *work;
    double *realWork(new double[std::max(1, 3*nG-2)]);  
    {
      complex optimalWorkCount;
      workCount = -1;
      zheev_(
        "V", "L", &nG, elements, &lda,
        eigenValues, &optimalWorkCount, &workCount, realWork, &info
      );
      workCount = static_cast<int>(optimalWorkCount.real());
      work = new complex[workCount];
    }
    // Solve eigenproblem
    zheev_(
      "V", "L", &nG, elements, &lda,
      eigenValues, work, &workCount, realWork, &info
    );
    // Check for convergence
    if (info > 0) throw new Exception("Failed to converge eigenvectors of energy matrix");
    delete[] work;
    delete[] realWork;
  } else {
    eigenValues = new double[0];
  }
}

void ReduceParticleHoleCoulombVertex::truncateUnitaryTransform() {
  if (EGH->wrld->rank == 0) {
    // allocate elements of transform. In principle, it can be as large as E
    transformElements = new complex[elementsCount];
    ng = 0;

    double energy(0.0);
    for (int i(0); i < nG; ++i) {
      LOG(2, "GREDUCE") <<
        "eigenvalue[" << i << "]=" << eigenValues[i] << std::endl;
      energy += eigenValues[i];
    }
    LOG(1, "GREDUCE") << "sum of eigenvalues=" << energy << std::endl;

    int bottom(0), top(nG-1), column;
    // approximated energy by truncation
    double e(0);
    double accuracy(getRealArgument("accuracy", DEFAULT_ACCURACY));
    while (bottom < top && std::abs(e-energy) > std::abs(energy) * accuracy) {
      if (std::abs(eigenValues[bottom]) > std::abs(eigenValues[top])) {
        // bottom value is larger in magnitude, take bottom
        column = bottom++;
        LOG(3, "GREDUCE") << "taking from bottom=" << column << std::endl;
      } else {
        // otherwise, take top
        column = top--;
        LOG(3, "GREDUCE") << "taking from top=" << column << std::endl;
      }
      // add selected eigenvalue to approximation
      e += eigenValues[column];
      // add selected eigenvector to transform
      for (int64_t i(0); i < nG; ++i) {
        transformElements[ng*nG+i] = elements[column*nG+i];
      }
      ng++;
    }
    LOG(1, "GREDUCE") <<
      "taking " << ng << " of " << nG << " for accuracy of " <<
      std::abs(e-energy)/energy << std::endl;
  } else {
    transformElements = new complex[0];
    ng = 0;
  }
  // broadcast ng to all ranks
  MPI_Bcast(&ng, 1, MPI_INT, 0, EGH->wrld->comm);

  delete[] eigenValues;
  delete[] elements;
}

void ReduceParticleHoleCoulombVertex::writeUnitaryTransform() {
  UGH = new Matrix<complex>(nG, ng, *EGH->wrld);
  int localNg(UGH->wrld->rank == 0 ? ng : 0);
  UGH->write(localNg*nG, indices, transformElements);
  delete[] transformElements;
  delete[] indices;
}

void ReduceParticleHoleCoulombVertex::reduceVertex() {
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  int lens[] = { UGH->lens[1], GammaGai->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *Gammagai = new Tensor<complex>(
    3, lens, syms, *UGH->wrld, "Gammagai"
  );
  allocatedTensorArgument<complex>(
    "ReducedParticleHoleCoulombVertex", Gammagai
  );
  (*Gammagai)["gai"] = (*GammaGai)["Gai"] * (*UGH)["Gg"];
}

