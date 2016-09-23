#include <algorithms/ReduceEnergyMatrix.hpp>
#include <util/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ReduceEnergyMatrix);

ReduceEnergyMatrix::ReduceEnergyMatrix(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ReduceEnergyMatrix::~ReduceEnergyMatrix() {
}

// FIXME: properly interface with distributed versions of LAPACK
extern "C" {
  void zheev_(
    const char *job, const char *uplo, int *n, complex *a, int *lda,
    double *w, complex *work, int *workCount, double *realWork, int *info
  );
};

void ReduceEnergyMatrix::run() {
  EGH = getTensorArgument<complex>("EnergyMatrix");
  shift = 0.0;
//  preconditionEnergyMatrix();
  readEnergyMatrix();
  diagonalizeEnergyMatrix();
  truncateUnitaryTransform();
  writeUnitaryTransform();
  if (isArgumentGiven("EnergySpectrum")) writeEnergySpectrum();
}

void ReduceEnergyMatrix::dryRun() {
  DryTensor<complex> *EGH(
    getTensorArgument<complex, DryTensor<complex>>("EnergyMatrix")
  );
  nG = EGH->lens[0];


  // calculate number of fieldVariables
  int64_t fieldVariables(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
  // if fieldVariables not given use reduction
  if (fieldVariables == -1) {
    double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
    fieldVariables = nG * reduction;
  }

  LOG(1, "GREDUCE") << "diagonalizing " <<
    nG << "x" << nG << " energy matrix..." << std::endl;

  LOG(1, "GREDUCE") <<
    "taking " << fieldVariables << " of " << nG << " eigenvectors" << std::endl;

  DryMatrix<complex> *UGg(new DryMatrix<complex>(nG, fieldVariables, NS));
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "EnergyMatrixTransform", UGg
  );
}


void ReduceEnergyMatrix::preconditionEnergyMatrix() {
  // evaluate trace
  Scalar<complex> e(*EGH->wrld);
  e[""] = (*EGH)["GG"];
  // shift by -Tr{E} to ensure that each eigenvalue is positive and away from 0
  (*EGH)["GG"] -= e[""];
  // remember shift to reconstruct the spectrum
  shift = std::real(e.get_val());
}

void ReduceEnergyMatrix::readEnergyMatrix() {
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

void ReduceEnergyMatrix::diagonalizeEnergyMatrix() {
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

void ReduceEnergyMatrix::truncateUnitaryTransform() {
  if (EGH->wrld->rank == 0) {
    // allocate elements of transform. In principle, it can be as large as E
    transformElements = new complex[elementsCount];
    ng = 0;

    double energy(0.0);
    for (int i(0); i < nG; ++i) {
      eigenValues[i] += shift;
      energy += eigenValues[i];
    }
    LOG(1, "GREDUCE") << "sum of eigenvalues=" << energy << std::endl;

    int bottom(0), top(nG-1), column(0);
    // approximated energy by truncation
    double e(0);

    // calculate number of fieldVariables
    int64_t fieldVariables(getIntegerArgument("fieldVariables", DEFAULT_FIELD_VARIABLES));
    // if fieldVariables not given use reduction
    if (fieldVariables == -1) {
      double reduction(getRealArgument("reduction", DEFAULT_REDUCTION));
      fieldVariables = nG * reduction;
    }

    //    while (bottom < top && ng < nG * reduction) {

    while (bottom < top && ng < fieldVariables) {
      if (std::abs(eigenValues[bottom]) > std::abs(eigenValues[top])) {
        // bottom value is larger in magnitude, take bottom
        column = bottom++;
      } else {
        // otherwise, take top
        column = top--;
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
      "taking " << ng << " of " << nG << " eigenvectors for a fit accuracy of " <<
      std::abs(e-energy)/energy << std::endl;
    LOG(1, "GREDUCE") <<
      "|smallest/largest eigenvalue|=" <<
      std::abs(eigenValues[column] / eigenValues[0]) << std::endl;
  } else {
    transformElements = new complex[0];
    ng = 0;
  }
  // broadcast ng to all ranks
  MPI_Bcast(&ng, 1, MPI_INT, 0, EGH->wrld->comm);

  delete[] elements;
}

void ReduceEnergyMatrix::writeUnitaryTransform() {
  Matrix<complex> *UGg(new Matrix<complex>(nG, ng, *EGH->wrld, "UGg"));
  int localNg(UGg->wrld->rank == 0 ? ng : 0);
  UGg->write(localNg*nG, indices, transformElements);
  allocatedTensorArgument<complex>(
    "EnergyMatrixTransform", UGg
  );
  delete[] transformElements;
  delete[] indices;
}

// TODO: maybe also write the reduced energy spectrum
// or sort them accordingly 
void ReduceEnergyMatrix::writeEnergySpectrum() {
  Vector<> *EG(new Vector<>(nG, *EGH->wrld, "EG"));
  int localNG(EGH->wrld->rank == 0 ? nG : 0);
  int64_t *eigenValueIndices(new int64_t[localNG]);
  for (int i(0); i < localNG; ++i) {
    eigenValueIndices[i] = i;
  }
  EG->write(localNG, eigenValueIndices, eigenValues);
  allocatedTensorArgument<>("EnergySpectrum", EG);
  delete[] eigenValueIndices;
  delete[] eigenValues;
}

