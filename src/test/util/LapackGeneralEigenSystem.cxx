#include <test/Test.hpp>

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>

#include <util/Log.hpp>

using namespace cc4s;

TEST_CASE( "LapackGeneralEigenSystem", "[util]" ) {
  LapackMatrix<complex> A(3,3);
  A(0,0) = 0.0; A(0,1) = 1.0; A(0,2) = 0.0;
  A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 1.0;
  A(2,0) = 1.0; A(2,1) = 0.0; A(2,2) = 0.0;
  std::vector<complex> expectedEigenValues(
    {{
      complex(-0.5,+0.8660254037844386), complex(-0.5,-0.8660254037844386), 1.0
    }}
  );

  LapackGeneralEigenSystem<complex> eigenSystem(A);

  REQUIRE(
    std::abs(eigenSystem.getEigenValues()[0] - expectedEigenValues[0]) < 1E-14
  );
  REQUIRE(
    std::abs(eigenSystem.getEigenValues()[1] - expectedEigenValues[1]) < 1E-14
  );
  REQUIRE(
    std::abs(eigenSystem.getEigenValues()[2] - expectedEigenValues[2]) < 1E-14
  );

  for (int k(0); k < eigenSystem.getRightEigenVectors().getColumns(); ++k) {
    std::stringstream rowStream;
    for (int i(0); i < eigenSystem.getRightEigenVectors().getRows(); ++i) {
      rowStream << "\t" << eigenSystem.getRightEigenVectors()(i,k);
    }
    LOG(0, "LapackGeneralEigenSystem") <<
      "EigenVector[" << k << "]=" << rowStream.str() << std::endl;
  };
}
