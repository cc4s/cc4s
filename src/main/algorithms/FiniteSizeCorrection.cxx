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

#include <algorithms/FiniteSizeCorrection.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <gte/TricubicInterpolation.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(FiniteSizeCorrection)

Ptr<MapNode> FiniteSizeCorrection::run(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));

  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    interpolation<TE>(arguments, result);
  }
  else {
    using TE = DefaultTensorEngine;
    interpolation<TE>(arguments, result);
  } 
  return result;
}


template <typename TE>
void FiniteSizeCorrection::interpolation(
  const Ptr<MapNode> &arguments,
  Ptr<MapNode> &result
) {
  using T = Tensor<Real<>, TE>;

  double sum3D(0.0), inter3D(0.0);
  int64_t countNO(0), countNOg(0);
  // hard coded resolution of the fine grid
  int N(20);
  // this is hard coded for vasp:
  // the coulomb Energy in eV using the lattice dimensions of vasp

  auto gridVectors(arguments->getMap("gridVectors"));
  auto volume(gridVectors->getValue<Real<>>("volume"));

// constant factor is still unclear to me
//  double factor(4.5835494674469/volume);
  double factor(0.00020880076563/volume);
  // READ THE MOMENTUM GRID
  auto grid(gridVectors->getValue<Ptr<T>>("data"));
  ASSERT_LOCATION(
    grid, "expecting the reciprocal Grid",
    gridVectors->sourceLocation
  );

  int64_t NG(grid->lens[1]);
  std::vector<Vector<>> cartesianGrid(NG);
  std::vector<Real<>> output(NG*3);
  output = grid->readAll();
  for (int64_t i(0); i < NG; i++){
    cartesianGrid[i][0] = output[3*i+0];
    cartesianGrid[i][1] = output[3*i+1];
    cartesianGrid[i][2] = output[3*i+2];
  }


  // READ THE RECIPROCAL CELL
  std::vector<Vector<>> B(3);
  auto Gix(gridVectors->getValue<Real<>>("Gix"));
  auto Giy(gridVectors->getValue<Real<>>("Giy"));
  auto Giz(gridVectors->getValue<Real<>>("Giz"));
  auto Gjx(gridVectors->getValue<Real<>>("Gjx"));
  auto Gjy(gridVectors->getValue<Real<>>("Gjy"));
  auto Gjz(gridVectors->getValue<Real<>>("Gjz"));
  auto Gkx(gridVectors->getValue<Real<>>("Gkx"));
  auto Gky(gridVectors->getValue<Real<>>("Gky"));
  auto Gkz(gridVectors->getValue<Real<>>("Gkz"));


  B[0][0] = Gix; B[0][1] = Giy; B[0][2] = Giz ;
  B[1][0] = Gjx; B[1][1] = Gjy; B[1][2] = Gjz ;
  B[2][0] = Gkx; B[2][1] = Gky; B[2][2] = Gkz ;

  // READ THE Structure Factor
  std::vector<Real<>> SofG;
  auto structureFactor(arguments->getValue<Ptr<T>>("structureFactor"));
//  ASSERT_LOCATION(
//    structureFactor, "expecting the structureFactor",
//    structureFactor->sourceLocation
//  );
  SofG = structureFactor->readAll();

  // the tricubic interpolation requires a rectangular grid
  // ---> transformation

  //construct the transformation matrix, which is the real cell
  std::vector<Vector<>> A(3);
  double Omega((B[0].cross(B[1])).dot(B[2]));
  A[0] = B[1].cross(B[2])/Omega;
  A[1] = B[2].cross(B[0])/Omega;
  A[2] = B[0].cross(B[1])/Omega;


  // determine bounding box in direct coordinates (in reciprocal space)
  Vector<> directMin, directMax;
  for (int g(0); g < NG; ++g)
  for (int d(0); d < 3; ++d) {
    double directComponent(A[d].dot(cartesianGrid[g]));
    directMin[d] = std::min(directMin[d], directComponent);
    directMax[d] = std::max(directMax[d], directComponent);
  }
  // build grid for the entire bounding box
  Vector<> boxDimensions, boxOrigin;
  int64_t boxSize(1);
  for (int d(0); d < 3; ++d) {
    boxSize *= boxDimensions[d] = std::floor(directMax[d] - directMin[d] + 1.5);
    boxOrigin[d] = std::floor(directMin[d] + 0.5);
  }

  // allocate and initialize direct-grid StructureFactor
  std::vector<Real<>> directSofG(boxSize, 0.0);

  OUT() << boxSize << std::endl;
  // enter known SG values
  for (int g(0); g < NG; ++g) {
    int64_t index(0);
    Vector<> directG;
    for (int d(2); d >= 0; --d) {
      directG[d] = A[d].dot(cartesianGrid[g]);
      index *= boxDimensions[d];
      index += std::floor(directG[d] + 0.5) - boxOrigin[d];
    }
    directSofG[index] = SofG[g];
  }

  // create trilinear or tricubic interpolator
  gte::IntpTricubic3<Real<>> interpolatedSofG(
    (int)boxDimensions[0], (int)boxDimensions[1], (int)boxDimensions[2],
    boxOrigin[0], 1.0, boxOrigin[1], 1.0, boxOrigin[2], 1.0,
    directSofG.data(),
    true
  );

  // double check: calculate the structure Factor on the regular grid
  for (int64_t i(0); i < NG; ++i) {
    if (cartesianGrid[i].length() < 1e-8) continue;
    sum3D += factor/cartesianGrid[i].sqrLength()*SofG[i];
  }


  MpiCommunicator communicator(
    Cc4s::world->getRank(), Cc4s::world->getProcesses(), Cc4s::world->getComm()
  );
  communicator.barrier();
  for (int a(0); a < N; ++a)
  for (int b(0); b < N; ++b)
  for (int c(0); c < N; ++c) {
    if ((std::abs(a+b+c)) % Cc4s::world->getProcesses() != Cc4s::world->getRank()) continue;
    Vector<double> ga(((B[0]/double(N))*double(a)));
    Vector<double> gb(((B[1]/double(N))*double(b)));
    Vector<double> gc(((B[2]/double(N))*double(c)));
    Vector<double> g(ga+gb+gc);
    countNOg++;
    // loop over coarse grid. i.e. shift fine grid onto every coarse grid-point
    for ( size_t i(0); i < NG; i++){
      g  = ga+gb+gc;
      g += cartesianGrid[i];
      Vector<double> directG;
      for (size_t d(0); d < 3; d++) directG[d] = A[d].dot(g);
      countNO++;
      if (g.length() < 1e-8) continue;
      double interpol(interpolatedSofG(directG[0], directG[1], directG[2]));
      inter3D += interpol * factor/g.sqrLength();
    }
  }
  double totalInter3D(0.0);
  communicator.allReduce(inter3D, totalInter3D);


  OUT() << "uncorrected: " << sum3D << "\n";
  OUT() << "corrected:   " << totalInter3D/N/N/N << "\n";


  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("corrected", inter3D);
  energy->setValue<Real<>>("uncorrected", sum3D);
  result->get("energy") = energy;
}
