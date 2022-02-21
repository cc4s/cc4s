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

#include <algorithms/UegVertexGenerator.hpp>

#include <tcc/Tcc.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

double UegVertexGenerator::evalMadelung(const double v){
  double kappa = pow(v,-1.0/3.0);
  double term2 = M_PI / (kappa*kappa*v);
  double term4 = 2 * kappa/sqrt(M_PI);
  double boxLength = 1.0/kappa;
  double recipsum = 0.0;
  double realsum = 0.0;
  for (int l1=-6; l1 <= 6; ++l1)
  for (int l2=-6; l2 <= 6; ++l2)
  for (int l3=-6; l3 <= 6; ++l3){
    int n2 = l1*l1 + l2*l2 + l3*l3;
    double modr = boxLength * sqrt((double)n2);
    double k2 = kappa*kappa*n2;
    if (n2 > 0){
     recipsum -= 1.0/(M_PI*k2)*exp(-M_PI*M_PI*k2/kappa/kappa)/v;
     realsum -= erfc(kappa*modr)/modr;
    }
  }
  return realsum + term2 + term4 + recipsum;
};

//define two functions which give the squared length of the grid-points
size_t sL(const ivec a)   {  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
double sL(const dvec a){  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
double UegVertexGenerator::Vijji(const dvec a, const dvec b, const double v){
  dvec q({a[0]-b[0], a[1]-b[1], a[2]-b[2]});
  if ( sL(q) < 1e-8 ) return madelung;
  return 4.0*M_PI/v/sL(q);
}



ALGORITHM_REGISTRAR_DEFINITION(UegVertexGenerator)

Ptr<MapNode> UegVertexGenerator::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (result = run<Real<>,TE>(arguments))
      || (result = run<Complex<>,TE>(arguments));
  } else {
    using TE = DefaultTensorEngine;
    (result = run<Real<>,TE>(arguments))
      || (result = run<Complex<>,TE>(arguments));
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> UegVertexGenerator::run(
  const Ptr<MapNode> &arguments
) {
  // We always use the HF reference.
  bool lhfref(true);
  No = arguments->getValue<int>("No");
  Nv = arguments->getValue<int>("Nv");
  rs = arguments->getValue<double> ("rs");
  NF = arguments->getValue<int>("NF",0);
  halfGrid = arguments->getValue<int>("halfGrid",0);
  madelung = arguments->getValue<double> ("madelung", -1.0);
  size_t Np(No+Nv);
  if (!No) THROW("No larger zero please");
  if (rs <= 0.0) THROW("Invalid rs");

  // setup the integer Grid.
  //  1) gather more than enough candidates
  //  2.) sort by length
  //  3.) split and cut
  int maxG = pow(5.0*Np,1.0/3.0);
  std::vector<ivec> iGrid;
  for (int g1(-maxG); g1 <= maxG; g1++)
  for (int g2(-maxG); g2 <= maxG; g2++)
  for (int g3(-maxG); g3 <= maxG; g3++)
    iGrid.push_back({g1, g2, g3});

  sort(iGrid.begin(), iGrid.end(), [](ivec a, ivec b){ return sL(a) < sL(b); });
  if (iGrid.size() < Np ) THROW("BUG related to Np & maxG\n");
  if (sL(iGrid[No]) == sL(iGrid[No-1])) THROW("No not valid\n");
  if (sL(iGrid[Np]) == sL(iGrid[Np-1])) THROW("Nv not valid\n");
  iGrid.resize(Np);

  // define volume, lattice Constant, and reciprocal lattice constant
  double v(rs*rs*rs/3.0*4.0*M_PI*No*2);
  double a(pow(v,1./3.));
  double b(2.0*M_PI/a);

  if (madelung < 0.0) madelung = evalMadelung(v);

  std::vector<dvec> dGrid;
  // here we can introduce a possible shift of the mesh
  for (auto i: iGrid)
    dGrid.push_back( { b*i[0], b*i[1], b*i[2], 0.0} );

  // now we can write the hartree fock energy in the 4th entry
  for (auto &d: dGrid){
    d[3] = 0.5*sL(d); // add the kinetic energy
    double exchE(0.0);
    for (size_t o(0); o < No; o++)
      exchE += Vijji(d, dGrid[o], v);
    if (lhfref) d[3] -= exchE;
  }
  double refE(0.0);
  for (size_t o(0); o < No; o++) {
    refE += dGrid[o][3];
    if (lhfref) refE += 0.5*sL(dGrid[o]);
  }


  std::vector<double> energies(dGrid.size());

  for (auto d(0); d < dGrid.size(); d++)
    energies[d] = dGrid[d][3];

  // Prepare eigenEnergies
  auto eigenEnergies(
    Tcc<TE>::template tensor<Real<>>({Np}, "EigenEnergies")
   );

  double fermiEnergy((energies[No] + energies[No-1])/2.0);

  auto metaData( New<MapNode>(SOURCE_LOCATION) );
  metaData->setValue("fermiEnergy", fermiEnergy);
  auto energiesNode( New<MapNode>(SOURCE_LOCATION) );
  for (size_t i(0); i < Np; ++i) {
    energiesNode->setValue(i, energies[i]);
  }
  metaData->get("energies") = energiesNode;
  eigenEnergies->getMetaData() = metaData;

  // construct the momentum transition grid
  // 1.) get the largest momentum vector between two states p - q
  // 2.) construct a full grid with a largest grid vec. of this size

  ivec maxMom({0,0,0});
  for (int p(0); p < Np; p++)
  for (int q(0); q < Np; q++){
    ivec d = { iGrid[p][0] - iGrid[q][0]
             , iGrid[p][1] - iGrid[q][1]
             , iGrid[p][2] - iGrid[q][2]
             };
    maxMom = max(maxMom, d, [](ivec a, ivec b) { return sL(a) < sL(b);});
  }

  int maxR = sL(maxMom);
  maxG = max( {maxMom[0], maxMom[1], maxMom[2]}
            , [](int a, int b){ return std::abs(a) < std::abs(b);});
  maxG = std::abs(maxG);
  std::map<ivec,int> momMap;
  int index(0);
  for (int g1(-maxG); g1 <= maxG; g1++)
  for (int g2(-maxG); g2 <= maxG; g2++)
  for (int g3(-maxG); g3 <= maxG; g3++){
    ivec t({g1,g2,g3});
    if ( sL(t) > maxR ) continue;
    momMap[t] = index++;
  }
  if (NF == 0) NF = momMap.size();

  if (NF != momMap.size() || halfGrid)
    OUT() << "WARNING: the Vertex will not be correct! Just for profiling!\n";

  double fac(4.0*M_PI/v);

  auto coulombVertex(
    Tcc<TE>::template tensor<Complex<>>({NF,Np,Np}, "CoulombVertex")
   );

  auto metaCoulomb( New<MapNode>(SOURCE_LOCATION) );
  metaCoulomb->setValue("halfGrid", halfGrid);
  coulombVertex->getMetaData() = metaCoulomb;

  OUT() << "System Information:\n";
  OUT() << std::setprecision(3) << "  rs " << rs
        << ", No " << No << ", Nv" << Nv << "\n";
  OUT() << std::setprecision(10) << "  Volume " << v
        << ", madelung " << madelung << "\n";
  OUT() << "  HOMO " << energies[No-1] << ", LUMO " << energies[No] << "\n";
  OUT() << "  Reference Energy per Electron/total "
        << refE/No/2 << "/" << refE << std::endl;

  // manually enter tensor dimension info in tcc
  auto stateDimension(New<TensorDimension>());
  stateDimension->name = "State";
  cc4s::TensorDimension::dimensions["State"] = stateDimension;
  auto auxFieldDimension(New<TensorDimension>());
  auxFieldDimension->name = "AuxiliaryField";
  cc4s::TensorDimension::dimensions["AuxiliaryField"] = auxFieldDimension;

  // manually enter tensor dimension in tensors
  eigenEnergies->dimensions = std::vector<Ptr<TensorDimension>>({stateDimension});
  coulombVertex->dimensions = std::vector<Ptr<TensorDimension>>(
    {auxFieldDimension, stateDimension, stateDimension}
  );

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setPtr("eigenEnergies", eigenEnergies);
  result->setPtr("coulombVertex", coulombVertex);

  // if we are in a dryRun there is nothing more to do
  if (Cc4s::dryRun) return result;

  //only rank 0 writes the data to the tensor
  std::vector<size_t> idx;
  if (!Cc4s::world->getRank()){
    idx.resize(energies.size());
    std::iota(idx.begin(), idx.end(), 0);
  }
  eigenEnergies->write(idx.size(), idx.data(), energies.data());



  // Writing CoulombVertex to buffer
  // We have to do it mpi-able...otherwise we will
  // not be able to write it to a ctf tensor
  size_t np = Cc4s::world->getProcesses();
  size_t rank = Cc4s::world->getRank();
  // We slice the number of states for all the mpi processes
  size_t slices(Np/np);
  std::vector<size_t> slicePerRank(np);
  for (size_t r(0); r < np; r++){
    size_t lslice(slices);
    for (size_t i(0); i < Np - slices*np; i++) if (r == i){
      lslice++;
    }
    slicePerRank[r] = lslice;
  }
  slices = slicePerRank[rank];
  //allocate only a buffer of needed size
  std::vector< complex<double> > out(NF*Np*slices,{0,0});
  // determine begin and end of the rank's slices
  auto sbegin( std::accumulate( slicePerRank.begin()
                              , slicePerRank.begin() + rank
                              , 0UL
                              , std::plus<size_t>()
                              )
             );
  auto send(sbegin+slices);

  for (size_t s(0); s < slices; s++)
  for (size_t q(0); q < Np; q++){
    auto p(s+sbegin);
    ivec d = { iGrid[q][0] - iGrid[p][0]
             , iGrid[q][1] - iGrid[p][1]
             , iGrid[q][2] - iGrid[p][2]
             };
    // This is a hack!
    // If NF is chosen by the user we will not have an overflow
    size_t ii = momMap[d] % NF;
    double res;
    (sL(d)) ? res = fac/( b*b*sL(d) ) : res = evalMadelung(v);
    out[ii+q*NF+s*NF*Np] = { sqrt(res), 0.0};
  }

  idx.resize(out.size());
  std::iota(idx.begin(), idx.end(), sbegin*Np*NF);
  coulombVertex->write(idx.size(), idx.data(), out.data());


  return result;
}

