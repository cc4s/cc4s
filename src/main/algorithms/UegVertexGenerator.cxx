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

Real<> UegVertexGenerator::evalMadelung(const Real<> v){
  Real<> kappa = pow(v,-1.0/3.0);
  Real<> term2 = Pi<>() / (kappa*kappa*v);
  Real<> term4 = 2 * kappa/sqrt(Pi<>());
  Real<> boxLength = 1.0/kappa;
  Real<> recipsum = 0.0;
  Real<> realsum = 0.0;
  for (Integer<> l1=-6; l1 <= 6; ++l1)
  for (Integer<> l2=-6; l2 <= 6; ++l2)
  for (Integer<> l3=-6; l3 <= 6; ++l3){
    Integer<> n2 = l1*l1 + l2*l2 + l3*l3;
    Real<> modr = boxLength * sqrt((Real<>)n2);
    Real<> k2 = kappa*kappa*n2;
    if (n2 > 0){
     recipsum -= 1.0/(Pi<>()*k2)*exp(-Pi<>()*Pi<>()*k2/kappa/kappa)/v;
     realsum -= erfc(kappa*modr)/modr;
    }
  }
  return realsum + term2 + term4 + recipsum;
}

//define two functions which give the squared length of the grid-points
Natural<> sL(const ivec a)   {  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
Real<> sL(const dvec a){  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
Real<> UegVertexGenerator::Vijji(const dvec a, const dvec b, const Real<> v){
  dvec q({a[0]-b[0], a[1]-b[1], a[2]-b[2]});
  if ( sL(q) < 1e-8 ) return madelung;
  return 4.0*Pi<>()/v/sL(q);
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
  // We use the HF reference by default.
  bool lhfref(arguments->getValue<bool>("hartreeFock", 1));
  bool lclosed(true);
  No = arguments->getValue<Natural<>>("No");
  Nv = arguments->getValue<Natural<>>("Nv");
  rs = arguments->getValue<Real<>> ("rs");
  NF = arguments->getValue<Natural<>>("NF",0);
  halfGrid = arguments->getValue<bool>("halfGrid",0);
  madelung = arguments->getValue<Real<>> ("madelung", -1.0);
  Natural<> Np(No+Nv);
  if (!No) THROW("No larger zero please");
  if (rs <= 0.0) THROW("Invalid rs");

  // setup the integer Grid.
  //  1) gather more than enough candidates
  //  2.) sort by length
  //  3.) split and cut
  Integer<> maxG = pow(5.0*Np,1.0/3.0);
  std::vector<ivec> iGrid;
  for (Integer<> g1(-maxG); g1 <= maxG; g1++)
  for (Integer<> g2(-maxG); g2 <= maxG; g2++)
  for (Integer<> g3(-maxG); g3 <= maxG; g3++)
    iGrid.push_back({g1, g2, g3});

  sort(iGrid.begin(), iGrid.end(), [](ivec a, ivec b){ return sL(a) < sL(b); });
  if (iGrid.size() < Np ) THROW("BUG related to Np & maxG\n");
  if (sL(iGrid[No]) == sL(iGrid[No-1])){
    OUT() << "WARNING: occupied orbitals form not a closed shell\n";
    if (!lhfref) 
      THROW("Zero gap system! Either change No or use Hartree-Fock!");
    lclosed = false;
  }
  if (sL(iGrid[Np]) == sL(iGrid[Np-1]))
    OUT() << "WARNING: virtual orbitals form not a closed shell\n";
  iGrid.resize(Np);

  // define volume, lattice Constant, and reciprocal lattice constant
  Real<> v(rs*rs*rs/3.0*4.0*Pi<>()*No*2);
  Real<> a(pow(v,1./3.));
  Real<> b(2.0*Pi<>()/a);

  if (madelung < 0.0) madelung = evalMadelung(v);

  std::vector<dvec> dGrid;
  // here we can introduce a possible shift of the mesh
  for (auto i: iGrid)
    dGrid.push_back( { b*i[0], b*i[1], b*i[2], 0.0} );

  // now we can write the hartree fock energy in the 4th entry
  for (auto &d: dGrid){
    d[3] = 0.5*sL(d); // add the kinetic energy
    Real<> exchE(0.0);
    for (Natural<> o(0); o < No; o++)
      exchE += Vijji(d, dGrid[o], v);
    if (lhfref) d[3] -= exchE;
  }
  Real<> refE(0.0);
  for (Natural<> o(0); o < No; o++) {
    refE += dGrid[o][3];
    if (lhfref) refE += 0.5*sL(dGrid[o]);
  }


  std::vector<Real<>> energies(dGrid.size());

  for (Natural<> d(0); d < dGrid.size(); d++)
    energies[d] = dGrid[d][3];

  // Prepare eigenEnergies
  auto eigenEnergies(
    Tcc<TE>::template tensor<Real<>>({Np}, "EigenEnergies")
   );

  Real<> fermiEnergy((energies[No] + energies[No-1])/2.0);

  auto metaData( New<MapNode>(SOURCE_LOCATION) );
  metaData->setValue("fermiEnergy", fermiEnergy);
  auto energiesNode( New<MapNode>(SOURCE_LOCATION) );
  for (Natural<> i(0); i < Np; ++i) {
    energiesNode->setValue(i, energies[i]);
  }
  metaData->get("energies") = energiesNode;
  eigenEnergies->getMetaData() = metaData;

  // construct the momentum transition grid
  // 1.) get the largest momentum vector between two states p - q
  // 2.) construct a full grid with a largest grid vec. of this size

  ivec maxMom({0,0,0});
  for (Natural<> p(0); p < Np; p++)
  for (Natural<> q(0); q < Np; q++){
    ivec d = { iGrid[p][0] - iGrid[q][0]
             , iGrid[p][1] - iGrid[q][1]
             , iGrid[p][2] - iGrid[q][2]
             };
    maxMom = max(maxMom, d, [](ivec a, ivec b) { return sL(a) < sL(b);});
  }

  Natural<> maxR = sL(maxMom);
  maxG = max( {maxMom[0], maxMom[1], maxMom[2]}
            , [](Integer<> a, Integer<> b){ return std::abs(a) < std::abs(b);});
  maxG = std::abs(maxG);
  std::map<ivec,Natural<>> momMap;
  Natural<> index(0);
  for (Integer<> g1(-maxG); g1 <= maxG; g1++)
  for (Integer<> g2(-maxG); g2 <= maxG; g2++)
  for (Integer<> g3(-maxG); g3 <= maxG; g3++){
    ivec t({g1,g2,g3});
    if ( sL(t) > maxR ) continue;
    momMap[t] = index++;
  }
  if (NF == 0) NF = momMap.size();

  if (NF != momMap.size() || halfGrid || !lclosed)
    OUT() << "WARNING: the Vertex will not be correct! Just for profiling!\n";

  Real<> fac(4.0*Pi<>()/v);

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
  std::vector<Natural<>> idx;
  if (!Cc4s::world->getRank()){
    idx.resize(energies.size());
    std::iota(idx.begin(), idx.end(), 0);
  }
  eigenEnergies->write(idx.size(), idx.data(), energies.data());



  // Writing CoulombVertex to buffer
  // We have to do it mpi-able...otherwise we will
  // not be able to write it to a ctf tensor
  Natural<> np = Cc4s::world->getProcesses();
  Natural<> rank = Cc4s::world->getRank();
  // We slice the number of states for all the mpi processes
  Natural<> slices(Np/np);
  std::vector<Natural<>> slicePerRank(np);
  for (Natural<> r(0); r < np; r++){
    Natural<> lslice(slices);
    for (Natural<> i(0); i < Np - slices*np; i++) if (r == i){
      lslice++;
    }
    slicePerRank[r] = lslice;
  }
  slices = slicePerRank[rank];
  //allocate only a buffer of needed size
  std::vector< Complex<> > out(NF*Np*slices,{0,0});
  // determine begin and end of the rank's slices
  auto sbegin( std::accumulate( slicePerRank.begin()
                              , slicePerRank.begin() + rank
                              , 0UL
                              , std::plus<Natural<>>()
                              )
             );

  for (Natural<> s(0); s < slices; s++)
  for (Natural<> q(0); q < Np; q++){
    auto p(s+sbegin);
    ivec d = { iGrid[q][0] - iGrid[p][0]
             , iGrid[q][1] - iGrid[p][1]
             , iGrid[q][2] - iGrid[p][2]
             };
    // This is a hack!
    // If NF is chosen by the user we will not have an overflow
    Natural<> ii = momMap[d] % NF;
    Real<> res;
    (sL(d)) ? res = fac/( b*b*sL(d) ) : res = evalMadelung(v);
    out[ii+q*NF+s*NF*Np] = { sqrt(res), 0.0};
  }

  idx.resize(out.size());
  std::iota(idx.begin(), idx.end(), sbegin*Np*NF);
  coulombVertex->write(idx.size(), idx.data(), out.data());


  return result;
}

