#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <complex>
#include <array>
#include <algorithm>
#include <fstream>
#include <map>
// this code produces a CoulombVertex for cc4s
// for the UEG for a single cubic lattice
// Choose the density (rs), and number of occupied (No)
// and virutal (Nv) orbitals.
// atomic units are used throughout!
// the size of the vertex is in principle of size Ng x Np x Np
// with Ng = Np x Np. However, using a smaller Ng is a good approximation

// default is the HF reference. However, free electron gas is also possible
bool lhfref = true;


// we have an integer mesh (as it is easier to work with int)
// and a double mesh (with possible shifts) and the HF/kin. energy as 4th entry
using namespace std;
using ivec  = array<int,3>;
using dvec  = array<double,4>;


double evalMadelung(const double v){
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
int sL(const ivec a)   {  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
double sL(const dvec a){  return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];}
double Vijji(const dvec a, const dvec b, const double v){
  dvec q({a[0]-b[0], a[1]-b[1], a[2]-b[2]});
  if ( sL(q) < 1e-8 ) return evalMadelung(v);
  return 4.0*M_PI/v/sL(q);
}


int main(int argc, char ** argv){

  if (argc != 4) {
    printf("Usage: $exe rs No Nv\n");
    return 1;
  }
  double rs(atof(argv[1]));
  int No(atoi(argv[2]));
  int Nv(atoi(argv[3]));
  size_t Np(No+Nv);
  std::vector<int> shell = { 1, 7, 19, 27, 33, 57, 81, 93, 123, 147, 171,\
                             179, 203, 257, 305, 341, 365, 389, 437, 461,\
                             485, 515, 587, 619};

if(No < 620 && std::find(shell.begin(), shell.end(), No) == shell.end()){
    printf("No not valid!\n Possible candidates: ");
    for (auto s: shell) std::cout << s << " ";
    printf("\n");
    return 1;
  }

  if(Np < 620 && std::find(shell.begin(), shell.end(), Np) == shell.end()){
    printf("Nv not valid!\n Possible candidates: ");
    for (auto s: shell) if (s-No > 0) std::cout << s-No << " ";
    printf("\n");
    return 1;
  }

  // setup the integer Grid.
  //  1) gather more than enough candidates
  //  2.) sort by length
  //  3.) split and cut
  int maxG = pow(5.0*Np,1.0/3.0);
  vector<ivec> iGrid;
  for (int g1(-maxG); g1 <= maxG; g1++)
  for (int g2(-maxG); g2 <= maxG; g2++)
  for (int g3(-maxG); g3 <= maxG; g3++)
    iGrid.push_back({g1, g2, g3});

  sort(iGrid.begin(), iGrid.end(), [](ivec a, ivec b){ return sL(a) < sL(b); });
  if (iGrid.size() < Np ) throw invalid_argument("BUG related to Np & maxG\n");
  if (sL(iGrid[No]) == sL(iGrid[No-1])) throw invalid_argument("No not valid\n");
  if (sL(iGrid[Np]) == sL(iGrid[Np-1])) throw invalid_argument("Nv not valid\n");
  iGrid.resize(Np);

  // define volume, lattice Constant, and reciprocal lattice constant
  double v(rs*rs*rs/3.0*4.0*M_PI*No*2);
  double a(pow(v,1./3.));
  double b(2.0*M_PI/a);

  vector<dvec> dGrid;
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
  printf("System Information:\n");
  printf("  rs %5.3lf, No %d, Nv %d\n", rs, No, Nv);
  printf("  Volume %14.10lf, madelung %14.10lf\n", v, evalMadelung(v));
  printf("  HOMO %14.10lf, LUMO %14.10lf\n", dGrid[No-1][3], dGrid[No][3]);
  printf("  Reference Energy per Electron/total %14.10lf / %14.10lf\n",
         refE/No/2, refE);
  cout << endl;


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
            , [](int a, int b){ return abs(a) < abs(b);});
  maxG = abs(maxG);
  map<ivec,int> momMap;
  int index(0);
  for (int g1(-maxG); g1 <= maxG; g1++)
  for (int g2(-maxG); g2 <= maxG; g2++)
  for (int g3(-maxG); g3 <= maxG; g3++){
    ivec t({g1,g2,g3});
    if ( sL(t) > maxR ) continue;
    momMap[t] = index++;
  }


  printf("Evaluate CoulombVertex ...\n");

  size_t NF = momMap.size();
  double fac(4.0*M_PI/v);

  // WRITING COMPLEX VERTEX TO FILE
  vector< complex<double> > out(NF*Np*Np,{0,0});

  for (size_t p(0); p < Np; p++)
  for (size_t q(0); q < Np; q++){
    ivec d = { iGrid[q][0] - iGrid[p][0]
             , iGrid[q][1] - iGrid[p][1]
             , iGrid[q][2] - iGrid[p][2]
             };
    size_t idx = momMap[d];
    double res;
    (sL(d)) ? res = fac/( b*b*sL(d) ) : res = evalMadelung(v);
    out[idx+p*NF+q*NF*Np] = { sqrt(res), 0.0};
  }

  string yamlout;
  yamlout += "version: 100\nscalarType: Complex64\ntype: Tensor\n";
  yamlout += "dimensions:\n  - length:     ";
  yamlout += to_string(NF);
  yamlout += "\n    type: AuxiliaryField\n  - length:     ";
  yamlout += to_string(Np);
  yamlout += "\n    type: State\n  - length:     ";
  yamlout += to_string(Np);
  yamlout += "\n    type: State\nelements:\n  type: IeeeBinaryFile\n";
  yamlout += "metaData:\n  halfGrid: 0\nunit: 1\n";

  ofstream vertyaml;
  vertyaml.open("CoulombVertex.yaml");
  vertyaml << yamlout;
  vertyaml.close();

  ofstream vertex;
  vertex.open("CoulombVertex.elements", ios::out | ios::binary);
  vertex.write((char*)&out[0], NF*Np*Np*sizeof(complex<double>));
  vertex.close();


  // Write Eigenenergies
  double fermiEnergy((dGrid[No][3]+dGrid[No-1][3])/2.0);

  printf("Writing EigenEnergies to file\n");
  char buf[50];
  yamlout.clear();
  yamlout += "version: 100\nscalarType: Real64\n";
  yamlout += "type: Tensor\n";
  yamlout += "dimensions:\n  - length:   ";
  yamlout += to_string(Np);
  yamlout += "\n    type: State\nelements:\n  type: TextFile\n";
  yamlout += "unit: 1\nmetaData:\n  fermiEnergy:  ";
  sprintf(buf, "%16.14lf\n", fermiEnergy);
  yamlout += buf;
  yamlout += "  energies:";
  for ( auto d: dGrid){
    sprintf(buf, "\n    -  %16.14lf", d[3]);
    yamlout += buf;
  }
  yamlout += "\n";

  // Write to yaml
  ofstream eigyaml;
  eigyaml.open("EigenEnergies.yaml");
  eigyaml << yamlout;
  eigyaml.close();

  ofstream eigdat;
  eigdat.open("EigenEnergies.elements");
  yamlout.clear();
  for ( auto d: dGrid){
    sprintf(buf, "%16.14lf\n", d[3]);
    yamlout += buf;
  }
  eigdat << yamlout;
  eigdat.close();


  return 0;
}
