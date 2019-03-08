#include <algorithms/BasisSetExtrapolation.hpp>

#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <tcc/DryTensor.hpp>
#include <math/Vector.hpp>
#include <util/SharedPointer.hpp>
#include <ctf.hpp>
#include <util/MpiCommunicator.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(BasisSetExtrapolation);

BasisSetExtrapolation::BasisSetExtrapolation(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

BasisSetExtrapolation::~BasisSetExtrapolation() {
}

void BasisSetExtrapolation::run() {

  int fQGG(getIntegerArgument("calculateQGG",0));
  if(fQGG == 1){
    // slice QGG evalution in case of memory bottleneck for the exchange term
    int slice(getIntegerArgument("slice",-1));
    LOG(0,"BasisSetExtrapolation:") << "evaluating QGG" << std::endl;
    evaluateQGG(slice);
  }

  int fFitF12(getIntegerArgument("fitF12",-1));
  if (fFitF12 == 1 || fFitF12 == 2){
    real minG(getRealArgument("minG",-1));
    real maxG(getRealArgument("maxG",-1));
    if ( minG > maxG || minG < 0. ) throw new EXCEPTION("need fitting range:minG and maxG");
    LOG(0,"BasisSetExtrapolation") << "fitting gamma" << std::endl;
    fitF12(fFitF12,minG,maxG);
  }

}



void BasisSetExtrapolation::evaluateQGG(int slice){
  
  PTR(Tensor<complex>) GammaGai;

  if (isArgumentGiven("ParticleHoleCoulombVertex")){
    GammaGai = NEW(
      Tensor<complex>, getTensorArgument<complex>("ParticleHoleCoulombVertex")
    );
  }
  else if (isArgumentGiven("CoulombVertex")){

    if (!isArgumentGiven("HoleEigenEnergies")){
      throw new EXCEPTION("Need HoleEigenEnergies for number of holes/particles");
    }

    auto epsi = NEW( Tensor<>, getTensorArgument<>("HoleEigenEnergies"));
    int No(epsi->lens[0]);
    auto GammaGqr = NEW(
      Tensor<complex>, getTensorArgument<complex>("CoulombVertex")
    );
    int NG(GammaGqr->lens[0]);
    int Np(GammaGqr->lens[1]);
    int Nv(Np-No);
    int aStart(Np-Nv), aEnd(Np);
    int iStart(0), iEnd(No);
    int GaiStart[] = {0, aStart, iStart};
    int GaiEnd[] = {NG, aEnd, iEnd};
    GammaGai = NEW(
      Tensor<complex>, GammaGqr->slice(GaiStart,GaiEnd)
    );
  }
  else {
    throw new EXCEPTION("Need Appropriate Coulomb Vertex");
  }

  // Gamma Gai --> Cai
  int NF(GammaGai->lens[0]);
  int No(GammaGai->lens[2]);
  int Nv(GammaGai->lens[1]);

  Tensor<> *ctfCoulombKernel(getTensorArgument<>("CoulombKernel"));

  int NFF[] = {NF};
  Tensor<complex> invSqrtVG(1, NFF);

  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
    [](real vG, complex &invVG){
      if ( std::isinf(vG) ){
        invVG = 0.;
      }
      else{
        invVG = std::sqrt(1./vG);
      }
    }
   )
  )(
    (*ctfCoulombKernel)["G"],invSqrtVG["G"]
  );

  (*GammaGai)["Gai"] *= invSqrtVG["G"];


  // determine if we are using full or half mesh

  auto momenta(NEW(Tensor<>, getTensorArgument<>("Momenta")));
  cc4s::Vector<> *cartesianMomenta(new cc4s::Vector<>[NF]);
  momenta->read_all(&cartesianMomenta[0][0]);

  cc4s::Vector<> check_grid;
  for (int g(0); g < NF; ++g){
    check_grid+=cartesianMomenta[g];
  }
  PTR(Tensor<complex>) CGai;
  if ( check_grid.length() > 1e-5){
    // maybe the easiest is to double Cia(G).
    // we need to rescale Cia(G)=>Cia(G)/sqrt(2)
    LOG(1,"Build up Q(G,F)") << "working with half mesh" << std::endl;
    int NGai[] = {2*NF-1, Nv, No};
    real invsqrt(1./std::sqrt(2.));
    CGai = NEW(Tensor<complex>,3, NGai);
    // First put all 'positive' G in the full Cia(G)
    CGai->slice(
      std::vector<int>({0,0,0}).data(),std::vector<int>({NF, Nv, No}).data(),1.0,*GammaGai,
      std::vector<int>({0,0,0}).data(),std::vector<int>({NF, Nv, No}).data(),invsqrt
    );
    // Now all 'negative' G
    conjugate(*GammaGai);
    CGai->slice(
      std::vector<int>({NF,0,0}).data(),std::vector<int>({2*NF-1,Nv,No}).data(),1.0,*GammaGai,
      std::vector<int>({ 1,0,0}).data(),std::vector<int>({  NF  ,Nv,No}).data(),invsqrt
    );
    NF = NF*2-1;
  }
  else{
    LOG(1,"Build up Q(G,F)") << "working with full mesh" << std::endl;
    CGai = NEW(Tensor<complex>,*GammaGai);
  }


  auto conjCGai(NEW(Tensor<complex>,*CGai));
  conjugate(*conjCGai);

  int NGG[] = {NF, NF};
  auto QGGs( new CTF::Tensor<complex>(2,NGG));
  auto QGGt( new CTF::Tensor<complex>(2,NGG));

  // Direct part.

  PTR(Tensor<complex>) FGone;
  PTR(Tensor<complex>) FGtwo;
  FGone = NEW(Tensor<complex>,2,NGG);
  FGtwo = NEW(Tensor<complex>,2,NGG);

  (*FGone)["GF"] =  (*conjCGai)["Gai"] * (*CGai)["Fai"];
  (*FGtwo)["GF"] =  (*conjCGai)["Fbj"] * (*CGai)["Gbj"];
  (*QGGs)["GF"]  = (0.5) * (*FGone)["GF"] * (*FGtwo)["GF"];
  (*QGGt)["GF"]  = (1.5) * (*QGGs)["GF"];

  // Exchange part, slicing.

  if (slice < 0 || slice > No) slice = No;

  LOG(0,"Number of bands treated simultaniously:") << slice << std::endl;

  int numberSlices(std::ceil(static_cast<double>(No)/static_cast<double>(slice)));

  for (int ii(0); ii<numberSlices; ++ii){

    int startBandSlice(ii*slice);
    int endBandSlice((ii+1)*slice);
    endBandSlice = std::min(endBandSlice,No);
    int NFG[] = {NF, NF, No, endBandSlice-startBandSlice};
    LOG(0,"Slicing for memory reasons:") << ii+1 << " From: " << startBandSlice
                                                 << " to: " << endBandSlice << std::endl;


    PTR(Tensor<complex>) FGij;
    PTR(Tensor<complex>) FGji;

    FGij = NEW(Tensor<complex>,4,NFG);
    FGji = NEW(Tensor<complex>,4,NFG);


    int CGajStart[] = {0, 0, startBandSlice};
    int CGajEnd[]   = {NF, Nv, endBandSlice};

    auto CGajSliced(NEW(Tensor<complex>, CGai->slice(CGajStart,CGajEnd)));
    auto conjCGajSliced(NEW(Tensor<complex>, conjCGai->slice(CGajStart,CGajEnd)));

    (*FGij)["GFij"] = (*conjCGai)["Gai"] * (*CGajSliced)["Faj"];
    (*FGji)["GFij"] = (*conjCGai)["Fbi"] * (*CGajSliced)["Gbj"];

    (*QGGs)["GF"]  += (0.5)* (*FGij)["GFij"] * (*FGji)["GFij"];
    (*QGGt)["GF"]  += (-0.75)* (*FGij)["GFij"] * (*FGji)["GFij"];


  }

  allocatedTensorArgument<complex>("QGGs",QGGs);
  allocatedTensorArgument<complex>("QGGt",QGGt);

}

void BasisSetExtrapolation::calculateNewSF(
  int type, real gamma, Tensor<> *coulombKernel, Tensor<> *newSF, Tensor<> *resNewSF){

  CTF::Tensor<complex> *QGG(getTensorArgument<complex>("QGG"));
  int NG(QGG->lens[0]);
  int NFF[] = {NG};
  
  auto cK(new Tensor<>(1,NFF));
  (*cK) = (*coulombKernel);

  real volume(getRealArgument("volume",-1.));

  auto reciprocalYC(new Tensor<complex>(1,NFF));

  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [&type, &volume, &gamma](real cK, complex &YC){
        if ( cK == 0 ){
          YC = 0.;
        }
        else{
          if (type ==1){
            YC = 1./cK/cK; //*volume/4.5835494674469; //Ha(eV)*a0(A)*4*pi/(2*pi)**2
  //        rYC = rYC/150.4121/(gamma*gamma+150.4121/rYC); hbar*hbar/2*me*10**20meter*(2pi)**2 in eV
  //        rYC = rYC*4.5835494674469/volume/(gamma*gamma+2.*150.4121/rYC);   // GAMMA[Energy]
            YC *= (-0.015237)/volume/(gamma*gamma+1./YC);               // GAMMA[1/A]
          }
          else if (type ==2){
            YC = (-0.015237)/volume/(gamma*gamma+cK*cK)/(gamma*gamma+cK*cK);
          }
       }
      }
    )
  )(
   (*cK)["G"],(*reciprocalYC)["G"]
  );

  auto dReciprocalYC(new Tensor<complex>(1,NFF));

  CTF::Transform<real, complex, complex >(
    std::function<void(real, complex, complex &)>(
      [&type, &volume, &gamma](real cK, complex YC, complex &dYC){
        if ( cK == 0){
          dYC = 0.; 
        }
        else{
          if(type==1){
            dYC = (-2.)*YC*gamma/(gamma*gamma+cK*cK);
          }
          else if(type==2){
            dYC = (-4.)*YC*gamma/(gamma*gamma+cK*cK);
          }
        }
       }
     ) 
  )(
   (*cK)["G"],(*reciprocalYC)["G"],(*dReciprocalYC)["G"]
  );
  auto complexNewSF(new Tensor<complex>(1,NFF));
  auto complexResNewSF(new Tensor<complex>(1,NFF));

  (*complexNewSF)["G"] = (*QGG)["GF"] * (*reciprocalYC)["F"];
  (*complexResNewSF)["G"] = (*QGG)["GF"] * (*dReciprocalYC)["F"];

  fromComplexTensor(*complexNewSF,*newSF);
  fromComplexTensor(*complexResNewSF,*resNewSF);

}

void BasisSetExtrapolation::fitF12(int type, real minG, real maxG){

  real volume(getRealArgument("volume",-1));
  if (volume < 0. ) throw new EXCEPTION("Set volume");

  CTF::Tensor<> *structureFactor(getTensorArgument<>("StructureFactor"));
  CTF::Tensor<> *coulombKernel(getTensorArgument<>("CoulombKernel"));

  int NG(coulombKernel->lens[0]);
  int NFF[] = {NG};

  auto resNewSF(new Tensor<>(1,NFF));
  auto newSF(new Tensor<>(1,NFF));
  auto absoluteG(new Tensor<>(1,NFF));

  //Take out infinity
  CTF::Transform<real>(
    std::function<void(real &)>(
      [&volume](real &cK){
        if( std::isinf(cK)){
          cK = 0.;
        }
      }
    )
  )(
    (*coulombKernel)["G"]
  );
  // Construct array: aboluteG

  CTF::Transform<real,real>(
    std::function<void(real , real &)>(
      [&volume](real cK, real &absG){
        if( cK == 0. ){
          absG = 0.;
        }
        else{
          absG = cK*volume/4.5835494674469;
          absG = 1./std::sqrt(absG);
        }
      }
    )
  )(
   (*coulombKernel)["G"],(*absoluteG)["G"]
  );

  real gamma(getRealArgument("gamma",1));

  for (int i(0); i<=iterations; ++i){

    calculateNewSF(type,gamma,absoluteG,newSF,resNewSF);

    auto dummy(new Tensor<>(1,NFF));
    int counter(0);
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&counter, &maxG, &minG](real absG, real res, real &dummy){
          if ( absG > maxG || absG < minG ){
            dummy = 0.;
          }
          else{
            dummy = res*res;
            counter++;
          }
        }
     )
    )(
     (*absoluteG)["G"],(*resNewSF)["G"],(*dummy)["G"]
    );

    CTF::Scalar<double> denom;
    denom[""] = (*dummy)["G"];
    
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [](real newsf, real sf, real &dummy){
          dummy = newsf - sf;
        }
      )
    )(
     (*newSF)["G"],(*structureFactor)["G"],(*dummy)["G"]
    );

    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&maxG, &minG](real absG, real resnewsf, real &dummy){
          if ( absG > maxG || absG < minG ){
            dummy = 0.;
          }
          else{
            dummy *= resnewsf;
          }
        }
      )
    )(
     (*absoluteG)["G"],(*resNewSF)["G"],(*dummy)["G"]
    );

    CTF::Scalar<double> numer;
    numer[""] = (*dummy)["G"];
    gamma -= numer.get_val()/denom.get_val();

    (*resNewSF) = (*structureFactor);
    CTF::Transform<real, real, real>(
      std::function<void(real, real, real &)>(
        [&maxG, &minG](real absG, real newsf, real &resnewsf){
          if ( absG > maxG || absG < minG){
            resnewsf = 0.;
          }
          else{
            resnewsf -= newsf;
          }
         }
      )
    )(
     (*absoluteG)["G"],(*newSF)["G"],(*resNewSF)["G"]
    );

    LOG(0,"gamma") << gamma << " norm: " << resNewSF->norm2() << std::endl;

  }
//  allocatedTensorArgument<>("ResNewSF",resNewSF);
  setRealArgument("gammaout",gamma);
  allocatedTensorArgument<>("NewSF",newSF);
// evaluate energy correction term E = v(G)*Q(G,G')*f12(G')
  Scalar<> f12EnergyCorrection(*Cc4s::world);
  f12EnergyCorrection[""] = (*coulombKernel)["G"] * (*newSF)["G"];
  setRealArgument("f12EnergyCorrection",f12EnergyCorrection);


}

