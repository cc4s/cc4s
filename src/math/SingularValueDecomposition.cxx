#include <math/SingularValueDecomposition.hpp>

#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <complex>
/*
#include <mkl.h>
#include <mkl_scalapack.h>
#include "mkl_lapacke.h"
#include <mkl_cblas.h>
*/
extern "C" {
  void pzgesvd_(
    const char *jobu, const char *jobvt,
    const int *m, const int *n,
    const complex *a, const int *ia, const int *ja, const int *desca,
    double *s, complex *u, const int *iu, const int *ju, const int *descu,
    complex *vt, const int *ivt, const int *jvt, const int *descvt,
    complex *work, const int *lwork, double *rwork, int *info
  );
};

using namespace cc4s;
using namespace CTF;


template <typename F>
SingularValueDecomposition<F>::SingularValueDecomposition(
  Matrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
Matrix<F> &SingularValueDecomposition<F>::get() {
  int n(inverse.lens[0]);
  int iA, jA;
  int descA[9];
  complex *a;
  double *sigma;
  complex *U, *VT;
  int iU, jU, iVT, jVT;
  int descU[9], descVT[9];
  complex *work;
  int workCount;
  double *realWork;
  int info;
  pzgesvd_(
    "V", "V", &n, &n, a, &iA, &jA, descA,
    sigma,
    U, &iU, &jU, descU,
    VT, &iVT, &jVT, descVT,
    work, &workCount, realWork,
    &info
  );
  return inverse;
}

// instantiate
template
SingularValueDecomposition<double>::SingularValueDecomposition(
  Matrix<double> const &matrix
);
template
Matrix<double> &SingularValueDecomposition<double>::get();

template
SingularValueDecomposition<complex>::SingularValueDecomposition(
  Matrix<complex> const &matrix
);
template
Matrix<complex> &SingularValueDecomposition<complex>::get();



template <typename F>
DrySingularValueDecomposition<F>::DrySingularValueDecomposition(
  DryMatrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
DryMatrix<F> &DrySingularValueDecomposition<F>::get() {
  return inverse;
}

// instantiate
template
DrySingularValueDecomposition<double>::DrySingularValueDecomposition(
  DryMatrix<double> const &matrix
);
template
DryMatrix<double> &DrySingularValueDecomposition<double>::get();

template
DrySingularValueDecomposition<complex>::DrySingularValueDecomposition(
  DryMatrix<complex> const &matrix
);
template
DryMatrix<complex> &DrySingularValueDecomposition<complex>::get();

/*
#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "Scalapack.h"
#include <mkl.h>
#include <mkl_scalapack.h>
#include "mkl_lapacke.h"
#include <mkl_cblas.h>

#define mat(matriz,coluna,i,j) (matriz[i*coluna+j])

#define p_of_i(i,bs,p) ( MKL_INT((i-1)/bs)%p)
#define l_of_i(i,bs,p) ( MKL_INT((i-1)/(p*bs)))
#define x_of_i(i,bs,p) (((i-1)%bs)+1)

#define   numroc_      NUMROC

using namespace std;

extern "C" 
{
    void Cblacs_pinfo(int* mypnum, int* nprocs);
    void Cblacs_get( MKL_INT context, MKL_INT request, MKL_INT* value);
    int  Cblacs_gridinit( MKL_INT* context, char * order, MKL_INT np_row, MKL_INT np_col);
    void Cblacs_gridinfo( MKL_INT context, MKL_INT*  np_row, MKL_INT* np_col, MKL_INT*  my_row,
    MKL_INT*  my_col);
    int  numroc_( MKL_INT *n, MKL_INT *nb, MKL_INT *iproc, MKL_INT *isrcproc, MKL_INT *nprocs);
    void Cblacs_gridexit(MKL_INT ictxt);
    void Cblacs_barrier(MKL_INT ictxt, char * order);
}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol);
int getIndex(MKL_INT row, MKL_INT col,MKL_INT NCOLS) {return row*NCOLS+col;}

CTEST_Scalapack::CTEST_Scalapack(void)
{
}

CTEST_Scalapack::~CTEST_Scalapack(void)
{
}

int CTEST_Scalapack::Scalapack(int argc, char ** argv) 
{

    int nprocs = 0;//MPI::COMM_WORLD.Get_size();
    int rank = 0;//MPI::COMM_WORLD.Get_rank();

    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::cout<<"Returned: "<<" ";
    std::cout << "Hello World! I am " << rank << " of " << nprocs <<
    std::endl;

    srand(1);
    MKL_INT myrow=0;
    MKL_INT mycol=0;
    MKL_INT ictxt=0;
    MKL_INT nprow=0,npcol=0;

    MKL_INT BLOCK_SIZE =2; //this gonna be tricky - should be 64, but cannot be larger than the original matrix

    MKL_INT locR=0, locC=0;
    MKL_INT block = BLOCK_SIZE;
    MKL_INT izero = 0;
    MKL_INT matrix_size = 9;
   
    MKL_INT myone = 1;
    
    MKL_INT nrhs = 1;
   
    MKL_INT info=0;
  
    int i=0,j=0;
    double mone=(-1.e0),pone=(1.e0);
    double AnormF=0.e0, XnormF=0.e0, RnormF=0.e0, BnormF=0.e0, residF=0.e0,eps=0.e0;

    find_nps(nprocs,nprow,npcol);

    Cblacs_pinfo( &rank, &nprocs ) ;
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    
    locR = numroc_(&matrix_size, &block, &myrow, &izero, &nprow);
    locC = numroc_(&matrix_size, &block, &mycol, &izero, &npcol);

   
    ////GLOBAL
    double * A = new double[matrix_size*matrix_size]();
    double * B = new double[matrix_size]();
    double * Acpy = new double[matrix_size*matrix_size]();
    double * Bcpy = new double[matrix_size]();
    
    //LOCAL
    double * local_know_vector = new double[locR]();
    double * local_matrix = new double[locR*locC]();
    
    MKL_INT* ipiv = new MKL_INT [locC*locR*block+1000000]();

    
    B[2] = 1;
    B[3] = 0;
    B[4] = 0;
    B[5] = 0;
    
    
    
    A[0] = 19;
    A[1] = 3;
    A[2] = 1;
    A[3] = 12;
    A[4] = 1;
    A[5] = 16;
    A[6] = 1;
    A[7] = 3;
    A[8] = 11;
    
    A[9] = -19;
    A[10] = 3;
    A[11] = 1;
    A[12] = 12;
    A[13] = 1;
    A[14] = 16;
    A[15] = 1;
    A[16] = 3;
    A[17] = 11;
    
    A[18] = -19;
    A[19] = -3;
    A[20] = 1;
    A[21] = 12;
    A[22] = 1;
    A[23] = 16;
    A[24] = 1;
    A[25] = 3;
    A[26] = 11;
    
    A[27] = -19;
    A[28] = -3;
    A[29] = -1;
    A[30] = 12;
    A[31] = 1;
    A[32] = 16;
    A[33] = 1;
    A[34] = 3;
    A[35] = 11;
    
    A[36] = -19;
    A[37] = -3;
    A[38] = -1;
    A[39] = -12;
    A[40] = 1;
    A[41] = 16;
    A[42] = 1;
    A[43] = 3;
    A[44] = 11;
    
    A[45] = -19;
    A[46] = -3;
    A[47] = -1;
    A[48] = -12;
    A[49] = -1;
    A[50] = 16;
    A[51] = 1;
    A[52] = 3;
    A[53] = 11;
    
    A[54] = -19;
    A[55] = -3;
    A[56] = -1;
    A[57] = -12;
    A[58] = -1;
    A[59] = -16;
    A[60] = 1;
    A[61] = 3;
    A[62] = 11;
    
    A[63] = -19;
    A[64] = -3;
    A[65] = -1;
    A[66] = -12;
    A[67] = -1;
    A[68] = -16;
    A[69] = -1;
    A[70] = 3;
    A[71] = 11;
    
    A[72] = -19;
    A[73] = -3;
    A[74] = -1;
    A[75] = -12;
    A[76] = -1;
    A[77] = -16;
    A[78] = -1;
    A[79] = -3;
    A[80] = 11;

    MKL_INT* descA  = new MKL_INT[9]();
    MKL_INT* descB  = new MKL_INT[9]();
   
    descA[0] = 1; // descriptor type
    descA[1] = ictxt; // blacs context
    descA[2] = matrix_size; // global number of rows
    descA[3] = matrix_size; // global number of columns
    descA[4] = block; // row block size
    descA[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descA[6] = 0; // initial process row(DEFINED 0)
    descA[7] = 0; // initial process column (DEFINED 0)
    descA[8] = locR; // leading dimension of local array

    descB[0] = 1; // descriptor type
    descB[1] = ictxt; // blacs context
    descB[2] = matrix_size; // global number of rows
    descB[3] = 1; // global number of columns
    descB[4] = block; // row block size
    descB[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descB[6] = 0; // initial process row(DEFINED 0)
    descB[7] = 0; // initial process column (DEFINED 0)
    descB[8] = locR; // leading dimension of local array

    int il=0, jl=0;
    for(i=1; i< matrix_size+1; i++) 
    {
       for(j=1; j< matrix_size+1; j++) 
       {
    
        int pi = p_of_i(i,block,nprow);
        
        int li = l_of_i(i,block,nprow);

        int xi = x_of_i(i,block,nprow);
        //printf("i = %d, j = %d, pi = %d, li = %d\n",i,j,pi,li);;fflush(stdout);
        int pj = p_of_i(j,block,npcol);
        
        int lj = l_of_i(j,block,npcol);
        
        int xj = x_of_i(j,block,npcol);
        //printf("i = %d, j = %d, pj = %d, lj = %d, xj = %d\n",i,j,pj,lj,xj);;fflush(stdout);

        if( (pi == myrow) && (pj == mycol)) 
        {
            il = li*block+xi;
            jl = lj*block+xj;
            local_matrix[getIndex(il-1, jl-1, locC)] = A[getIndex(i-1,j-1,matrix_size)];
        }
    
        if(  (pi == myrow) &&(mycol==0)  )
        {
            local_know_vector[il-1] = B[i-1];
        }

       }
    
    }
      
    ////STARTING PDGESV
    pdgesv_(&matrix_size, &nrhs, local_matrix, &myone, &myone, descA, ipiv, local_know_vector, &myone, &myone, descB, &info);
    
    if(rank==0)
      {
        if(info != 0) cout <<"PDGESV problem! Info "<<info<<endl;
      }
    
    
    for(i=0; i< locR; i++)
    {
      cout<<"**\n"<<"rank "<<rank<<"  answer: "<<local_know_vector[i]<<endl;
    }

    if(NULL!=descA)                        {delete [] descA; descA=NULL;} 
    if(NULL!=descB)                        {delete [] descB; descB=NULL;} 
    if(NULL!=local_know_vector)            {delete [] local_know_vector; local_know_vector=NULL;} 
    if(NULL!=local_matrix)                {delete [] local_matrix; local_matrix=NULL;} 
    if(NULL!=Acpy)                        {delete [] Acpy; Acpy=NULL;} 
    if(NULL!=Bcpy)                        {delete [] Bcpy; Bcpy=NULL;} 
    if(NULL!=A)                            {delete [] A; A=NULL;} 
    if(NULL!=B)                            {delete [] B; B=NULL;} 
    

    Cblacs_gridexit(ictxt);

    return 0;

}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol) 
{

MKL_INT min_nprow=100000;
MKL_INT min_npcol=100000;

nprow = np;
npcol = np;

while(1) {

   npcol--;
  if(np%2==0   ) {
  if(npcol ==1){
   nprow --;
   npcol = nprow;
  }
  }else {
  if(npcol ==0){
   nprow --;
   npcol = nprow;
  }

  }

  if(nprow*npcol == np) {
    min_npcol = npcol;
    if(nprow < min_nprow)    min_nprow = nprow;
  }

    if(nprow ==1 ) break;

}

nprow = min_nprow;
npcol = min_npcol;

}
RSS
Top
10 posts / 0 new
Last post
For more complete information about compiler optimizations, see our Optimization Notice.
[Status Points: 6,610 (Total Points: 39,385)]

    Log in to post comments

Ying H. (Intel)
Tue, 12/30/2014 - 21:27

Hi 

I’ve corrected the example to properly define process grid and init MPI.

Now it works fine for any number of MPI processes. Attach here. 

and the article :  https://software.intel.com/en-us/articles/using-cluster-mkl-pblasscalapack-fortran-routine-in-your-c-program/  had pdgemv sample code https://software.intel.com/sites/default/files/article/165948/pdgemv.c

#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
//#include "Scalapack.h"
#include <mkl.h>
#include <mkl_scalapack.h>
#include "mkl_lapacke.h"
#include <mkl_cblas.h>

#define mat(matriz,coluna,i,j) (matriz[i*coluna+j])

#define p_of_i(i,bs,p) ( MKL_INT((i-1)/bs)%p)
#define l_of_i(i,bs,p) ( MKL_INT((i-1)/(p*bs)))
#define x_of_i(i,bs,p) (((i-1)%bs)+1)

//#define   numroc_      NUMROC

using namespace std;

extern "C"
{
    void Cblacs_pinfo(int* mypnum, int* nprocs);
    void Cblacs_get( MKL_INT context, MKL_INT request, MKL_INT* value);
    int  Cblacs_gridinit( MKL_INT* context, char * order, MKL_INT np_row, MKL_INT np_col);
    void Cblacs_gridinfo( MKL_INT context, MKL_INT*  np_row, MKL_INT* np_col, MKL_INT*  my_row,
    MKL_INT*  my_col);
    int  numroc_( MKL_INT *n, MKL_INT *nb, MKL_INT *iproc, MKL_INT *isrcproc, MKL_INT *nprocs);
    void Cblacs_gridexit(MKL_INT ictxt);
    void Cblacs_barrier(MKL_INT ictxt, char * order);
}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol);
int getIndex(MKL_INT row, MKL_INT col,MKL_INT NCOLS) {return row*NCOLS+col;}

int main(int argc, char ** argv)
{

    int nprocs = 0;//MPI::COMM_WORLD.Get_size();
    int rank = 0;//MPI::COMM_WORLD.Get_rank();

    MPI_Init(&argc,&argv);   
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::cout<<"Returned: "<<" ";
    std::cout << "Hello World! I am " << rank << " of " << nprocs <<
    std::endl;

    srand(1);
    MKL_INT myrow=0;
    MKL_INT mycol=0;
    MKL_INT ictxt=0;
    MKL_INT nprow=0,npcol=0;

    MKL_INT BLOCK_SIZE =2; //this gonna be tricky - should be 64, but cannot be larger than the original matrix

    MKL_INT locR=0, locC=0;
    MKL_INT block = BLOCK_SIZE;
    MKL_INT izero = 0;
    MKL_INT matrix_size = 9;
  
    MKL_INT myone = 1;
   
    MKL_INT nrhs = 1;
  
    MKL_INT info=0;
 
    int i=0,j=0;
    double mone=(-1.e0),pone=(1.e0);
    double AnormF=0.e0, XnormF=0.e0, RnormF=0.e0, BnormF=0.e0, residF=0.e0,eps=0.e0;

    find_nps(nprocs,nprow,npcol);

    Cblacs_pinfo( &rank, &nprocs ) ;
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
   
    locR = numroc_(&matrix_size, &block, &myrow, &izero, &nprow);
    locC = numroc_(&matrix_size, &block, &mycol, &izero, &npcol);

  
    ////GLOBAL
    double * A = new double[matrix_size*matrix_size]();
    double * B = new double[matrix_size]();
    double * Acpy = new double[matrix_size*matrix_size]();
    double * Bcpy = new double[matrix_size]();
   
    //LOCAL
    double * local_know_vector = new double[locR]();
    double * local_matrix = new double[locR*locC]();
   
    MKL_INT* ipiv = new MKL_INT [locC*locR*block+1000000]();

 B[0] = 0;
    B[1] = 0;
    B[2] = 1;
    B[3] = 0;
    B[4] = 0;
    B[5] = 0;
 B[6] = 0;
    B[7] = 0;
   
   
   
    A[0] = 19;
    A[1] = 3;
    A[2] = 1;
    A[3] = 12;
    A[4] = 1;
    A[5] = 16;
    A[6] = 1;
    A[7] = 3;
    A[8] = 11;
   
    A[9] = -19;
    A[10] = 3;
    A[11] = 1;
    A[12] = 12;
    A[13] = 1;
    A[14] = 16;
    A[15] = 1;
    A[16] = 3;
    A[17] = 11;
   
    A[18] = -19;
    A[19] = -3;
    A[20] = 1;
    A[21] = 12;
    A[22] = 1;
    A[23] = 16;
    A[24] = 1;
    A[25] = 3;
    A[26] = 11;
   
    A[27] = -19;
    A[28] = -3;
    A[29] = -1;
    A[30] = 12;
    A[31] = 1;
    A[32] = 16;
    A[33] = 1;
    A[34] = 3;
    A[35] = 11;
   
    A[36] = -19;
    A[37] = -3;
    A[38] = -1;
    A[39] = -12;
    A[40] = 1;
    A[41] = 16;
    A[42] = 1;
    A[43] = 3;
    A[44] = 11;
   
    A[45] = -19;
    A[46] = -3;
    A[47] = -1;
    A[48] = -12;
    A[49] = -1;
    A[50] = 16;
    A[51] = 1;
    A[52] = 3;
    A[53] = 11;
   
    A[54] = -19;
    A[55] = -3;
    A[56] = -1;
    A[57] = -12;
    A[58] = -1;
    A[59] = -16;
    A[60] = 1;
    A[61] = 3;
    A[62] = 11;
   
    A[63] = -19;
    A[64] = -3;
    A[65] = -1;
    A[66] = -12;
    A[67] = -1;
    A[68] = -16;
    A[69] = -1;
    A[70] = 3;
    A[71] = 11;
   
    A[72] = -19;
    A[73] = -3;
    A[74] = -1;
    A[75] = -12;
    A[76] = -1;
    A[77] = -16;
    A[78] = -1;
    A[79] = -3;
    A[80] = 11;

    MKL_INT* descA  = new MKL_INT[9]();
    MKL_INT* descB  = new MKL_INT[9]();
  
    descA[0] = 1; // descriptor type
    descA[1] = ictxt; // blacs context
    descA[2] = matrix_size; // global number of rows
    descA[3] = matrix_size; // global number of columns
    descA[4] = block; // row block size
    descA[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descA[6] = 0; // initial process row(DEFINED 0)
    descA[7] = 0; // initial process column (DEFINED 0)
    descA[8] = locR; // leading dimension of local array

    descB[0] = 1; // descriptor type
    descB[1] = ictxt; // blacs context
    descB[2] = matrix_size; // global number of rows
    descB[3] = 1; // global number of columns
    descB[4] = block; // row block size
    descB[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descB[6] = 0; // initial process row(DEFINED 0)
    descB[7] = 0; // initial process column (DEFINED 0)
    descB[8] = locR; // leading dimension of local array

    int il=0, jl=0;
    for(i=1; i< matrix_size+1; i++)
    {
       for(j=1; j< matrix_size+1; j++)
       {
   
        int pi = p_of_i(i,block,nprow);
       
        int li = l_of_i(i,block,nprow);

        int xi = x_of_i(i,block,nprow);
        //printf("i = %d, j = %d, pi = %d, li = %d\n",i,j,pi,li);;fflush(stdout);
        int pj = p_of_i(j,block,npcol);
       
        int lj = l_of_i(j,block,npcol);
       
        int xj = x_of_i(j,block,npcol);
        //printf("i = %d, j = %d, pj = %d, lj = %d, xj = %d\n",i,j,pj,lj,xj);;fflush(stdout);

        if( (pi == myrow) && (pj == mycol))
        {
            il = li*block+xi;
            jl = lj*block+xj;
            local_matrix[getIndex(il-1, jl-1, locC)] = A[getIndex(i-1,j-1,matrix_size)];
        }
   
        if(  (pi == myrow) &&(mycol==0)  )
        {
            local_know_vector[il-1] = B[i-1];
        }

       }
   
    }
     
//  above initialization code have some issue. in our sample, we use   
for(i=0;i<locR;++i) for(j=0;j<locC;++j) {
  int gi = i%block + block*myrow + (i/block)*block*nprow;
  int gj = j%block + block*mycol + (j/block)*block*npcol;
  local_matrix[i + locR*j] = A[gi*matrix_size + gj]; // note: col-major <- row-major
  local_know_vector[i] = B[gi];
}

    ////STARTING PDGESV
    pdgesv_(&matrix_size, &nrhs, local_matrix, &myone, &myone, descA, ipiv, local_know_vector, &myone, &myone, descB, &info);
   
    if(rank==0)
      {
        if(info != 0) cout <<"PDGESV problem! Info "<<info<<endl;
      }
   
    for (i = 0; i < locR; i++)
{
  int gi = i%block + block*myrow + (i/block)*block*nprow;
  if (mycol == 0) printf("res[ %i ]=%lg\n",gi,local_know_vector[i]);
}

  //  for(i=0; i< locR; i++)
  //  {
  //    cout<<"**\n"<<"rank "<<rank<<"  answer: "<<local_know_vector[i]<<endl;
  //  }

    if(NULL!=descA)                        {delete [] descA; descA=NULL;}
    if(NULL!=descB)                        {delete [] descB; descB=NULL;}
    if(NULL!=local_know_vector)            {delete [] local_know_vector; local_know_vector=NULL;}
    if(NULL!=local_matrix)                {delete [] local_matrix; local_matrix=NULL;}
    if(NULL!=Acpy)                        {delete [] Acpy; Acpy=NULL;}
    if(NULL!=Bcpy)                        {delete [] Bcpy; Bcpy=NULL;}
    if(NULL!=A)                            {delete [] A; A=NULL;}
    if(NULL!=B)                            {delete [] B; B=NULL;}
   

    Cblacs_gridexit(ictxt);

    return 0;

}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol)
{
#if 1
    nprow = (int)sqrt( np );
    npcol = np / nprow;
    return;

#else

MKL_INT min_nprow=100000;
MKL_INT min_npcol=100000;

nprow = np;
npcol = np;

while(1) {

   npcol--;
  if(np%2==0   ) {
  if(npcol ==1){
   nprow --;
   npcol = nprow;
  }
  }else {
  if(npcol ==0){
   nprow --;
   npcol = nprow;
  }

  }

  if(nprow*npcol == np) {
    min_npcol = npcol;
    if(nprow < min_nprow)    min_nprow = nprow;
  }

    if(nprow ==1 ) break;

}

nprow = min_nprow;
npcol = min_npcol;
#endif
}

 
Top

    Log in to post comments

Yang Y.
Wed, 01/13/2016 - 21:59

Hi Ying,

I ran your code on NERSC Cori supercomputer, but I got the error information:

Returned:  Hello World! I am 0 of 2

Returned:  Hello World! I am 1 of 2

Rank 0 [Wed Jan 13 12:23:10 2016] [c0-0c0s8n3] Fatal error in MPI_Send: Invalid tag, error stack:

MPI_Send(186): MPI_Send(buf=0x139d640, count=2, MPI_INT, dest=1, tag=5000000, comm=0x84000001) failed

MPI_Send(111): Invalid tag, value is 5000000

Rank 1 [Wed Jan 13 12:23:10 2016] [c0-0c0s9n0] Fatal error in MPI_Recv: Invalid tag, error stack:

MPI_Recv(199): MPI_Recv(buf=0x139b5c0, count=2, MPI_INT, src=0, tag=5000000, comm=0x84000001, status=0x7fffffff5618) failed

MPI_Recv(118): Invalid tag, value is 5000000

srun: error: nid00035: task 0: Aborted

srun: Terminating job step 929053.0

srun: Job step aborted: Waiting up to 32 seconds for job step to finish.

srun: error: nid00036: task 1: Aborted

Could you help me?
Top
[Status Points: 6,610 (Total Points: 39,385)]

    Log in to post comments

Ying H. (Intel)
Wed, 01/13/2016 - 23:55

Hi Yang, 

Could you also show your compile command and run command, MPI version, MKL version etc information? 

There was some MPI sample in MKL install directory, could you please try them first and see if mkl work fine on the supercomputer? 

Best Regards,

Ying 
Top

    Log in to post comments

Yang Y.
Thu, 01/14/2016 - 16:53

Hi Ying,

I'm using the NERSC cori supercomputer, all the information is at:

http://www.nersc.gov/users/computational-systems/cori/

The compile command is:

source /opt/intel/bin/compilervars.sh intel64

source /opt/intel/impi/5.0.2.044/bin64/mpivars.sh

CC -o execute pdgesv.cpp -mkl:cluster

The run command is sbatch test.sl

where test.sl is

#!/bin/bash -l

#SBATCH --partition debug

#SBATCH --nodes 2

#SBATCH --time=00:03:00

cd $SLURM_SUBMIT_DIR

srun -n 2 ./execute;
Top
[Status Points: 440 (Total Points: 6,005)]

    Log in to post comments

Dmitry Baksheev (Intel)
Thu, 01/14/2016 - 22:32

Hi Yang,

You are probably using Cray MPI that does not support tags having large value. MKL used to use such tags internally. Now the issue is fixed, and if you use a recent MKL version (e.g. MKL 11.3.1) the test should pass.

Alternatively, it should work with Intel MPI.

Thanks
Dima
 
Top

    Log in to post comments

Yang Y.
Sat, 02/20/2016 - 10:41

If I change the number of processes, the solution of linear equation will be different. Also, the solutions of this code is different from the Matlab solution. Let us solve the linear equation Ax = b

A is 

    19     3     1    12     1    16     1     3    11
   -19     3     1    12     1    16     1     3    11
   -19    -3     1    12     1    16     1     3    11
   -19    -3    -1    12     1    16     1     3    11
   -19    -3    -1   -12     1    16     1     3    11
   -19    -3    -1   -12    -1    16     1     3    11
   -19    -3    -1   -12    -1   -16     1     3    11
   -19    -3    -1   -12    -1   -16    -1     3    11
   -19    -3    -1   -12    -1   -16    -1    -3    11

b is

     0
     0
     1
     0
     0
     0
     0
     0
     0

If we use Matlab, the solution is:

         0
   -0.1667
    0.5000
         0
         0
         0
         0
         0
         0

If we use this code for one process, the solution is:

answer[0]: 0.000000

answer[1]: 0.000000

answer[2]: 0.500000

answer[3]: -0.500000

answer[4]: -0.000000

answer[5]: -0.000000

answer[6]: 0.000000

answer[7]: -0.000000

answer[8]: -0.000000

If we use this code for two processes, the solution is:

answer[0]: 0.003434

answer[1]: 0.016138

answer[2]: 0.019397

answer[3]: 0.033449

answer[4]: -0.025858

answer[5]: 0.013383

answer[6]: 0.022018

answer[7]: -0.060255

answer[8]: 0.025895

If we use this code for three processes, the solution is:

answer[0]: -0.023551

answer[1]: 0.017566

answer[2]: 0.004775

answer[3]: 0.004377

answer[4]: 0.023704

answer[5]: 0.027967

answer[6]: -0.023534

answer[7]: 0.017509

answer[8]: 0.000398

Can you tell me where went wrong? How can I get the correct solution by this code?
Top
[Status Points: 440 (Total Points: 6,005)]

    Log in to post comments

Dmitry Baksheev (Intel)
Wed, 02/24/2016 - 02:38

The matrix should be distributed correctly. For example, consider distribution of rows. Let N be global number of rows, LOCR the local number of rows as obtained with a call to numroc(), B the block size, p my processor row, P the number of processor rows. Then for locn = 0...LOCR-1 the global index n=0...N-1  is computed thus:

n = locn % B // offset within current block (block locn/B)
   + B*p   // offset of the current block within current cycle of block-cyclic distribution
   + P*B*(locn/B)   // offset of the current cycle on the grid of P processors
   ;
 

Once you know (i,j) for (loci,locj) you can initialize local_matrix[loci + locR*locj] = A[i + matrix_size*j]. Please note that ScaLAPACK supports only column-major ordering of matrices. That is, local matrix element (i,j) should be located at offset i+locR*j, not at i*locC + j,

I hope this will help you.

Thanks
Dima

 

 

 
Top
[Status Points: 6,610 (Total Points: 39,385)]

    Log in to post comments

Ying H. (Intel)
Wed, 02/24/2016 - 22:59

Hi Dima,

Thank you much for the correction. 

The proper code for initialization of local matrix should be inserted before the call of pdgesv:  ( I will reedit the above of the code)

for(i=0;i<locR;++i) for(j=0;j<locC;++j) {

  int gi = i%block + block*myrow + (i/block)*block*nprow;

  int gj = j%block + block*mycol + (j/block)*block*npcol;

  local_matrix[i + locR*j] = A[gi*matrix_size + gj]; // note: col-major <- row-major

  local_know_vector[i] = B[gi];

}

 

The solution can be gathered to some rank, or it can be just printed:

 

for (i = 0; i < locR; i++)

{

  int gi = i%block + block*myrow + (i/block)*block*nprow;

  if (mycol == 0) printf("res[ %i ]=%lg\n",gi,local_know_vector[i]);

}

 

 

Here are a few runs:

 

$ mpirun -np 1 --prepend-rank ./a.out | sort -k2

[0] res[ 0 ]=0

[0] res[ 1 ]=-0.166667

[0] res[ 2 ]=0.5

[0] res[ 3 ]=0

[0] res[ 4 ]=0

[0] res[ 5 ]=0

[0] res[ 6 ]=0

[0] res[ 7 ]=0

[0] res[ 8 ]=0

 

 

$ mpirun -np 4 --prepend-rank ./a.out | sort -k2

[0] res[ 0 ]=0

[0] res[ 1 ]=-0.166667

[2] res[ 2 ]=0.5

[2] res[ 3 ]=0

[0] res[ 4 ]=0

[0] res[ 5 ]=0

[2] res[ 6 ]=0

[2] res[ 7 ]=0

[0] res[ 8 ]=0

 

Thanks

Ying
Top

    Log in to post comments

Yang Y.
Wed, 03/16/2016 - 19:41

Thank you, Dima and Ying, it works very well :-)
Top
Back to original post
Leave a Comment
Please sign in to add a comment. Not a member? Join today

    Support Terms of Use *Trademarks Privacy Cookies Publications 

Look for us on:

    Facebook
    Twitter
    Google+
    LinkedIn
    YouTube

    English
*/
