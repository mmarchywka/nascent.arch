#ifndef MJM_MATRIX_OPS_H__
#define MJM_MATRIX_OPS_H__


/***

A collection of methods for dealing with points in space such as
atoms in molecule.

*/

#include "mjm_globals.h"
#include "mjm_conflicts.h" //  #include <cblas.h>

// lapack sigs
extern "C"
{       
/*
void zheevr_(char* JOBZ, char* RANGE, char* UPLO, int * N, CoTy* A, int * LDA,
                double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
                double* W, CoTy* Z, int* LDZ, int* ISUPPZ, CoTy* WORK, int* LWORK,
                double* RWORK, int* LRWORK, int* IWORK, int* LIWORK, int* INFO);
*/

//void cgesv_(int * N, int * NRHS, CoTy * A, int* LDA, int * IPIV, CoTy * B, int * LDB, int* INFO);

// DGEBRD real double SVD
void dgesvd_(char * JOBU, char * JOBVT,int * M, int * N, double * A,
 int *LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK, 
int * INFO);


}; // lapacks
namespace mjm_misc_mat
{
// miscellaneous tings
// composite a square xform mat for repeated dirs
template <class Ts, class Td> void composite(Td & d, const Ts & s, const int dim, const int levels, const int skips=0)
{
int reps=0;
int dimf=dim;
// levels is "1" for 2x, with skips at 1 this makes 
// one iter through the d[i]j[] loop
const int sl=levels-skips;
while ( reps<levels) { dimf*=dim; ++reps; } 
for ( int i=0; i<dimf; ++i)
{
for ( int j=0; j<dimf; ++j)
{
d[i][j]=1;
int df=1;
for ( int k=0; k<=levels; ++k)
{
int il=(i/df)%dim;
int jl=(j/df)%dim;
const bool ident=(k>sl);
 d[i][j]*= ident?((il==jl)?1:0):s(il,jl);
df*=dim;
}

} //j 
} // i 

}


template <class Ts, class Td,class Tv> void mul(Td & d, const Ts & s, const Tv & v, const int col, const int rows)
{
for ( int j=0; j<col; ++j)
{
d[j]=0; 
for ( int i=0; i<rows; ++i) { d[j]+=s[i][j]*v[i]; }
}
}

}; // namespace mjm_misc_mat


//typedef mjm_generic_traits mjm_node_traits;
namespace mjm_functors
{
typedef double DefaultD;
typedef unsigned int IdxTy;
class average
{

typedef average Myt;
typedef DefaultD D;
enum { DIM=3};
// just dump in other crap too LOL.
D m[DIM];
D m_min[DIM],m_max[DIM];
IdxTy m_n;
 
public:
average() { 
for (unsigned int i=0; i<DIM; ++i) m[i]=0; 
for (unsigned int i=0; i<DIM; ++i) m_min[i]=1e100; 
for (unsigned int i=0; i<DIM; ++i) m_max[i]=-1e100; 
}

~average() {  } 

template < class Ts> Myt & operator()(const Ts & loc)
{
m[0]+=loc[0];
m[1]+=loc[1];
m[2]+=loc[2];
Maxx(m_max[0],loc[0]);
Minn(m_min[0],loc[0]);
Maxx(m_max[1],loc[1]);
Minn(m_min[1],loc[1]);
Maxx(m_max[2],loc[2]);
Minn(m_min[2],loc[2]);



++m_n;
return *this;
}

template < class Td > void solve(Td & result)
{
// should not throw on zero n
if (m_n==0) m_n=1; // counters should be initzlied, waste the div 
result[0]=m[0]/m_n;
result[1]=m[1]/m_n;
result[2]=m[2]/m_n;

}
private:
template <class Td> void Maxx(Td & d, const Td &s) { if ( s>d) d=s; } 
template <class Td> void Minn(Td & d, const Td &s) { if ( s<d) d=s; } 

}; // average
















class moments
{
typedef moments Myt;
typedef DefaultD D;
typedef D * MatrixTy;

       MatrixTy m;
const D m_x,m_y,m_z;
 
public:
//moments() { m= new D[9]; for (unsigned int i=0; i<9; ++i) m[i]=0; }
moments(const D x, const D y, const D z):m_x(x),m_y(y),m_z(z) 
{ m= new D[9]; for (unsigned int i=0; i<9; ++i) m[i]=0; }

~moments() { delete[] m; } 

template < class Ts> Myt & operator()(const Ts & loc)
{
// hopefully compiler can optimize the mults LOL. 
const D x=loc[0]-m_x;
const D y=loc[1]-m_y;
const D z=loc[2]-m_z;
const D xx=x*x;
const D yy=y*y;
const D zz=z*z;

m[0]+=yy+zz;
        m[1*3+1]+=xx+zz;
        m[2*3+2]+=xx+yy;
// this is dumb, just do it on "solve" ..... 
        m[0*3+1]-=x*y; //m[1][0]=m[0][1];
        m[0*3+2]-=x*z; //m[2][0]=m[0][2];
        m[1*3+2]-=y*z; //m[2][1]=m[1][2];

return *this;
}

// do svd to invert 

//void dgesvd_(char * JOBU, char * JOBVT,int * M, int * N, double * A,
// int *LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK, 
//int * INFO);
template < class Td > void solve(Td & result)
{
       m[1*3+0]= m[0*3+1]; // -=loc[0]*loc[1]; //m[1][0]=m[0][1];
       m[2*3+0]= m[0*3+2]; // -=loc[0]*loc[2]; //m[2][0]=m[0][2];
       m[2*3+1]= m[1*3+2]; // -=loc[1]*loc[2]; //m[2][1]=m[1][2];
int d=3;
int dw=d*5;
int info=0;
D S[d];
D U[d*d];
D VT[d*d];
D W[d*5];
char x='A';
std::cout<<MM_MARK<<" svd input "<<info<<" SVD:"<<CRLF; // S[0]<<" "<<S[1]<<" "<<S[2];
dgesvd_(&x,&x,&d,&d,m,&d,S,U,&d,VT,&d,W,&dw,&info);

std::cout<<MM_MARK<<" svd output "<<info<<" SVD:"<<S[0]<<" "<<S[1]<<" "<<S[2];
std::cout<<CRLF;
result[0]=S[0];
result[1]=S[1];
result[2]=S[2];
int dim1=3;
double alpha=1;
double beta=0;
D * A =VT;
D * B=U;
D C[9];
cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim1, dim1, dim1, alpha, A, dim1, B, dim1, beta, C, dim1);
	

}


}; //moments








}; // functors ns
#endif

